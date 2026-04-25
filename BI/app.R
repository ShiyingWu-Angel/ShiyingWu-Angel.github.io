# SET UP
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
library(readr)
library(teal.modules.clinical)
library("lmodel2")
library("mixtools")
library("rsconnect")
data_root <- file.path("resource", "data")
study_ids <- c("MOCK_001", "MOCK_002")

#HELPERS
read_csv_safe <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

parse_full_date_time <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("T", " ", x)
  
  out <- suppressWarnings(
    as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  )
  
  idx <- is.na(out)
  out[idx] <- suppressWarnings(
    as.POSIXct(x[idx], format = "%Y-%m-%d %H:%M", tz = "UTC")
  )
  
  idx <- is.na(out)
  out[idx] <- suppressWarnings(
    as.POSIXct(paste0(substr(x[idx], 1, 10), " 00:00:00"),
               format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  )
  
  out
}

parse_full_date <- function(x) {
  x <- as.character(x)
  as.Date(substr(x, 1, 10))
}

get_labels <- function(dataset) {
  names(dataset) |>
    purrr::set_names() |>
    purrr::map_chr(function(name) {
      label <- attributes(dataset[[name]])$label
      if (is.null(label)) "No label" else label
    })
}

set_labels <- function(dataset, labels) {
  labels <- labels[names(labels) %in% names(dataset)]
  dataset |>
    dplyr::mutate(dplyr::across(
      names(labels),
      ~ structure(.x, label = labels[[dplyr::cur_column()]])
    ))
}

sort_visits <- function(df) {
  df |>
    dplyr::select(dplyr::all_of(c("VISIT", "VISITNUM"))) |>
    dplyr::distinct() |>
    dplyr::arrange(VISITNUM) |>
    dplyr::pull(VISIT)
}

smallest_ref_range <- function(df, low, high) {
  lbstnrlo_max_tmp <- c()
  lbstnrhi_min_tmp <- c()
  var <- unique(df$LBTESTCD)
  
  for (v in var) {
    llo <- as.numeric(unlist(unique(df[df$LBTESTCD == v, low])))
    lhi <- as.numeric(unlist(unique(df[df$LBTESTCD == v, high])))
    
    lbstnrlo_max_tmp <- c(lbstnrlo_max_tmp, if (any(!is.na(llo))) max(llo, na.rm = TRUE) else NA)
    lbstnrhi_min_tmp <- c(lbstnrhi_min_tmp, if (any(!is.na(lhi))) min(lhi, na.rm = TRUE) else NA)
  }
  
  out <- data.frame(
    LBTESTCD = var,
    lbstnrlo_max = lbstnrlo_max_tmp,
    lbstnrhi_min = lbstnrhi_min_tmp
  )
  
  attributes(out$lbstnrlo_max)$label <- "Maximum Reference Range Lower Limit"
  attributes(out$lbstnrhi_min)$label <- "Minimum Reference Range Upper Limit"
  out
}

get_arm_defaults <- function(dm) {
  arm <- dplyr::coalesce(dm$ACTARM, dm$ARM)
  unique(arm[!is.na(arm) & arm != ""])
}

#DERIVED DATA FUNCTIONS
getRS <- function(sdtm_set) {
  rs_std_cols <- c("VISIT", "SP_RECIST_OVRLRESP", "SP_RECIST_TRGRESP", "SP_RECIST_NTRGRESP")
  
  choose_latest_rs <- function(df) {
    df %>%
      dplyr::group_by(STUDYID, USUBJID, VISITNUM, VISIT, RSTESTCD) %>%
      dplyr::arrange(RSDTC, .by_group = TRUE) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::ungroup()
  }
  
  rs_r <- sdtm_set$rs %>%
    dplyr::mutate(
      RSCAT = stringr::str_replace(RSCAT, stringr::fixed("RECIST1.1"), "RECIST 1.1"),
      RSSTAT = dplyr::coalesce(as.character(RSSTAT), "")
    ) %>%
    dplyr::filter(RSCAT == "RECIST 1.1", RSSTAT != "NOT DONE") %>%
    choose_latest_rs() %>%
    tidyr::pivot_wider(
      names_from = RSTESTCD,
      values_from = RSSTRESC,
      names_prefix = "SP_RECIST_",
      id_cols = c(STUDYID, USUBJID, VISITNUM, VISIT)
    )
  
  rs_i <- sdtm_set$rs %>%
    dplyr::mutate(RSSTAT = dplyr::coalesce(as.character(RSSTAT), "")) %>%
    dplyr::filter(RSCAT == "iRECIST", RSSTAT != "NOT DONE") %>%
    choose_latest_rs() %>%
    tidyr::pivot_wider(
      names_from = RSTESTCD,
      values_from = RSSTRESC,
      names_prefix = "SP_iRECIST_",
      id_cols = c(STUDYID, USUBJID, VISITNUM, VISIT)
    )
  
  rs_r %>%
    dplyr::left_join(rs_i, by = c("STUDYID", "USUBJID", "VISITNUM", "VISIT")) %>%
    dplyr::select(dplyr::any_of(rs_std_cols), dplyr::everything())
}

getBORs <- function(rs, confirmed = TRUE) {
  bestResponse <- function(time_point_responses, confirmed) {
    time_point_responses <- time_point_responses[!is.na(time_point_responses) & time_point_responses != ""]
    if (length(time_point_responses) == 0) return("NE")
    
    if (any(time_point_responses == "PD")) {
      time_point_responses <- time_point_responses[seq_len(grep("PD", time_point_responses)[1])]
    }
    
    runs <- rle(time_point_responses)
    best_response <- if (any(runs$values %in% c("CR", "PR", "SD"))) "SD" else "PD"
    
    min_run_length <- ifelse(confirmed, 2, 1)
    confirmed_CR <- any(runs$values == "CR" & runs$lengths >= min_run_length)
    confirmed_PR <- any(runs$values == "PR" & runs$lengths >= min_run_length)
    
    runs_PRCR <- rle(gsub("CR", "PR", time_point_responses))
    confirmed_PR_CR <- any(runs_PRCR$values == "PR" & runs_PRCR$lengths >= min_run_length)
    
    if (confirmed_CR) {
      best_response <- "CR"
    } else if (confirmed_PR || confirmed_PR_CR) {
      best_response <- "PR"
    }
    
    best_response
  }
  
  subj_ids <- unique(rs$USUBJID)
  
  bors <- sapply(subj_ids, function(subj_id) {
    time_point_resp <- unlist(
      dplyr::filter(rs, USUBJID == subj_id) %>%
        dplyr::select(dplyr::contains("OVRLRESP"))
    )
    bestResponse(time_point_resp, confirmed)
  })
  
  tibble::tibble(USUBJID = names(bors), BOR = bors)
}

getRespList <- function(sdtm_set) {
  dm <- sdtm_set$dm
  fa <- sdtm_set$fa
  cm <- sdtm_set$cm
  ds <- sdtm_set$ds
  tr <- sdtm_set$tr
  rs <- getRS(sdtm_set)
  rs2 <- sdtm_set$rs
  
  rt_set <- dm %>%
    dplyr::select(USUBJID, RFXSTDTC, ACTARMCD) %>%
    dplyr::filter(!stringr::str_detect(ACTARMCD, stringr::fixed("SCFAIL")) &
                    !stringr::str_detect(ACTARMCD, stringr::fixed("SCRNFAIL"))) %>%
    dplyr::left_join(
      rs %>%
        dplyr::select(USUBJID, SP_RECIST_OVRLRESP) %>%
        dplyr::filter(!is.na(SP_RECIST_OVRLRESP)) %>%
        getBORs(.),
      by = "USUBJID"
    )
  
  if (nrow(fa %>% dplyr::filter(FATESTCD == "PDL1DCPS")) != 0) {
    rt_set <- rt_set %>%
      dplyr::left_join(
        fa %>%
          dplyr::filter(FATESTCD == "PDL1DCPS") %>%
          dplyr::select(USUBJID, FASTRESC) %>%
          dplyr::rename(PDL1DCPS = FASTRESC),
        by = "USUBJID"
      )
  } else {
    rt_set <- rt_set %>% dplyr::mutate(PDL1DCPS = "NA")
  }
  
  if (nrow(cm %>% dplyr::filter(CMCAT == "SYSTEMIC ANTI-CANCER THERAPY", CMTRT == "CETUXIMAB")) != 0) {
    rt_set <- rt_set %>%
      dplyr::left_join(
        cm %>%
          dplyr::filter(CMCAT == "SYSTEMIC ANTI-CANCER THERAPY", CMTRT == "CETUXIMAB") %>%
          dplyr::mutate(CETUXIMAB_PRETRT = "Yes") %>%
          dplyr::distinct(USUBJID, CETUXIMAB_PRETRT) %>%
          dplyr::select(USUBJID, CETUXIMAB_PRETRT),
        by = "USUBJID"
      ) %>%
      dplyr::mutate(CETUXIMAB_PRETRT = tidyr::replace_na(CETUXIMAB_PRETRT, "No"))
  } else {
    rt_set <- rt_set %>% dplyr::mutate(CETUXIMAB_PRETRT = "NA")
  }
  
  rt_set <- rt_set %>%
    dplyr::left_join(
      ds %>%
        dplyr::filter(grepl("COMPLETION", DSSCAT) | grepl("WITHDRAWAL", DSTERM)) %>%
        dplyr::mutate(ONGOING = "No") %>%
        dplyr::select(USUBJID, ONGOING),
      by = "USUBJID"
    )
  
  completed <- ds %>%
    dplyr::filter(grepl("TRIAL COMPLETION", DSSCAT)) %>%
    dplyr::select(USUBJID) %>%
    dplyr::mutate(ONGOING2 = "No")
  
  rt_set <- rt_set %>%
    dplyr::left_join(completed, by = "USUBJID") %>%
    dplyr::mutate(
      ONGOING = ifelse(!is.na(ONGOING2), ONGOING2, ONGOING),
      ONGOING = ifelse(is.na(ONGOING), "Yes", ONGOING)
    ) %>%
    dplyr::select(-ONGOING2)
  
  if (("VIRLSTAT" %in% fa$FATESTCD) & ("FATSTDTL" %in% names(fa))) {
    rt_set <- rt_set %>%
      dplyr::left_join(
        fa %>%
          dplyr::filter(FATSTDTL == "HNSCC", FATESTCD == "VIRLSTAT") %>%
          dplyr::select(USUBJID, FAORRES) %>%
          dplyr::rename(HPV = FAORRES),
        by = "USUBJID"
      )
  } else {
    rt_set <- rt_set %>% dplyr::mutate(HPV = "NA")
  }
  
  rt_set <- rt_set %>%
    dplyr::left_join(
      tr %>%
        dplyr::filter(TRBLFL == "Y", TRTESTCD == "SUMDIAM") %>%
        dplyr::select(USUBJID, TRSTRESN) %>%
        dplyr::rename(Sum_Diam_BL = TRSTRESN),
      by = "USUBJID"
    ) %>%
    dplyr::left_join(
      fa %>%
        dplyr::filter(FATESTCD == "REGPRSIT") %>%
        dplyr::select(USUBJID, FASTRESC) %>%
        dplyr::rename(Prim_Tum_Loc = FASTRESC),
      by = "USUBJID"
    ) %>%
    dplyr::left_join(
      rs2 %>%
        dplyr::group_by(USUBJID) %>%
        dplyr::filter(!is.na(RSSTRESC), RSSTRESC != "") %>%
        dplyr::mutate(recistdate = as.Date(RSDTC)) %>%
        dplyr::summarise(Last_Recist = max(recistdate), .groups = "drop"),
      by = "USUBJID"
    ) %>%
    dplyr::left_join(
      fa %>%
        dplyr::group_by(USUBJID) %>%
        dplyr::filter(FAOBJ == "METASTATIC SITE", FASTRESC == "Y", FALOC != "OTHER") %>%
        dplyr::summarise(Num_Metast = dplyr::n(), .groups = "drop"),
      by = "USUBJID"
    ) %>%
    dplyr::mutate(
      Prim_Tum_Loc = ifelse(is.na(Prim_Tum_Loc), "NA", Prim_Tum_Loc),
      HPV = ifelse(is.na(HPV), "NA", HPV)
    )
  
  rt_set %>%
    dplyr::select(
      USUBJID, RFXSTDTC, Last_Recist, ACTARMCD, BOR, Sum_Diam_BL,
      PDL1DCPS, Prim_Tum_Loc, HPV, Num_Metast,
      CETUXIMAB_PRETRT, ONGOING
    )
}

getOH <- function(sdtm_set) {
  mh_tum <- sdtm_set$mh
  
  if ("MHSCAT" %in% colnames(mh_tum)) {
    mh_tum <- mh_tum %>% dplyr::filter(MHSCAT == "PRIMARY DIAGNOSIS")
  } else if ("MHCAT" %in% colnames(mh_tum)) {
    mh_tum <- mh_tum %>% dplyr::filter(MHCAT == "DISEASE HISTORY")
  }
  
  if (nrow(mh_tum) == 0) {
    return(tibble::tibble(
      USUBJID = character(), MHTERM = character(), D_TNM = character(),
      SP_D_TSTAGE = character(), SP_S_TSTAGE = character()
    ))
  }
  
  mh_tum %>%
    dplyr::mutate(
      MHTERM = stringr::str_trim(gsub("CLASSIFIC. OF TUMOR-", "", as.character(MHTERM))),
      D_TNM = sample(c("T1N1M0", "T2N1M0", "T3N2M0", "T4N2M1"), dplyr::n(), replace = TRUE),
      SP_D_TSTAGE = sample(c("T1", "T2", "T3", "T4"), dplyr::n(), replace = TRUE),
      SP_S_TSTAGE = sample(c("Stage II", "Stage III", "Stage IV"), dplyr::n(), replace = TRUE)
    ) %>%
    dplyr::select(USUBJID, MHTERM, D_TNM, SP_D_TSTAGE, SP_S_TSTAGE)
}

getBiomarkerSet <- function(sdtm_set) {
  if (!"yc" %in% names(sdtm_set)) {
    return(tibble::tibble(
      USUBJID = character(), PARACAT1 = character(), AVAL = numeric(),
      PARAM = character(), VISIT = character(), VISITDY = numeric(),
      LLOQ = numeric(), ULOQ = numeric()
    ))
  }
  
  sdtm_set$yc %>%
    dplyr::mutate(
      PARACAT1 = dplyr::coalesce(YCMETHOD, "UNKNOWN METHOD"),
      AVAL = YCSTRESN,
      PARAM = YCTEST,
      PARACAT1 = dplyr::if_else(PARACAT1 == "FLOW CYTOMETRY", "FACS Cellular Biomarker", PARACAT1),
      PARACAT1 = dplyr::if_else(PARACAT1 == "FLUORESCENT IMMUNOASSAY", "Soluble Biomarker", PARACAT1),
      PARACAT1 = dplyr::if_else(PARACAT1 == "ELECTROCHEMILUMINESCENCE IMMUNOASSAY", "Soluble Biomarker", PARACAT1)
    ) %>%
    dplyr::filter(!is.na(AVAL)) %>%
    dplyr::transmute(
      USUBJID = as.factor(USUBJID),
      PARACAT1 = as.factor(PARACAT1),
      AVAL = as.numeric(AVAL),
      PARAM = as.factor(PARAM),
      VISIT = as.factor(VISIT),
      VISITDY = as.numeric(VISITDY),
      LLOQ = as.numeric(YCLLOQ),
      ULOQ = as.numeric(YCULOQ)
    ) %>%
    dplyr::group_by(USUBJID, PARACAT1, PARAM, VISIT, VISITDY) %>%
    dplyr::summarise(
      AVAL = mean(AVAL, na.rm = TRUE),
      LLOQ = mean(LLOQ, na.rm = TRUE),
      ULOQ = mean(ULOQ, na.rm = TRUE),
      .groups = "drop"
    )
}

#READ AND PREPROCESS ONE STUDY
read_local_study <- function(study_name, data_root = file.path("resource", "data")) {
  study_path <- file.path(data_root, study_name)
  
  data_list <- list(
    ae = read_csv_safe(file.path(study_path, "ae.csv")),
    cm = read_csv_safe(file.path(study_path, "cm.csv")),
    dm = read_csv_safe(file.path(study_path, "dm.csv")),
    ds = read_csv_safe(file.path(study_path, "ds.csv")),
    ex = read_csv_safe(file.path(study_path, "ex.csv")),
    lb = read_csv_safe(file.path(study_path, "lb.csv")),
    mh = read_csv_safe(file.path(study_path, "mh.csv")),
    fa = read_csv_safe(file.path(study_path, "fa.csv")),
    sv = read_csv_safe(file.path(study_path, "sv.csv")),
    vs = read_csv_safe(file.path(study_path, "vs.csv")),
    rs = read_csv_safe(file.path(study_path, "rs.csv")),
    tr = read_csv_safe(file.path(study_path, "tr.csv")),
    yc = read_csv_safe(file.path(study_path, "yc.csv"))
  )
  
  data_list <- purrr::map(data_list, ~ {
    if ("STUDYID" %in% names(.x)) .x$STUDYID <- study_name
    .x
  })
  
  # Fill expected columns
  if (!"AETERM" %in% names(data_list$ae)) data_list$ae$AETERM <- data_list$ae$AEDECOD
  if (!"AEBODSYS" %in% names(data_list$ae)) data_list$ae$AEBODSYS <- data_list$ae$AESOC
  if (!"AETOXGR" %in% names(data_list$ae)) data_list$ae$AETOXGR <- 1L
  if (!"AECONTRT" %in% names(data_list$ae)) data_list$ae$AECONTRT <- "NONE"
  if (!"AEACN" %in% names(data_list$ae)) data_list$ae$AEACN <- "NONE"
  
  if (!"CMDECOD" %in% names(data_list$cm)) data_list$cm$CMDECOD <- data_list$cm$CMTRT
  if (!"CMCLAS" %in% names(data_list$cm)) data_list$cm$CMCLAS <- data_list$cm$CMCAT
  
  if (!"FAORRES" %in% names(data_list$fa)) data_list$fa$FAORRES <- data_list$fa$FASTRESC
  if (!"FATSTDTL" %in% names(data_list$fa)) data_list$fa$FATSTDTL <- ""
  if (!"FAOBJ" %in% names(data_list$fa)) data_list$fa$FAOBJ <- ""
  if (!"FALOC" %in% names(data_list$fa)) data_list$fa$FALOC <- ""
  
  if (!"ARM" %in% names(data_list$dm)) data_list$dm$ARM <- data_list$dm$ACTARM
  if (!"ARMCD" %in% names(data_list$dm)) data_list$dm$ARMCD <- data_list$dm$ACTARMCD
  if (!"SITEID" %in% names(data_list$dm)) data_list$dm$SITEID <- "001"
  if (!"ETHNIC" %in% names(data_list$dm)) data_list$dm$ETHNIC <- "NOT HISPANIC OR LATINO"
  if (!"COUNTRY" %in% names(data_list$dm)) data_list$dm$COUNTRY <- "USA"
  if (!"RFPENDTC" %in% names(data_list$dm)) data_list$dm$RFPENDTC <- data_list$dm$RFXENDTC
  
  if (!"DSSTDY" %in% names(data_list$ds)) data_list$ds$DSSTDY <- NA_real_
  
  if (!"LBDY" %in% names(data_list$lb)) data_list$lb$LBDY <- data_list$lb$VISITDY
  if (!"LBORRESU" %in% names(data_list$lb)) data_list$lb$LBORRESU <- ""
  if (!"LBSTRESU" %in% names(data_list$lb)) data_list$lb$LBSTRESU <- data_list$lb$LBORRESU
  if (!"LBCAT" %in% names(data_list$lb)) data_list$lb$LBCAT <- "CHEMISTRY"
  if (!"LBSCAT" %in% names(data_list$lb)) data_list$lb$LBSCAT <- "LAB"
  if (!"LBSPEC" %in% names(data_list$lb)) data_list$lb$LBSPEC <- "SERUM"
  if (!"LBDTC" %in% names(data_list$lb)) data_list$lb$LBDTC <- as.character(as.Date("2025-01-01") + data_list$lb$VISITDY)
  if (!"LBNRIND" %in% names(data_list$lb)) {
    data_list$lb$LBNRIND <- ifelse(data_list$lb$LBSTRESN < data_list$lb$LBSTNRLO, "LOW",
                                   ifelse(data_list$lb$LBSTRESN > data_list$lb$LBSTNRHI, "HIGH", "NORMAL"))
  }
  
  if (!"VSORRES" %in% names(data_list$vs) && "VSSTRESN" %in% names(data_list$vs)) {
    data_list$vs$VSORRES <- as.character(data_list$vs$VSSTRESN)
  }
  if (!"VSDY" %in% names(data_list$vs)) data_list$vs$VSDY <- data_list$vs$VISITDY
  if (!"VSSTRESU" %in% names(data_list$vs) && "VSORRESU" %in% names(data_list$vs)) {
    data_list$vs$VSSTRESU <- data_list$vs$VSORRESU
  }
  
  if (!"MHSTDY" %in% names(data_list$mh)) data_list$mh$MHSTDY <- NA_integer_
  if (!"MHENRTPT" %in% names(data_list$mh)) data_list$mh$MHENRTPT <- ""
  
  data_list_labels <- purrr::map(data_list, get_labels)
  
  # AE
  data_list$ae <- data_list$ae %>%
    dplyr::mutate(
      AETERM = dplyr::coalesce(AETERM, AEDECOD),
      AEBODSYS = dplyr::coalesce(AEBODSYS, AESOC),
      aestdtc_dt = parse_full_date_time(AESTDTC),
      aeendtc_dt = parse_full_date_time(AEENDTC),
      aeendtc_dt = dplyr::if_else(aestdtc_dt > aeendtc_dt, aestdtc_dt, aeendtc_dt),
      AETOXGRC = dplyr::coalesce(as.character(AETOXGR), as.character(1L)),
      aeser = (AESER == "Y")
    )
  
  # CM
  data_list$cm <- data_list$cm %>%
    dplyr::mutate(
      CMDECOD = dplyr::coalesce(CMDECOD, CMTRT),
      CMCAT = dplyr::coalesce(CMCAT, "CONCOMITANT"),
      CMCLAS = dplyr::coalesce(CMCLAS, CMCAT),
      cmstdtc_dt = parse_full_date_time(CMSTDTC),
      cmendtc_dt = parse_full_date_time(CMENDTC),
      cmendtc_dt = dplyr::if_else(cmstdtc_dt > cmendtc_dt, cmstdtc_dt, cmendtc_dt)
    )
  
  # DM
  data_list$dm <- data_list$dm %>%
    dplyr::mutate(
      USUBJID = factor(USUBJID),
      rfxstdtc_dt = parse_full_date_time(RFXSTDTC),
      rfxendtc_dt = parse_full_date_time(RFXENDTC),
      rficdtc_dt = parse_full_date_time(RFICDTC),
      rfxstdtc_d = parse_full_date(RFXSTDTC),
      rfxendtc_d = parse_full_date(RFXENDTC),
      rficdtc_d = parse_full_date(RFICDTC),
      rfpendtc_d = parse_full_date(RFPENDTC),
      ARM_DEFAULT = dplyr::coalesce(ACTARM, ARM),
      ARM = dplyr::if_else(is.na(ARM) | ARM == "", "Screen Failure", ARM),
      ARMCD = dplyr::if_else(is.na(ARMCD) | ARMCD == "", "SCRNFAIL", ARMCD),
      ACTARMCD = dplyr::if_else(is.na(ACTARMCD) | ACTARMCD == "", "SCRNFAIL", ACTARMCD),
      ACTARM = dplyr::if_else(is.na(ACTARM) | ACTARM == "", "Screen Failure", ACTARM),
      SEX = dplyr::if_else(is.na(SEX) | SEX == "", "NA", SEX),
      RACE = dplyr::if_else(is.na(RACE) | RACE == "", "NA", RACE),
      ETHNIC = dplyr::if_else(is.na(ETHNIC) | ETHNIC == "", "NA", ETHNIC),
      COUNTRY = dplyr::if_else(is.na(COUNTRY) | COUNTRY == "", "NA", COUNTRY)
    )
  
  # DS
  data_list$ds <- data_list$ds %>%
    dplyr::mutate(
      random_dt = dplyr::if_else(
        DSDECOD == "RANDOMIZED",
        parse_full_date_time(DSSTDTC),
        as.POSIXct(NA, tz = "UTC")
      ),
      random_dy = dplyr::if_else(
        DSDECOD == "RANDOMIZED",
        as.numeric(DSSTDY),
        NA_real_
      )
    )
  attributes(data_list$ds$random_dt)$label <- "Randomization Date/Time"
  attributes(data_list$ds$random_dy)$label <- "Randomization Day"
  
  # EX
  data_list$ex <- data_list$ex %>%
    dplyr::mutate(
      exstdtc_dt = parse_full_date_time(EXSTDTC),
      exendtc_dt = parse_full_date_time(EXENDTC)
    )
  
  # LB
  sorted_visits <- sort_visits(data_list$sv)
  
  data_list$lb <- data_list$lb %>%
    dplyr::left_join(
      smallest_ref_range(data_list$lb, "LBSTNRLO", "LBSTNRHI"),
      by = c("LBTESTCD" = "LBTESTCD")
    ) %>%
    dplyr::mutate(
      USUBJID = factor(USUBJID),
      LBCAT = factor(LBCAT),
      LBSCAT = factor(LBSCAT),
      LBTESTCD = factor(LBTESTCD),
      VISIT = factor(VISIT, levels = sorted_visits),
      lbdtc_d = as.Date(LBDTC),
      lborresnum = as.numeric(LBSTRESN),
      lbornrlonum = as.numeric(LBSTNRLO),
      lbornrhinum = as.numeric(LBSTNRHI),
      param = dplyr::if_else(LBSTRESU == "", LBTEST, paste0(LBTEST, " (", LBORRESU, ")")),
      VISIT = factor(as.character(VISIT), levels = c("SCREENING", setdiff(unique(as.character(VISIT)), "SCREENING")))
    )
  
  data_list$lb_lineplot <- data_list$lb %>%
    dplyr::filter(!is.na(VISITDY)) %>%
    dplyr::mutate(
      sort_cat = paste(LBTESTCD, LBSPEC, sep = "_"),
      param_lp = paste(LBTESTCD, LBSPEC, sep = "_")
    ) %>%
    dplyr::group_by(USUBJID, sort_cat, param_lp, VISIT) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  data_list$lb_lineplot <- dv.explorer.parameter::prefix_repeat_parameters(
    data_list$lb_lineplot,
    cat_var = "LBSCAT",
    par_var = "param_lp"
  )
  data_list_labels[["lb_lineplot"]] <- data_list_labels[["lb"]]
  
  # VS
  data_list$vs <- data_list$vs %>%
    dplyr::mutate(
      USUBJID = factor(USUBJID),
      VSTESTCD = factor(VSTESTCD),
      VISIT = factor(VISIT, levels = sorted_visits),
      vsorresnum = as.numeric(VSORRES),
      vsdtc_d = as.Date(VSDTC)
    )
  
  # restore labels before derived data
  data_list <- purrr::map2(data_list, data_list_labels, set_labels)
  
  # Derived data required by app
  data_list$resp <- getRespList(sdtm_set = data_list)
  data_list$oh <- getOH(sdtm_set = data_list)
  data_list$bm <- getBiomarkerSet(sdtm_set = data_list)
  
  data_list$bm <- dv.explorer.parameter::prefix_repeat_parameters(
    data_list$bm,
    cat_var = "PARACAT1",
    par_var = "PARAM"
  )
  
  attr(data_list, "arm_defaults") <- get_arm_defaults(data_list$dm)
  
  data_list <- purrr::map(data_list, function(x) {
    if (is.data.frame(x)) attr(x, "meta") <- list(mtime = Sys.time())
    x
  })
  
  data_list
}

#BUILD MULTI-STUDY LIST
study_list <- purrr::map(study_ids, ~ read_local_study(.x, data_root = data_root))
names(study_list) <- study_ids

#APP CONFIG
def_vars <- list(
  resp = c("USUBJID", "RFXSTDTC", "Last_Recist", "ACTARMCD", "BOR", "Sum_Diam_BL",
           "PDL1DCPS", "Prim_Tum_Loc", "HPV", "Num_Metast",
           "CETUXIMAB_PRETRT", "ONGOING"),
  ae = c("USUBJID", "AETERM", "AEDECOD", "AESER", "AEREL", "AESTDY", "AEENDY"),
  cm = c("USUBJID", "CMTRT", "CMINDC", "CMDOSE", "CMDOSU", "CMDOSFRQ", "CMROUTE"),
  dm = c("USUBJID", "RFXSTDTC", "AGE", "SEX", "RACE", "ACTARM", "COUNTRY"),
  ds = c("USUBJID", "DSDECOD", "DSSTDTC"),
  ex = c("USUBJID", "EXTRT", "EXDOSE", "EXSTDY", "EXENDY"),
  lb = c("USUBJID", "LBTESTCD", "LBSTRESN", "LBSTNRLO", "LBSTNRHI", "LBNRIND", "VISITNUM"),
  mh = c("USUBJID", "MHDECOD", "MHSTDY", "MHENRTPT"),
  oh = c("USUBJID", "MHTERM", "D_TNM", "SP_D_TSTAGE", "SP_S_TSTAGE"),
  vs = c("USUBJID", "VSTEST", "VSSTRESN", "VSORRESU", "VISITNUM")
)

all_visits <- unique(unlist(lapply(study_list, function(x) x$sv$VISIT)))
default_visits <- all_visits[!grepl("UNSCHED", all_visits)]
ref_study <- names(study_list)[1]

#MODULES
j_keys <- join_keys(
  join_key("dm", keys = c("STUDYID", "USUBJID")),
  join_key("ae", keys = c("STUDYID", "USUBJID")),
  join_key("dm", "ae", keys = c("STUDYID", "USUBJID"))
)

tm_t_events_mod <- tm_t_events(
  label = "AE Hierarchy Table",
  dataname = "ae",
  parentname = "dm",
  arm_var = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$dm, c("ARM", "ACTARM")),
    selected = "ACTARM"
  ),
  llt = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$ae, c("AETERM", "AEDECOD")),
    selected = c("AEDECOD")
  ),
  hlt = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$ae, c("AEBODSYS", "AESOC")),
    selected = "AEBODSYS"
  )
)

ae_hier <- dv.teal::mod_teal(
  module_id = "ae_hier",
  teal_module = tm_t_events_mod,
  j_keys = j_keys
)

j_keys <- join_keys(
  j_keys,
  join_key("cm", keys = c("STUDYID", "USUBJID")),
  join_key("dm", "cm", keys = c("STUDYID", "USUBJID"))
)

tm_t_events_cm <- tm_t_events(
  label = "CM Hierarchy Table",
  dataname = "cm",
  parentname = "dm",
  arm_var = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$dm, c("ARM", "ACTARM")),
    selected = "ACTARM"
  ),
  llt = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$cm, c("CMDECOD", "CMTRT")),
    selected = c("CMDECOD")
  ),
  hlt = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$cm, c("CMCAT", "CMCLAS")),
    selected = "CMCAT"
  )
)

cm_hier <- dv.teal::mod_teal(
  module_id = "cm_hier",
  teal_module = tm_t_events_cm,
  j_keys = j_keys
)

clinlines <- dv.clinlines::mod_clinical_timelines(
  module_id = "clinlines",
  basic_info = list(
    subject_level_dataset_name = "dm",
    trt_start_var = "rfxstdtc_dt",
    trt_end_var = "rfxendtc_dt",
    icf_date_var = "rficdtc_dt"
  ),
  mapping = list(
    ae = list(
      "Adverse Events" = list(
        start_dt_var = "aestdtc_dt",
        end_dt_var = "aeendtc_dt",
        detail_var = "AEDECOD"
      )
    )
  ),
  drug_admin = list(
    dataset_name = "ex",
    trt_var = "EXTRT",
    start_var = "exstdtc_dt",
    end_var = "exendtc_dt",
    detail_var = "EXTRT",
    label = "Drug Exposure",
    dose_var = "EXDOSE",
    dose_unit_var = "EXDOSU"
  ),
  default_plot_settings = list(
    x_param = "day", start_day = -5, boxheight_val = 60
  ),
  receiver_id = "papo"
)

papo <- dv.papo::mod_patient_profile(
  module_id = "papo",
  subject_level_dataset_name = "dm",
  subjid_var = "USUBJID",
  sender_ids = c("clinlines", "lb_lineplot", "bm_lineplot"),
  summary = list(
    vars = c("ARM", "SITEID", "RACE", "SEX", "RFICDTC", "RFXSTDTC", "RFXENDTC", "RFPENDTC"),
    column_count = 4
  ),
  listings = list(
    "Adverse Events" = list(dataset = "ae", default_vars = def_vars[["ae"]]),
    "Concomitant Medication" = list(dataset = "cm", default_vars = def_vars[["cm"]]),
    "Demographics" = list(dataset = "dm", default_vars = def_vars[["dm"]]),
    "Disposition" = list(dataset = "ds", default_vars = def_vars[["ds"]]),
    "Drug Exposure" = list(dataset = "ex", default_vars = def_vars[["ex"]]),
    "Laboratory" = list(dataset = "lb", default_vars = def_vars[["lb"]]),
    "Medical History" = list(dataset = "mh", default_vars = def_vars[["mh"]]),
    "Vital Signs" = list(dataset = "vs", default_vars = def_vars[["vs"]])
  ),
  plots = list(
    timeline_info = c(
      icf_date = "rficdtc_d",
      trt_start_date = "rfxstdtc_d",
      trt_end_date = "rfxendtc_d",
      part_end_date = "rfpendtc_d"
    ),
    range_plots = list(
      "Adverse Events Plot" = list(
        dataset = "ae",
        vars = c(
          start_date = "aestdtc_dt",
          end_date = "aeendtc_dt",
          decode = "AEDECOD",
          grading = "AETOXGRC",
          serious_ae = "aeser"
        ),
        tooltip = c(
          "Preferred Term: " = "AEDECOD",
          "Primary SOC: " = "AESOC",
          "Serious Event: " = "AESER",
          "<br><br>AE Start Date: " = "AESTDTC",
          "AE Stop Date: " = "AEENDTC",
          "AE Start Day: " = "AESTDY",
          "AE Stop Day: " = "AEENDY",
          "<br>CM or add. Treatment: " = "AECONTRT",
          "Action Taken with Study Treatment: " = "AEACN"
        )
      ),
      "Concomitant Medication" = list(
        dataset = "cm",
        vars = c(
          start_date = "cmstdtc_dt",
          end_date = "cmendtc_dt",
          decode = "CMDECOD",
          grading = "CMINDC"
        ),
        tooltip = c(
          "Standardized Medication Name: " = "CMDECOD",
          "Indication: " = "CMINDC",
          "<br>CM Dose: " = "CMDOSE",
          "CM Dose Unit: " = "CMDOSU",
          "Route of Administration: " = "CMROUTE",
          "<br>CM Start Date: " = "CMSTDTC",
          "CM End Date: " = "CMENDTC",
          "CM Start Day: " = "CMSTDY",
          "CM End Day: " = "CMENDY"
        )
      )
    ),
    value_plots = list(
      "Lab Plot" = list(
        dataset = "lb",
        vars = c(
          analysis_param = "param",
          analysis_val = "lborresnum",
          analysis_date = "lbdtc_d",
          analysis_indicator = "LBNRIND",
          range_low_limit = "lbornrlonum",
          range_high_limit = "lbornrhinum"
        ),
        tooltip = c(
          "Lab Parameter:" = "param",
          "Lab Test Day: " = "LBDY",
          "Lab Test Visit: " = "VISIT",
          "<br>High Limit: " = "lbornrhinum",
          "Lab Original Value: " = "lborresnum",
          "Lower Limit: " = "lbornrlonum",
          "Analysis Indicator: " = "LBNRIND"
        )
      )
    ),
    vline_vars = c(),
    vline_day_numbers = c("Study Treatment Start Day: Day 1" = 1)
  )
)

j_keys <- join_keys(
  j_keys,
  join_key("dm", keys = c("STUDYID", "USUBJID"))
)

tm_t_sum <- tm_t_summary(
  label = "Demographics Explorer",
  dataname = "dm",
  parentname = "dm",
  arm_var = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$dm, c("ARM", "ACTARM")),
    selected = "ACTARM"
  ),
  summarize_vars = choices_selected(
    choices = variable_choices(study_list[[ref_study]]$dm, c("AGE", "SEX", "RACE", "ETHNIC", "COUNTRY")),
    selected = "AGE"
  ),
  add_total = TRUE,
  useNA = "ifany"
)

dm_expl <- dv.teal::mod_teal(
  module_id = "dm_expl",
  teal_module = tm_t_sum,
  j_keys = j_keys
)

lb_lineplot <- dv.explorer.parameter::mod_lineplot(
  module_id = "lb_lineplot",
  bm_dataset_name = "lb_lineplot",
  group_dataset_name = "dm",
  cat_var = "LBSCAT",
  par_var = "param_lp",
  value_vars = "LBSTRESN",
  visit_vars = c("VISIT", "VISITDY"),
  default_visit_val = list("VISIT" = default_visits),
  subjid_var = "USUBJID",
  ref_line_vars = c("lbstnrlo_max", "lbstnrhi_min"),
  receiver_id = "papo"
)

bm_lineplot <- dv.explorer.parameter::mod_lineplot(
  module_id = "bm_lineplot",
  bm_dataset_name = "bm",
  group_dataset_name = "dm",
  cat_var = "PARACAT1",
  par_var = "PARAM",
  value_vars = "AVAL",
  visit_vars = c("VISIT", "VISITDY"),
  default_visit_val = list("VISIT" = default_visits),
  subjid_var = "USUBJID",
  receiver_id = "papo"
)

listing <- dv.listings::mod_listings(
  module_id = "listing",
  dataset_names = c("resp", "ae", "cm", "dm", "ds", "ex", "lb", "mh", "oh", "vs"),
  default_vars = def_vars
)

mod_list <- list(
  "AE Hierarchy Table" = ae_hier,
  "CM Hierarchy Table" = cm_hier,
  "Clinical Timelines" = clinlines,
  "Patient Profiles" = papo,
  "Demographic Explorer" = dm_expl,
  "Lab Line Plot" = lb_lineplot,
  "BM Line Plot" = bm_lineplot,
  "Listings" = listing
)

#LAUNCH APP
dv.manager::run_app(
  data = study_list,
  module_list = mod_list,
  filter_data = "dm",
  title = "Interactive Clinical Trial Data Platform"
)
