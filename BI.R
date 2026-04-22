# 1) Install packages 
#if (!require("remotes")) install.packages("remotes")
d#avinci_packages <- c(
#  "dv.loader", 
#  "dv.manager",
#  "dv.bookman", 
#  "dv.clinlines", 
##  "dv.edish", 
##  "dv.explorer.parameter", 
#  "dv.listings",
#  "dv.papo",
#  "dv.spiderplot",
#  "dv.swimmerplot",
#  "dv.tables",
#  "dv.teal"
#)

#for (pkg in davinci_packages) {
#  remotes::install_github(paste0("Boehringer-Ingelheim/", pkg), upgrade = TRUE)
#}

# 2) Load libraries
library(dplyr)
library(tibble)
library(lubridate)
library(stringr)
library(tidyr)
library(dv.manager)
library(dv.listings)
library(dv.edish)
library(dv.teal)
library(dv.clinlines)

# 3) Create mock SDTM-like data
set.seed(123)
LBSTRESN = case_when(
  LBTEST == "ALT" ~ round(rlnorm(n(), log(30), 0.45), 1),
  LBTEST == "AST" ~ round(rlnorm(n(), log(28), 0.45), 1),
  LBTEST == "BILI" ~ round(rlnorm(n(), log(0.8), 0.35), 2),
  LBTEST == "ALP" ~ round(rlnorm(n(), log(85), 0.25), 1),
  TRUE ~ NA_real_
),
LBNRIND = case_when(
  LBSTRESN < LBSTNRLO ~ "LOW",
  LBSTRESN > LBSTNRHI ~ "HIGH",
  TRUE ~ "NORMAL"
)
) %>%
  mutate(
    LBORRES = as.character(LBSTRESN),
    LBSTRESU = case_when(
      LBTEST %in% c("ALT", "AST", "ALP") ~ "U/L",
      LBTEST == "BILI" ~ "mg/dL",
      TRUE ~ NA_character_
    )
  )

# simple date variables for timeline module
sv <- expand.grid(
  USUBJID = usubjid,
  VISIT = visits,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(dm %>% select(USUBJID, STUDYID, RFXSTDTC), by = "USUBJID") %>%
  left_join(visit_lookup, by = "VISIT")

# -----------------------------
# 4) Convert and label data
# -----------------------------
dm <- dv.listings::convert_data(dm)
ae <- dv.listings::convert_data(ae)
cm <- dv.listings::convert_data(cm)
lb <- dv.listings::convert_data(lb)
sv <- dv.listings::convert_data(sv)

attributes(dm)$label <- "Subject Level"
attributes(ae)$label <- "Adverse Events"
attributes(cm)$label <- "Concomitant Medications"
attributes(lb)$label <- "Laboratory Data"
attributes(sv)$label <- "Subject Visits"

# For clinical timelines, create date versions if needed
# Adjust to match your real module expectations
if (!"rfxstdtc_dt" %in% names(dm)) {
  dm <- dm %>% mutate(
    rfxstdtc_dt = as.Date(RFXSTDTC),
    rfxendtc_dt = as.Date(RFXENDTC),
    rficdtc_dt = as.Date(RFICDTC)
  )
}

selected = "ACTARM"
),
llt = teal.widgets::choices_selected(
  choices = c("AETERM", "AEDECOD"),
  selected = "AEDECOD"
),
hlt = teal.widgets::choices_selected(
  choices = c("AESOC", "AEBODSYS"),
  selected = "AESOC"
)
),
j_keys = j_keys
)

# CM hierarchy table
module_list[["CM Hierarchy"]] <- dv.teal::mod_teal(
  module_id = "cm_hier",
  teal_module = teal.modules.clinical::tm_t_events(
    label = "CM Hierarchy Table",
    dataname = "cm",
    parentname = "dm",
    arm_var = teal.widgets::choices_selected(
      choices = c("ARM", "ACTARM"),
      selected = "ACTARM"
    ),
    llt = teal.widgets::choices_selected(
      choices = c("CMDECOD", "CMTRT"),
      selected = "CMDECOD"
    ),
    hlt = teal.widgets::choices_selected(
      choices = c("CMCAT", "CMCLAS"),
      selected = "CMCAT"
    )
  ),
  j_keys = j_keys
)

# Optional: Clinical timelines
# This is the module most likely to need adjustment depending on installed package version.
# Uncomment and adapt if your environment supports it.
# module_list[["Clinical Timelines"]] <- dv.clinlines::mod_clinical_timelines(
#   module_id = "clinlines",
#   basic_info = list(
#     subject_level_dataset_name = "dm",
#     trt_start_var = "rfxstdtc_dt",
#     trt_end_var = "rfxendtc_dt",
#     icf_date_var = "rficdtc_dt"
#   )
# )

# -----------------------------
# 8) Run app
# -----------------------------
mock_data <- list(
  dm = dm,
  ae = ae,
  cm = cm,
  lb = lb,
  sv = sv
)

dv.manager::run_app(
  data = list("Mock_Study" = mock_data),
  module_list = module_list,
  filter_data = "dm"
)
