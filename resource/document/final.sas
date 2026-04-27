ods pdf
file = "/home/u64113323/sw3455_final.pdf"
startpage=yes;

/*Import data*/
proc import datafile="/home/u64113323/cognition.xlsx"
    out=cognition
    dbms=xlsx
    replace;
    sheet="Data";
    getnames=yes;
run;

data cognition_mets; 
    set cognition; 
    /* Elevated blood pressure */
    bp_high = (meansbp >= 130 or meandbp >= 85 or bp_med = 1); 
    /* Elevated fasting glucose */ 
    glucose_high = (glucose >= 100 or diabetes_med = 1); 
    /* High waist circumference */ 
    if sex = 2 then waist_high = (waist >= 89); 
    else if sex = 1 then waist_high = (waist >= 102); 
    /* Elevated triglycerides level */ 
    trig_high = (trig >= 150 or fibrate_med = 1); 
    /* Abnormal HDL cholesterol */ 
    if sex = 2 then hdl_low = (hdl < 50 or lipid_med = 1); 
    else if sex = 1 then hdl_low = (hdl < 40 or lipid_med = 1); 

run; 

/*Get risk score*/ 
data cognition_mets; 
	set cognition_mets; 
	risk_score = sum( bp_high, glucose_high, waist_high, trig_high, hdl_low); 
	run;

/*classified MetS*/ 
data cognition_mets;
    set cognition_mets;
    if risk_score >= 3 then MetS = 1;
    else MetS = 0;
run;

/*label variable*/ 
data cognition_mets;
    set cognition_mets;
    if IPAQtotal in (7777, 8888, 9999) then IPAQtotal = .;
    label
        srt_tot = "Selective Reminding Test"
        cfl = "C-F-L Letter Fluency Test"
        cat = "Category Fluency Test"
        ctt1 = "Color Trails Test 1"
        ctt2= "Color Trails Test 2"
        MetS = "Metabolic Syndrome Status"
        IPAQtotal = "Physical Activity"
        promis = "Depression";
run;

/*check missing value*/ 
proc means data=cognition_mets n nmiss; 
	var 
	age sex education IPAQtotal promis srt_tot cfl cat ctt1 ctt2 MetS bp_high waist_high trig_high hdl_low glucose_high; 
	run; 


/*Distribution outcome*/ 
proc univariate data=cognition_mets; 
	var srt_tot cfl cat ctt1 ctt2; 
	histogram / normal; 
	qqplot / normal(mu=est sigma=est);
	run;

data cognition_long;
    set cognition_mets;
    length outcome $30 value 8;

    outcome = "SRT Total Score";        value = srt_tot; output;
    outcome = "C-F-L Letter Fluency";   value = cfl;     output;
    outcome = "Category Fluency";       value = cat;     output;
    outcome = "Color Trails Test 1";    value = ctt1;    output;
    outcome = "Color Trails Test 2";    value = ctt2;    output;

    keep outcome value;
run;
ods graphics on;

title "Figure 1. Distribution of Neuropsychological Test Scores";

proc sgpanel data=cognition_long;
    panelby outcome / columns=3 rows=2 spacing=10;
    histogram value / transparency=0.2;
    density value / type=normal;
    colaxis label="Score";
    rowaxis label="Percent";
run;

title;
ods graphics off;

/*make format for later analyze*/ 
proc format;
    value sexfmt 1 = "Male" 2 = "Female";
    value yesnofmt 0 = "No" 1 = "Yes";
run;

data cognition_mets;
    set cognition_mets;
    format 
        sex sexfmt.
        IPAQtotal promis MetS bp_high waist_high trig_high hdl_low glucose_high yesnofmt.;
run;

%let cont_vars = age promis education IPAQtotal srt_tot cfl cat ctt1 ctt2;

%let cat_vars  = sex bp_high glucose_high waist_high trig_high hdl_low;

/*Continuous variables*/
proc tabulate data=cognition_mets missing;
    class MetS;
    var age promis education IPAQtotal srt_tot cfl cat ctt1 ctt2;
    format MetS yesnofmt.;
    table
        (age promis education IPAQtotal srt_tot cfl cat ctt1 ctt2),
        (all MetS) *
        (n mean std median q1 q3);
run;

/* P-values for continuous vars (Wilcoxon rank-sum) */
%macro wilcox_pvals(varlist);
  %local i v;
  %let i=1;

  %do %while(%scan(&varlist,&i) ne );
    %let v=%scan(&varlist,&i);

    ods output WilcoxonTest=wil_&v;
    proc npar1way data=cognition_mets wilcoxon;
      class MetS;
      var &v;
    run;
    ods output close;

    proc print data=wil_&v noobs;
    run;

    %let i=%eval(&i+1);
  %end;
%mend;

%wilcox_pvals(&cont_vars);



/*Categorical variables*/
proc tabulate data=cognition_mets missing;
    class 
        MetS
        sex
        bp_high
        glucose_high
        waist_high
        trig_high
        hdl_low;

    format 
        sex sexF.
        MetS bp_high glucose_high waist_high trig_high hdl_low yesnoF.;

    table
        (sex
         bp_high
         glucose_high
         waist_high
         trig_high
         hdl_low),
        (all MetS) *
        (n colpctn) * f=8.1;
run;

/*p value for categorical,chisquare */
%macro chisq_pvals(varlist);
  %local i v;
  %let i=1;

  %do %while(%scan(&varlist,&i) ne );
    %let v=%scan(&varlist,&i);

    ods select none;
    ods output ChiSq=chisq_&v;
    proc freq data=cognition_mets;
      tables &v * MetS / chisq;
    run;
    ods select all;

    proc print data=chisq_&v noobs;
    run;

    %let i=%eval(&i+1);
  %end;
%mend;

%chisq_pvals(&cat_vars);

/*1. Are there differences in SRT scores by MetS status?*/
/*Unadjusted Analysis*/
/* Because SRT is not normally distributed, so I choice to use Wilcoxon rank-sum test*/
proc npar1way data=cognition_mets wilcoxon;
    class MetS;
    var srt_tot;
run;
/*adjusted Analysis*/
/* I used linear regression to estimate adjusted associations between MetS and SRT performance, controlling for age, sex, education, depressive symptoms (PROMIS), and physical activity (IPAQ total score).*/
proc glm data=cognition_mets;
    class MetS;
    model srt_tot = MetS age sex education promis IPAQtotal;
run;
quit;

/*I tried stepwise*/
proc reg data=cognition_mets;
    model srt_tot = MetS age sex education promis IPAQtotal/ selection=STEPWISE;
run;
quit;


/*2.Are there differences in other neuropsychological exams by MetS status?*/

/*Define outcome list for secondary aim so macro would make easier*/
%let neuro_outcomes = cfl cat ctt1 ctt2;

/*UNADJUSTED*/
/*t-tests (means by group, assumes approx normal) */
proc ttest data=cognition_mets;
  class MetS;
  var &neuro_outcomes;
run;

/*Wilcoxon rank-sum tests (nonparametric) */
%macro wilcox_list(varlist);
  %local i v;
  %let i=1;
  %do %while(%scan(&varlist,&i) ne );
    %let v=%scan(&varlist,&i);
    proc npar1way data=cognition_mets wilcoxon;
      class MetS;
      var &v;
    run;

    %let i=%eval(&i+1);
  %end;
%mend;

%wilcox_list(&neuro_outcomes);


/*ADJUSTED LINEAR REGRESSION Adjust for: age sex education promis IPAQtotal*/

%macro adj_glm_list(varlist);
  %local i v;
  %let i=1;

  %do %while(%scan(&varlist,&i) ne );
    %let v=%scan(&varlist,&i);

    title "Adjusted linear regression (PROC GLM): &v ~ MetS + age + sex + education + promis + IPAQtotal";

    proc glm data=cognition_mets;
      class MetS sex;           /* sex and MetS are categorical */
      model &v = MetS age sex education promis IPAQtotal / solution;
      /* Type III p-values are standard for adjusted inference */
    run;
    quit;

    %let i=%eval(&i+1);
  %end;
%mend;

%adj_glm_list(&neuro_outcomes);

 /* 3.Are your results consistent when using MetS risk score as a continuous exposure
rather than categorical MetS status?*/
 /*unadjust*/
proc glm data=cognition_mets;
    model srt_tot = risk_score;
run;
quit;
 /*adjust*/
proc glm data=cognition_mets;
    class sex;
    model srt_tot = risk_score age sex education promis IPAQtotal;
run;
quit;


ods graphics on;

title "Figure 2. Selective Reminding Test Scores Across Metabolic Syndrome Risk Score Levels";

proc sgplot data=cognition_mets;
  scatter x=risk_score y=srt_tot / markerattrs=(symbol=circle) transparency=0.2;
  reg x=risk_score y=srt_tot / cli clm;  /* clm = mean CI band; cli = prediction limits */
  xaxis label="risk_score";
  yaxis label="Selective Reminding Test";
run;

title;
ods graphics off;


 /* 4.Are there different associations with your outcomes when testing your MetS
scores individually? Are some more strongly related to your primary outcome
than others?*/
 /*unadjust*/
proc glm data=cognition_mets;
  model srt_tot = bp_high glucose_high waist_high trig_high hdl_low;
run;
quit;


 /*adjust*/
proc glm data=cognition_mets;
  class sex;
  model srt_tot = bp_high glucose_high waist_high trig_high hdl_low age sex education promis IPAQtotal;
run;
quit;


ods pdf close;
