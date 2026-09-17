######################################################################################
## TITLE: Data preparation of CCU073 datasets
##
## Author: Genevieve Cezard, Wen Shi
######################################################################################
install.packages("R.utils", repos = "https://packages.sde.digital.nhs.uk/repository/cran-mirror/", dependencies = TRUE)
library(R.utils)
library(data.table)
library(tidyverse)

setwd("~/collab/CCU073_01")

# Read in data
df = fread('~data/ccu073_01_year2022_20241126.csv.gz')
#df = fread('~data/ccu073_01_year2022_20241213.csv.gz') # For HES APC pos1


##---------------------------------------------------------##
## 1. Create dataset with all required formatted variables ##
##---------------------------------------------------------##

#############################################
# 1.a. Reformat/Refactor existing variables #
#############################################

df[,sex:=factor(SEX,levels=c(1,2),labels=c('Male','Female'))]
df[,ethnicity:=factor(ETHNIC_CAT,levels=c("White","Asian or Asian British","Black, Black British, Caribbean or African", "Mixed or multiple ethnic groups","Other ethnic group",""),
                              labels=c('White','Asian','Black','Mixed','Other','Missing'))]
df[,deprivation:=factor(IMD_2019_DECILES,levels=10:1,
                          labels=rep(c("9_10_least_deprived","7_8","5_6","3_4","1_2_most_deprived"),each=2))]
df[,region:=factor(region,levels=c('London','South East','South West','East Midlands',
                                   'West Midlands',
                                   'East of England','North East', 'North West',
                                   'Yorkshire and The Humber'))]
df[,Smoking_status:=factor(Smoking_status,levels=c('Never','Ex','Current','Missing'))]
df[,vaccination_num:=factor(vaccination_status,levels=c('0','1','2','3+'))]
df[,covid:=factor(covid_flag,levels=c(0,1),labels=c('No','Yes'))]

# Select variables
df = df[,!c('CHUNK','SEX','ETHNIC_CAT','IMD_2019_DECILES','vaccination_status','covid_flag')]

############################
# 1.b. Create new variable #
############################

# Add study end date
df[,data_end := as.IDate('2024-06-30')]

############################################################
# 1.c. Create/derive new variables from existing variables #
############################################################

# Add BMI categories - once data is transferred
# cat option 1: underweight BMI < 18.5, normal weight 18-24, overweight 25-29, obese 30-34, severe obesity 35-39, morbid obesity 40+
# cat option 2: underweight BMI < 18.5, normal weight 18-24, overweight 25-29, obese 30+

df[cov_bmi_value<18.5,bmi_cat:='Underweight']
df[cov_bmi_value>=18.5&cov_bmi_value<25,bmi_cat:='Healthy']
df[cov_bmi_value>=25&cov_bmi_value<30,bmi_cat:='Overweight']
df[cov_bmi_value>=30,bmi_cat:='Obese']
df[is.na(cov_bmi_value),bmi_cat:='Missing']
df[,bmi_cat:=factor(bmi_cat,levels=c('Healthy','Underweight','Overweight','Obese','Missing'))]

# Calculate baseline age and add age categories
df[,baseline_age:=floor((baseline_date-DOB)/365.25)]
df[baseline_age<50,base_age_cat:='40_49']
df[baseline_age>=50 & baseline_age<60,base_age_cat:='50_59']
df[baseline_age>=60 & baseline_age<70,base_age_cat:='60_69']
df[baseline_age>=70 & baseline_age<80,base_age_cat:='70_79']
df[baseline_age>=80 & baseline_age<90,base_age_cat:='80_89']
df[baseline_age>=90,base_age_cat:='90_over']
df[,base_age_cat:=factor(base_age_cat,levels=c('40_49','50_59','60_69','70_79','80_89','90_over'))]

# Add egfr categories - from egfr records
df[cov_egfr_value<60,egfr_cat:='<60']
df[cov_egfr_value>=60&cov_egfr_value<90,egfr_cat:='60_89']
df[cov_egfr_value>=90,egfr_cat:='>=90']
df[is.na(cov_egfr_value),egfr_cat:='Missing']
df[,egfr_cat:=factor(egfr_cat,levels=c('<60','60_89','>=90','Missing'))]

# Convert Creatinine - Serum creatinine should be in mg/dL
df$Creatinine_converted <- df$cov_Creatinine_value / 88.4 #converts umol/L to mg/dL

# Calculate eGFR from Creatinine (cov_Creatinine_value) 
# (using the 2021 CKD-Epi equation)
df$cov_egfr_calculated_value <- ifelse(df$sex == "Female",
                         #women
                         142 * (ifelse(df$Creatinine_converted/0.7 <1, df$Creatinine_converted/0.7, 1))^(-0.329) *
                           (ifelse(df$Creatinine_converted/0.7 >1, df$Creatinine_converted/0.7, 1))^(-1.200) * (0.9938^df$baseline_age) * 1.012,
                         #men
                         142 * (ifelse(df$Creatinine_converted/0.9 <1, df$Creatinine_converted/0.9, 1))^(-0.411) *
                           (ifelse(df$Creatinine_converted/0.9 >1, df$Creatinine_converted/0.9, 1))^(-1.200) * (0.9938^df$baseline_age))
df[which(df$cov_egfr_calculated_value > 175),"cov_egfr_calculated_value"] <- NA
df$cov_egfr_calculated_value <- round(df$cov_egfr_calculated_value,1)
df$cov_egfr_calculated_date <- df$cov_Creatinine_date

# Note: eGFR in range 90-120 mL/min/1.73M2 is considered normal (stage 1), 60-89 early stage CKD (stage 2), 15-59 CKD (stage 3-4), <15 kidney failure (stage 5) 
df[cov_egfr_calculated_value<60,egfr_calc_cat:='<60']
df[cov_egfr_calculated_value>=60 & cov_egfr_calculated_value<90,egfr_calc_cat:='60_89']
df[cov_egfr_calculated_value>=90,egfr_calc_cat:='>=90']
df[is.na(cov_egfr_calculated_value),egfr_calc_cat:='Missing']
df[,egfr_calc_cat:=factor(egfr_calc_cat, levels=c('<60','60_89','>=90','Missing'))]

# Add glucose categories
df[cov_Glucose_level_value<7,glucose_cat:='<7']
df[cov_Glucose_level_value>=7&cov_Glucose_level_value<11,glucose_cat:='7_10']
df[cov_Glucose_level_value>=11,glucose_cat:='>=11']
df[is.na(cov_Glucose_level_value),glucose_cat:='Missing']
df[,glucose_cat:=factor(glucose_cat,levels=c('<7','7_10','>=11','Missing'))]


######################################################################################
# 1.d. Create outcome variables combining fatal and non-fatal for the 7 CVD outcomes #
######################################################################################

# create out_._date for 7 major cvd outcomes including nonfatal_._date and fatal_._date
for (i in c('Angina','CHD','ARR','HF','stroke','PVD','DVT')){
  df[,c(paste0('out_',i,'_date')):=pmin(get(paste0('nonfatal_',i,'_date')), get(paste0('fatal_',i,'_date')),na.rm=TRUE)]
  }

# create composite cvd
s1_prev = grep('^hist_.+_date$',colnames(df),value=TRUE)[!grep('^hist_.*_date$',colnames(df),value=TRUE)=='hist_Sudden_death_date']
s2_prev = s1_prev[!s1_prev %in% c("hist_Hypertensive_disease_date","hist_OtherCVD_date")]
df[,c('hist_s1_date','hist_s2_date','hist_s3_date'):= .(as.IDate(do.call(pmin,c(mget(s1_prev),na.rm=TRUE))),
                                                        as.IDate(do.call(pmin,c(mget(s2_prev),na.rm=TRUE))),
                                                        as.IDate(do.call(pmin,c(mget(s2_prev),na.rm=TRUE))))]

s1_inci = c(grep('^nonfatal_.+_date$',colnames(df),value=TRUE),grep('^fatal_.+_date$',colnames(df),value=TRUE))
s1_inci = s1_inci[!s1_inci %in% c("nonfatal_Hypertensive_disease_date",
                                  "nonfatal_OtherHD_date",
                                  "nonfatal_Arter_date",
                                  "nonfatal_OtherCVD_date",
                                  "nonfatal_Sudden_death_date")]
s2_inci = s1_inci[!s1_inci %in% c("nonfatal_Angina_date","fatal_Angina_date")]
s3_inci = s2_inci[!s2_inci %in% c("nonfatal_ARR_date","fatal_ARR_date","nonfatal_PVD_date","fatal_PVD_date","nonfatal_DVT_date","fatal_DVT_date")]

for (i in c('s1','s2','s3')){
  df[, c(paste0('out_',i,'_date')):= .(as.IDate(do.call(pmin,c(mget(get(paste0(i,'_inci'))),na.rm=TRUE))))]
}

# create history flags for 7 major cvd outcomes and 3 composite outcomes
for (i in c('Angina','CHD','ARR','HF','stroke','PVD','DVT','s1','s2','s3')){
  df[,paste0('hist_',i,'_flag'):=0]
  df[!is.na(get(paste0('hist_',i,'_date'))),paste0('hist_',i,'_flag'):=1]
}

##------------------------------------------------------------------##
## 2. Link to diabetes data and create exposure formatted variables ##
##------------------------------------------------------------------##

#####################################################
# 2.1. Link formatted dataset to DDSC diabetes data #
#####################################################

df_DDSC_DM_type <- data.table::fread("data/ccu073_01_ddsc_diabetes_type_2022_01_20241202.csv.gz", data.table = FALSE)

df_DDSC_DM_type_select <- df_DDSC_DM_type %>% select(PERSON_ID,date_of_diagnosis,age_at_diagnosis,out_diabetes)
df <- left_join(df, df_DDSC_DM_type_select, by = "PERSON_ID")

df$diabetes_status_DDSC <- ifelse(!is.na(df$out_diabetes), 
                                  ifelse(df$out_diabetes == "Type 1","T1DM",
                                         ifelse(df$out_diabetes == "Type 2","T2DM",
                                                ifelse(df$out_diabetes == "Other","Other","No diabetes"))
                                  ),
                                  "No diabetes")
df[,diabetes_status_DDSC:=factor(diabetes_status_DDSC,levels=c('T1DM','T2DM','Other','No diabetes'))]


df$age_T1DM_DDSC <- ifelse(df$out_diabetes == "Type 1",df$age_at_diagnosis,NA)
df$age_T2DM_DDSC <- ifelse(df$out_diabetes == "Type 2",df$age_at_diagnosis,NA)

df <- df %>% 
  mutate(age_DM_DDSC_10y_cat = cut(age_at_diagnosis,
                                   breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, Inf),
                                   labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
                                   right = FALSE),
         age_T1DM_DDSC_10y_cat = cut(age_T1DM_DDSC,
                                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, Inf),
                                     labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
                                     right = FALSE),
         age_T2DM_DDSC_10y_cat = cut(age_T2DM_DDSC,
                                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, Inf),
                                     labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
                                     right = FALSE),
  )

##################################
# 2.2. Define exposure variables #
##################################

#  DDSC T2DM age at diagnosis categorical
df[diabetes_status_DDSC == "No diabetes",expo_T2DM_DDSC_10y_age_cat:='No diabetes']
df[diabetes_status_DDSC == "T1DM",expo_T2DM_DDSC_10y_age_cat:='Missing']
df[diabetes_status_DDSC == "Other",expo_T2DM_DDSC_10y_age_cat:='Missing'] 
df[age_T2DM_DDSC_10y_cat %in% c("0-9","10-19","20-29"),expo_T2DM_DDSC_10y_age_cat:='0-29']
df[age_T2DM_DDSC_10y_cat %in% c("30-39","40-49","50-59","60-69","70-79"),expo_T2DM_DDSC_10y_age_cat:=age_T2DM_DDSC_10y_cat]
df[age_T2DM_DDSC_10y_cat %in% c("80-89","90+"),expo_T2DM_DDSC_10y_age_cat:='80+']
df$expo_T2DM_DDSC_10y_age_cat = fct_relevel(df$expo_T2DM_DDSC_10y_age_cat, "No diabetes")

#  DDSC T2DM age at diagnosis categorical - by sex
df[sex == "Male" & diabetes_status_DDSC == "No diabetes",expo_T2DM_DDSC_10y_age_cat_sex:='Male - No diabetes']
df[sex == "Male" & diabetes_status_DDSC == "T1DM",expo_T2DM_DDSC_10y_age_cat_sex:='Male - Missing']
df[sex == "Male" & diabetes_status_DDSC == "Other",expo_T2DM_DDSC_10y_age_cat_sex:='Male - Missing'] 
df[sex == "Male" & age_T2DM_DDSC_10y_cat %in% c("0-9","10-19","20-29"),expo_T2DM_DDSC_10y_age_cat_sex:='Male - 0-29']
df[sex == "Male" & age_T2DM_DDSC_10y_cat %in% c("30-39","40-49","50-59","60-69","70-79"),expo_T2DM_DDSC_10y_age_cat_sex:=paste0("Male - ",age_T2DM_DDSC_10y_cat)]
df[sex == "Male" & age_T2DM_DDSC_10y_cat %in% c("80-89","90+"),expo_T2DM_DDSC_10y_age_cat_sex:='Male - 80+']

df[sex == "Female" & diabetes_status_DDSC == "No diabetes",expo_T2DM_DDSC_10y_age_cat_sex:='Female - No diabetes']
df[sex == "Female" & diabetes_status_DDSC == "T1DM",expo_T2DM_DDSC_10y_age_cat_sex:='Female - Missing']
df[sex == "Female" & diabetes_status_DDSC == "Other",expo_T2DM_DDSC_10y_age_cat_sex:='Female - Missing'] 
df[sex == "Female" & age_T2DM_DDSC_10y_cat %in% c("0-9","10-19","20-29"),expo_T2DM_DDSC_10y_age_cat_sex:='Female - 0-29']
df[sex == "Female" & age_T2DM_DDSC_10y_cat %in% c("30-39","40-49","50-59","60-69","70-79"),expo_T2DM_DDSC_10y_age_cat_sex:=paste0("Female - ",age_T2DM_DDSC_10y_cat)]
df[sex == "Female" & age_T2DM_DDSC_10y_cat %in% c("80-89","90+"),expo_T2DM_DDSC_10y_age_cat_sex:='Female - 80+']

df$expo_T2DM_DDSC_10y_age_cat_sex = fct_relevel(df$expo_T2DM_DDSC_10y_age_cat_sex, "Male - No diabetes")

#  DDSC T2DM age at diagnosis continuous
df[diabetes_status_DDSC %in% c("No diabetes","T1DM","Other"),expo_T2DM_DDSC_age10c50:=0]
df[diabetes_status_DDSC == "T2DM",expo_T2DM_DDSC_age10c50:=(age_T2DM_DDSC-50)/10]

df[diabetes_status_DDSC %in% c("No diabetes","T1DM","Other"),expo_T2DM_DDSC_age10c40:=0]
df[diabetes_status_DDSC == "T2DM",expo_T2DM_DDSC_age10c40:=(age_T2DM_DDSC-40)/10]

df[diabetes_status_DDSC %in% c("No diabetes","T1DM","Other"),expo_T2DM_DDSC_age10c60:=0]
df[diabetes_status_DDSC == "T2DM",expo_T2DM_DDSC_age10c60:=(age_T2DM_DDSC-60)/10]

#  DDSC T2DM history
df[diabetes_status_DDSC %in% c("No diabetes","T1DM","Other"),expo_T2DM_DDSC_hist:="No"]
df[diabetes_status_DDSC == "T2DM",expo_T2DM_DDSC_hist:="Yes"]

# For sensitivity analyses with T1DM
#  DDSC T1DM age at diagnosis categorical
df[diabetes_status_DDSC == "No diabetes",expo_T1DM_DDSC_10y_age_cat:='No diabetes']
df[diabetes_status_DDSC == "T2DM",expo_T1DM_DDSC_10y_age_cat:='Missing']
df[diabetes_status_DDSC == "Other",expo_T1DM_DDSC_10y_age_cat:='Missing'] 
df[age_T1DM_DDSC_10y_cat %in% c("0-9","10-19","20-29"),expo_T1DM_DDSC_10y_age_cat:='0-29']
df[age_T1DM_DDSC_10y_cat %in% c("30-39","40-49","50-59","60-69","70-79"),expo_T1DM_DDSC_10y_age_cat:=age_T1DM_DDSC_10y_cat]
df[age_T1DM_DDSC_10y_cat %in% c("80-89","90+"),expo_T1DM_DDSC_10y_age_cat:='80+']
df$expo_T1DM_DDSC_10y_age_cat = fct_relevel(df$expo_T1DM_DDSC_10y_age_cat, "No diabetes")

#  DDSC T1DM age at diagnosis categorical - version 2 after discussion on Dec 10 2024
df[diabetes_status_DDSC == "No diabetes",expo_T1DM_DDSC_10y_age_cat_50plus:='No diabetes']
df[diabetes_status_DDSC == "T2DM",expo_T1DM_DDSC_10y_age_cat_50plus:='Missing']
df[diabetes_status_DDSC == "Other",expo_T1DM_DDSC_10y_age_cat_50plus:='Missing'] 
df[age_T1DM_DDSC_10y_cat %in% c("0-9","10-19","20-29","30-39","40-49"),expo_T1DM_DDSC_10y_age_cat_50plus:=age_T1DM_DDSC_10y_cat]
df[age_T1DM_DDSC_10y_cat %in% c("50-59","60-69","70-79","80-89","90+"),expo_T1DM_DDSC_10y_age_cat_50plus:='50+']
df$expo_T1DM_DDSC_10y_age_cat_50plus = fct_relevel(df$expo_T1DM_DDSC_10y_age_cat_50plus, "No diabetes")

#  DDSC T1DM age at diagnosis continuous
df[diabetes_status_DDSC %in% c("No diabetes","T2DM","Other"),expo_T1DM_DDSC_age10c50:=0]
df[diabetes_status_DDSC == "T1DM",expo_T1DM_DDSC_age10c50:=(age_T1DM_DDSC-50)/10]

#  DDSC T1DM history
df[diabetes_status_DDSC %in% c("No diabetes","T2DM","Other"),expo_T1DM_DDSC_hist:="No"]
df[diabetes_status_DDSC == "T1DM",expo_T1DM_DDSC_hist:="Yes"]

######################################################
# 2.3. Correct HbA1c using age at diabetes diagnosis #
######################################################

# Convert HbA1C from % to mmol/mol
df <- df %>% mutate(
    hba1c_value = ifelse( cov_hba1c_value <=20, 
                          ifelse( cov_hba1c_value <=15, round(cov_hba1c_value*10.93 - 23.5, 1), # i.e. cov_hba1c_value in ]0,15]
                                  ifelse(!is.na(age_at_diagnosis),round(cov_hba1c_value*10.93 - 23.5, 1),round(cov_hba1c_value,1)) # i.e. cov_hba1c_value in ]15,20]
                          ),
                          round(cov_hba1c_value,1)), # i.e. cov_hba1c_value in ]20;195]
    # Correction applied after conversion based on Carmen's approach
    hba1c_value = ifelse( (hba1c_value>= 10 | hba1c_value <= 195),hba1c_value, NA)
  )

# Add HbA1c categories based on corrected HbA1c values and with more granular HbA1c categories and labels
df[hba1c_value<39,hba1c_cat:='<39']
df[hba1c_value>=39&hba1c_value<48,hba1c_cat:='39_47']
df[hba1c_value>=48,hba1c_cat:=ifelse(hba1c_value<53,'48_52',
                      ifelse(hba1c_value<70,'53_69','>=70'))]
df[is.na(hba1c_value),hba1c_cat:='Missing']
df[,hba1c_cat:=factor(hba1c_cat,levels=c('<39','39_47','48_52','53_69','>=70','Missing'))]

###################### Recalculate age at diagnosis of diabetes #######################################
df[, age_at_diagnosis_2:=floor((as.IDate(date_of_diagnosis)-as.IDate(DOB))/365.25)]
df[,age_T2DM_DDSC_2:=ifelse(out_diabetes == "Type 2",age_at_diagnosis_2,NA)]
df[, 'age_T2DM_DDSC_10y_cat_2' := cut(age_T2DM_DDSC_2,
                                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, Inf),
                                     labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
                                     right = FALSE)]
df[diabetes_status_DDSC == "No diabetes",expo_T2DM_DDSC_10y_age_cat_2:='No diabetes']
df[diabetes_status_DDSC %in% c("T1DM","Other"),expo_T2DM_DDSC_10y_age_cat_2:='Missing']
df[age_T2DM_DDSC_10y_cat_2 %in% c("0-9","10-19","20-29"),expo_T2DM_DDSC_10y_age_cat_2:='0-29']
df[age_T2DM_DDSC_10y_cat_2 %in% c("30-39","40-49","50-59","60-69","70-79"),expo_T2DM_DDSC_10y_age_cat_2:=age_T2DM_DDSC_10y_cat_2]
df[age_T2DM_DDSC_10y_cat_2 %in% c("80-89","90+"),expo_T2DM_DDSC_10y_age_cat_2:='80+']
df[,expo_T2DM_DDSC_10y_age_cat_2 := factor(df$expo_T2DM_DDSC_10y_age_cat_2, levels=c("No diabetes",'0-29','30-39','40-49','50-59','60-69','70-79','80+','Missing'))]

  
# Add duration of diabetes
df[,'duration_T2DM':=.(baseline_age-age_T2DM_DDSC_2)]
df[,'duration_T2DM_cat':=.(ifelse(duration_T2DM<5,'<5',
                                ifelse(duration_T2DM<10,'5_9',
                                      ifelse(duration_T2DM<15,'10_14',
                                             ifelse(duration_T2DM<20,'15_19',ifelse(duration_T2DM>=20,'>=20',NA))))))]
df[is.na(duration_T2DM_cat)&diabetes_status_DDSC == "No diabetes",duration_T2DM_cat:='No diabetes']
df[is.na(duration_T2DM_cat)&diabetes_status_DDSC %in% c("T1DM","Other"),duration_T2DM_cat:='Missing']
df[,duration_T2DM_cat:=factor(duration_T2DM_cat,levels=c('No diabetes','<5','5_9','10_14','15_19','>=20','Missing'))]

#  Duration T2DM continuous per decade
df[diabetes_status_DDSC %in% c("No diabetes","T1DM","Other"),duration_T2DM_10y:=0]
df[diabetes_status_DDSC == "T2DM",duration_T2DM_10y:=duration_T2DM/10]

##---------------------------------------##
## 3. Save final dataset(s) for analysis ##
##---------------------------------------##
# Save out
#saveRDS(df, '~data/cohort_2022_20250630_pos1_duration.rds')
saveRDS(df, '~data/cohort_2022_20251017.rds')
