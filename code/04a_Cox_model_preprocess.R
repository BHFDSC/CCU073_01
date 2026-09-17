######################################################################################
## TITLE: Prepare datasets for Cox model
##
## Author: Genevieve Cezard
######################################################################################
library(data.table)
library(tidyverse)
library(broom)

output_dir = "~output/"

#-------------------------------------------------------------#
# Preprocess for Cox including duration of diabetes and HbA1C #
#-------------------------------------------------------------#
# Dataset with duration of diabetes and HbA1c new categorisation
df_duration = readRDS('~data/cohort_2022_20250509_duration.rds')
df_duration_pos1 = readRDS('~data/cohort_2022_20250630_pos1_duration.rds')

cox_preprocess_duration = function(df,df_name){
  
  # Reduce dataset from 178 variables to 49 variables used in Cox
  df <- df %>%
    dplyr::select('PERSON_ID','DOB','DOD','baseline_date','data_end',
                  'out_Angina_date','out_CHD_date','out_ARR_date','out_HF_date','out_stroke_date','out_PVD_date','out_DVT_date',
                  'out_s1_date','out_s2_date','out_s3_date',
                  'hist_Angina_date','hist_CHD_date','hist_ARR_date','hist_HF_date','hist_stroke_date','hist_PVD_date','hist_DVT_date',
                  'hist_s1_date','hist_s2_date','hist_s3_date',
                  'sex','baseline_age','base_age_cat','ethnicity','region','deprivation','Smoking_status','bmi_cat','vaccination_num',
                  'diabetes_status_DDSC','out_diabetes','age_at_diagnosis', 
                  'age_T1DM_DDSC', 'age_T1DM_DDSC_10y_cat','expo_T1DM_DDSC_10y_age_cat','expo_T1DM_DDSC_10y_age_cat_50plus','expo_T1DM_DDSC_hist','expo_T1DM_DDSC_age10c50',
                  'age_T2DM_DDSC', 'age_T2DM_DDSC_10y_cat','expo_T2DM_DDSC_10y_age_cat','expo_T2DM_DDSC_hist','expo_T2DM_DDSC_age10c50',
                  'covid',
                  'duration_T2DM_cat','duration_T2DM_10y','hba1c_cat')
  
  
  # Change format for all dates
  df$DOD <-as.Date(df$DOD, format = "%d/%m/%Y")
  df$data_end <-as.Date(df$data_end, format = "%d/%m/%Y")
  df$baseline_date <-as.Date(df$baseline_date, format = "%d/%m/%Y")
  
  df$out_Angina_date <-as.Date(df$out_Angina_date, format = "%d/%m/%Y")
  df$out_CHD_date <-as.Date(df$out_CHD_date, format = "%d/%m/%Y")
  df$out_ARR_date <-as.Date(df$out_ARR_date, format = "%d/%m/%Y")
  df$out_HF_date <-as.Date(df$out_HF_date, format = "%d/%m/%Y")
  df$out_stroke_date <-as.Date(df$out_stroke_date, format = "%d/%m/%Y")
  df$out_PVD_date <-as.Date(df$out_PVD_date, format = "%d/%m/%Y")
  df$out_DVT_date <-as.Date(df$out_DVT_date, format = "%d/%m/%Y")
  df$out_s1_date <-as.Date(df$out_s1_date, format = "%d/%m/%Y")
  df$out_s2_date <-as.Date(df$out_s2_date, format = "%d/%m/%Y")
  df$out_s3_date <-as.Date(df$out_s3_date, format = "%d/%m/%Y")
  
  df$hist_Angina_date <-as.Date(df$hist_Angina_date, format = "%d/%m/%Y")
  df$hist_CHD_date <-as.Date(df$hist_CHD_date, format = "%d/%m/%Y")
  df$hist_ARR_date <-as.Date(df$hist_ARR_date, format = "%d/%m/%Y")
  df$hist_HF_date <-as.Date(df$hist_HF_date, format = "%d/%m/%Y")
  df$hist_stroke_date <-as.Date(df$hist_stroke_date, format = "%d/%m/%Y")
  df$hist_PVD_date <-as.Date(df$hist_PVD_date, format = "%d/%m/%Y")
  df$hist_DVT_date <-as.Date(df$hist_DVT_date, format = "%d/%m/%Y")
  df$hist_s1_date <-as.Date(df$hist_s1_date, format = "%d/%m/%Y")
  df$hist_s2_date <-as.Date(df$hist_s2_date, format = "%d/%m/%Y")
  df$hist_s3_date <-as.Date(df$hist_s3_date, format = "%d/%m/%Y")
  
  # Save dataset ready for Cox analysis
  saveRDS(df, paste0("~data/",df_name,"_forCox.rds"))
}
cox_preprocess_duration(df_duration,"cohort_2022_20250509_duration")
cox_preprocess_duration(df_duration_pos1,"cohort_2022_20250630_pos1_duration")

#-----------------------------------------------------------------------------#
# Preprocess for Cox including extra age at T2DM diagnosis exposure variables #
#-----------------------------------------------------------------------------#
# Dataset with duration of diabetes and HbA1c new categorisation
df_v3 = readRDS('~data/cohort_2022_20251017.rds')

cox_preprocess_v3 = function(df,df_name){
  
  # Reduce dataset from 178 variables to 49 variables used in Cox
  df <- df %>%
    dplyr::select('PERSON_ID','DOB','DOD','baseline_date','data_end',
                  'out_Angina_date','out_CHD_date','out_ARR_date','out_HF_date','out_stroke_date','out_PVD_date','out_DVT_date',
                  'out_s1_date','out_s2_date','out_s3_date',
                  'hist_Angina_date','hist_CHD_date','hist_ARR_date','hist_HF_date','hist_stroke_date','hist_PVD_date','hist_DVT_date',
                  'hist_s1_date','hist_s2_date','hist_s3_date',
                  'sex','baseline_age','base_age_cat','ethnicity','region','deprivation','Smoking_status','bmi_cat','vaccination_num',
                  'diabetes_status_DDSC','out_diabetes','age_at_diagnosis', 
                  'age_T1DM_DDSC', 'age_T1DM_DDSC_10y_cat','expo_T1DM_DDSC_10y_age_cat','expo_T1DM_DDSC_10y_age_cat_50plus','expo_T1DM_DDSC_hist','expo_T1DM_DDSC_age10c50',
                  'age_T2DM_DDSC', 'age_T2DM_DDSC_10y_cat','expo_T2DM_DDSC_10y_age_cat','expo_T2DM_DDSC_hist','expo_T2DM_DDSC_age10c50',
                  'expo_T2DM_DDSC_age10c40', 'expo_T2DM_DDSC_age10c60', 'expo_T2DM_DDSC_10y_age_cat_sex',
                  'covid',
                  'duration_T2DM_cat','duration_T2DM_10y','hba1c_cat')
  
  # Add HbA1c flag
  df$hba1c_flag <- ifelse(df$hba1c_cat == "Missing", "No","Yes")
  
  # Change format for all dates
  df$DOD <-as.Date(df$DOD, format = "%d/%m/%Y")
  df$data_end <-as.Date(df$data_end, format = "%d/%m/%Y")
  df$baseline_date <-as.Date(df$baseline_date, format = "%d/%m/%Y")
  
  df$out_Angina_date <-as.Date(df$out_Angina_date, format = "%d/%m/%Y")
  df$out_CHD_date <-as.Date(df$out_CHD_date, format = "%d/%m/%Y")
  df$out_ARR_date <-as.Date(df$out_ARR_date, format = "%d/%m/%Y")
  df$out_HF_date <-as.Date(df$out_HF_date, format = "%d/%m/%Y")
  df$out_stroke_date <-as.Date(df$out_stroke_date, format = "%d/%m/%Y")
  df$out_PVD_date <-as.Date(df$out_PVD_date, format = "%d/%m/%Y")
  df$out_DVT_date <-as.Date(df$out_DVT_date, format = "%d/%m/%Y")
  df$out_s1_date <-as.Date(df$out_s1_date, format = "%d/%m/%Y")
  df$out_s2_date <-as.Date(df$out_s2_date, format = "%d/%m/%Y")
  df$out_s3_date <-as.Date(df$out_s3_date, format = "%d/%m/%Y")
  
  df$hist_Angina_date <-as.Date(df$hist_Angina_date, format = "%d/%m/%Y")
  df$hist_CHD_date <-as.Date(df$hist_CHD_date, format = "%d/%m/%Y")
  df$hist_ARR_date <-as.Date(df$hist_ARR_date, format = "%d/%m/%Y")
  df$hist_HF_date <-as.Date(df$hist_HF_date, format = "%d/%m/%Y")
  df$hist_stroke_date <-as.Date(df$hist_stroke_date, format = "%d/%m/%Y")
  df$hist_PVD_date <-as.Date(df$hist_PVD_date, format = "%d/%m/%Y")
  df$hist_DVT_date <-as.Date(df$hist_DVT_date, format = "%d/%m/%Y")
  df$hist_s1_date <-as.Date(df$hist_s1_date, format = "%d/%m/%Y")
  df$hist_s2_date <-as.Date(df$hist_s2_date, format = "%d/%m/%Y")
  df$hist_s3_date <-as.Date(df$hist_s3_date, format = "%d/%m/%Y")
  
  # Save dataset ready for Cox analysis
  saveRDS(df, paste0("~data/",df_name,"_forCox.rds"))
}
cox_preprocess_v3(df_v3,"cohort_2022_20251017")