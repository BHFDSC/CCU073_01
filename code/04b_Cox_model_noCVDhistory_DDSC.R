######################################################################################
## TITLE: Cox models
##
## Description: Cox model of CVD outcomes with diabetes (type 1 and type 2)
##              stratified by sex and region
##              adjusted for
##              baseline_age, ethnicity,  deprivation, 
##              bmi category and smoking, number of COVID-19 vaccines
##
## For those without CVD history
##
## Author: Genevieve Cezard
#############################################################################################
library(data.table)
library(survival)
library(tidyverse)
library(broom)

output_dir = "~output/"

#--------------#
# 1. Load data #
#--------------#
# With outcomes defined using main diagnosis position only from hospitalisation records
df_pos1 <- readRDS('~data/cohort_2022_20250630_pos1_duration_forCox.rds')
df_pos1_covid <- df_pos1 [df_pos1$covid == 'Yes']
df_pos1_non_covid <- df_pos1 [df_pos1$covid == 'No']

# With diabetes type 2 duration and HbA1c new categories
df <- readRDS('~data/cohort_2022_20250509_duration_forCox.rds')
df_covid <- df [df$covid == 'Yes']
df_non_covid <- df [df$covid == 'No']

# With extra exposure variables: "age at T2DM diagnosis - sex"
df <- readRDS('~data/cohort_2022_20251017_forCox.rds')


#-----------------------------------------------------------------------------------------#
# 2. Cox models for full cohort and COVID-19/non-COVID-19 sub-cohorts with stratification #
#-----------------------------------------------------------------------------------------#
cox_analysis_nocvdhist_agecont = function(df, i, inclusion_DMcat, exposure_var, ind_vars_cont, ind_vars_cat, strata_var, all_vars, output_name){

  tmp_df <- df
  
  # Apply inclusion/exclusion 
  # 1. Exclude those with history of the outcome
  tmp_df[,hist_flag:=0]
  tmp_df[pmin(baseline_date,get(paste0('hist_s1_date')))<=baseline_date, hist_flag:=1]
  print(paste0("Prevalent events: ", nrow(tmp_df[tmp_df$hist_flag == 1])))
  tmp_df <- tmp_df[tmp_df$hist_flag == 0]
  # 2. Exclude the other type of diabetes
  print(paste0("Dataset - Number of rows before exclusion of other type of diabetes : ",nrow(tmp_df)))
  print(table(tmp_df$diabetes_status_DDSC))
  tmp_df <- tmp_df[tmp_df$diabetes_status_DDSC %in% inclusion_DMcat]
  print(paste0("Dataset - Number of rows after exclusion of other type of diabetes : ",nrow(tmp_df)))

  # Define outcome flag and associated time
  tmp_df[,out_flag:=0]
  tmp_df[pmin(DOD,data_end,get(paste0('out_',i,'_date')),na.rm = TRUE)==get(paste0('out_',i,'_date')),out_flag:=1]
  tmp_df[out_flag==1,time:=get(paste0('out_',i,'_date'))-baseline_date]
  tmp_df[out_flag==0,time:=pmin(DOD,data_end,na.rm = TRUE)-baseline_date]
  tmp_df[,time:=time+1]
  print(paste0("Number of events: ", nrow(tmp_df[tmp_df$out_flag == 1])))
  
  # reduce dataset to complete cases based on exposure and covariates
  tmp_df <- subset(tmp_df,!is.na(exposure_var))
  #print(nrow(tmp_df)) 
  for (j in all_vars){
    tmp_df <- subset(tmp_df,!is.na(get(j)))
  }
  print(nrow(tmp_df))
  
  # Count Number of events and Number of people for each category of exposure and covariates
  df_temp <- tmp_df %>%
    group_by_at(all_of(exposure_var)) %>%
    summarise(
      total_n = n(),
      events = sum(out_flag == "1"),
      pyears = sum(time)/365.25
    )%>%
    ungroup() %>%
    rename(category = exposure_var) %>%
    mutate(variable = exposure_var)
  
  for (j in all_vars){
    df_temp2 <- tmp_df %>%
      group_by_at(vars(one_of(j))) %>%
      summarise(
        total_n = n(),
        events = sum(out_flag == "1"),
        pyears = sum(time)/365.25
      )%>%
      ungroup() %>%
      rename(category = one_of(j)) %>%
      mutate(variable = j)
    df_temp <- rbind(df_temp, df_temp2)
  }
  df_temp <- df_temp %>% 
    mutate(
      N_rounded5 = round(total_n/5.0)*5,
      N_event_rounded5 = round(events/5.0)*5,
      PY = round(as.numeric(gsub(" days","",pyears)), 0),
      term = paste0(variable,category)
    )%>%
    select(term, variable, category,N_event_rounded5,N_rounded5,PY)
  #return(df_temp)
  
  # Set formula, run model and get results
  formula = as.formula(paste(c(paste0("Surv(tstart, tstop, out_flag) ~ ",exposure_var), ind_vars_cont, ind_vars_cat, paste0("strata(",strata_var,")")), collapse = " + "))
  print(formula)
  
  dmerge = tmerge(tmp_df,tmp_df,id=PERSON_ID,tstop=time)
  model <- coxph(formula,dmerge)
  print(model)
  
  #print("Print covariance matrix - with cov()")
  #cov = vcov(model)
  #print(cov)
  #write.csv(cov, paste0(output_dir, "HRs_cov_matrix_",i,"_",output_name,".csv"), row.names = TRUE)

  hr_results <- tidy(model) %>%
    mutate(
      HR = round(exp(estimate), 2),
      ucl =  round(exp(estimate + 1.96 * std.error), 2),
      lcl =  round(exp(estimate - 1.96 * std.error), 2),
      p_value = round(p.value, 3)
    ) %>%
    select(term, HR, lcl, ucl, p_value) %>%
    mutate(Results = paste0(HR, " (", lcl, ", ", ucl, ")")) %>%
    data.frame()
  
  df_temp3 <- left_join(df_temp,hr_results, by="term")%>%
    mutate(
      HR = ifelse(is.na(HR),1,HR),
      Results = ifelse(is.na(Results),"Reference",Results)
    )
  #print(df_temp3)
  
  # Get results for continuous variables
  df_temp4 <- subset(right_join(df_temp,hr_results, by="term"), is.na(variable) & !is.na(HR))
  #print(df_temp4)
  
  df_n_model <- data.frame(term = "N_model", variable="",category="", 
                           N_event_rounded5 = round(sum(tmp_df$out_flag)/5.0)*5,
                           N_rounded5 = round(model$n/5.0)*5,
                           PY = round(as.numeric(gsub(" days","",sum(tmp_df$time)/365.25)), 0),
                           HR = NA, ucl = NA, lcl = NA, p_value = NA, Results = NA
  )
  df_results <- rbind(df_n_model,df_temp3,df_temp4)
  
  write.csv(df_results, paste0(output_dir, "HRs_",i,"_",output_name,".csv"), row.names = FALSE)
}


# Define sets of adjustment without stratification variables (e.g. sex and region):
ind_vars_cont = c("baseline_age")

set1_ind_vars_cat = c("")
set2_ind_vars_cat = c("ethnicity", "deprivation")
set3_ind_vars_cat = c("ethnicity", "deprivation", "Smoking_status", "bmi_cat")
set4_ind_vars_cat = c("ethnicity", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num")
set5_ind_vars_cat = c("ethnicity", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num",'hba1c_cat')

set6_ind_vars_cat = c("ethnicity", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num",'hba1c_flag')

# Define sets of variables to include in output file (i.e. adjustment and stratification variables):
set1_all_vars = c("sex", "region")
set2_all_vars = c("sex", "region", "ethnicity", "deprivation")
set3_all_vars = c("sex", "region", "ethnicity", "deprivation","Smoking_status","bmi_cat")
set4_all_vars = c("sex", "ethnicity", "region", "deprivation","Smoking_status","bmi_cat","vaccination_num")
set5_all_vars = c("sex", "ethnicity", "region", "deprivation","Smoking_status","bmi_cat","vaccination_num",'hba1c_cat')

set6_all_vars = c("sex", "ethnicity", "region", "deprivation","Smoking_status","bmi_cat","vaccination_num",'hba1c_flag')

# For specific outcome
outcomes = c('Angina','CHD','ARR','HF','stroke','PVD','DVT','s1','s2','s3')

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC)
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
  
  # Adjusted for set 6
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set6_ind_vars_cat,"sex, region",set5_all_vars,"all_expo_DDSC_T2DM_noCVDhist_adj_set6")
  gc()
} 

# Models with T2DM stratified by sex and region - for a different definition of the outcomes (1st position only in HES APC)
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_expo_DDSC_T2DM_pos1_adj_set1")
  cox_analysis_nocvdhist_agecont(df_pos1_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_expo_DDSC_T2DM_pos1_adj_set1")
  cox_analysis_nocvdhist_agecont(df_pos1_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_expo_DDSC_T2DM_pos1_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_expo_DDSC_T2DM_pos1_adj_set2")
  cox_analysis_nocvdhist_agecont(df_pos1_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_expo_DDSC_T2DM_pos1_adj_set2")
  cox_analysis_nocvdhist_agecont(df_pos1_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_expo_DDSC_T2DM_pos1_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_expo_DDSC_T2DM_pos1_adj_set3")
  cox_analysis_nocvdhist_agecont(df_pos1_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_expo_DDSC_T2DM_pos1_adj_set3")
  cox_analysis_nocvdhist_agecont(df_pos1_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_expo_DDSC_T2DM_pos1_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_expo_DDSC_T2DM_pos1_adj_set4")
  cox_analysis_nocvdhist_agecont(df_pos1_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_expo_DDSC_T2DM_pos1_adj_set4")
  cox_analysis_nocvdhist_agecont(df_pos1_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_expo_DDSC_T2DM_pos1_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_expo_DDSC_T2DM_pos1_adj_set5")
  cox_analysis_nocvdhist_agecont(df_pos1_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_expo_DDSC_T2DM_pos1_adj_set5")
  cox_analysis_nocvdhist_agecont(df_pos1_non_covid, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_expo_DDSC_T2DM_pos1_adj_set5")
  gc()
} 

# Models with T1DM stratified by sex and region - with a different categorisation of age at diagnosis - all 50+ combined
for (i in outcomes){
  print(i)

  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_expo_DDSC_T1DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_expo_DDSC_T1DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_expo_DDSC_T1DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_expo_DDSC_T1DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_expo_DDSC_T1DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_expo_DDSC_T1DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_expo_DDSC_T1DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_expo_DDSC_T1DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_expo_DDSC_T1DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_expo_DDSC_T1DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_expo_DDSC_T1DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_expo_DDSC_T1DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_expo_DDSC_T1DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_expo_DDSC_T1DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_10y_age_cat_50plus", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_expo_DDSC_T1DM_noCVDhist_adj_set5")
  gc()
}

#----------------------------------------------#
# Models with duration of diabetes as exposure #
#----------------------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC)
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_expo_T2DM_duration_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_expo_T2DM_duration_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_expo_T2DM_duration_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "duration_T2DM_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_expo_T2DM_duration_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_expo_T2DM_duration_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_expo_T2DM_duration_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_expo_T2DM_duration_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_expo_T2DM_duration_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_expo_T2DM_duration_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_expo_T2DM_duration_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_expo_T2DM_duration_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_expo_T2DM_duration_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_expo_T2DM_duration_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_expo_T2DM_duration_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid, i, c("T2DM","No diabetes"), "duration_T2DM_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_expo_T2DM_duration_noCVDhist_adj_set5")
  gc()
} 

#--------------------------------------------------------------#
# Models with age at T2DM diagnosis and sex groups as exposure #
#--------------------------------------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC)
set1_all_vars = c("region") # Without sex
set5_all_vars = c("ethnicity", "region", "deprivation","Smoking_status","bmi_cat","vaccination_num","hba1c_cat") # Without sex
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat_sex", ind_vars_cont, set1_ind_vars_cat,"region",set1_all_vars,"all_expo_T2DM_age_sex_noCVDhist_adj_set1")

  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat_sex", ind_vars_cont, set5_ind_vars_cat,"region",set5_all_vars,"all_expo_T2DM_age_sex_noCVDhist_adj_set5")
  gc()
} 

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - TEST STRATIFIED BY SEX
set1_all_vars = c("sex", "region")
set5_all_vars = c("sex", "region", "ethnicity","deprivation","Smoking_status","bmi_cat","vaccination_num","hba1c_cat")
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat_sex", ind_vars_cont, set1_ind_vars_cat,"sex,region",set1_all_vars,"all_expo_T2DM_age_sex_noCVDhist_stratasex_adj_set1")
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat_sex", ind_vars_cont, set5_ind_vars_cat,"sex,region",set5_all_vars,"all_expo_T2DM_age_sex_noCVDhist_stratasex_adj_set5")
  gc()
} 

#--------------------------#
# Subgroup analyses by sex #
#--------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Male
df_Male <- df [df$sex == 'Male']
df_covid_Male <- df_covid [df_covid$sex == 'Male']
df_non_covid_Male <- df_non_covid [df_non_covid$sex == 'Male']

print('Male')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"region",set1_all_vars,"all_Male_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"region",set1_all_vars,"covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"region",set1_all_vars,"non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"region",set2_all_vars,"all_Male_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"region",set2_all_vars,"covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"region",set2_all_vars,"non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"region",set3_all_vars,"all_Male_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"region",set3_all_vars,"covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"region",set3_all_vars,"non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"region",set4_all_vars,"all_Male_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"region",set4_all_vars,"covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"region",set4_all_vars,"non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"region",set5_all_vars,"all_Male_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"region",set5_all_vars,"covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"region",set5_all_vars,"non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Female
df_Female <- df [df$sex == 'Female']
df_covid_Female <- df_covid [df_covid$sex == 'Female']
df_non_covid_Female <- df_non_covid [df_non_covid$sex == 'Female']

print('Female')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"region",set1_all_vars,"all_Female_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"region",set1_all_vars,"covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"region",set1_all_vars,"non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"region",set2_all_vars,"all_Female_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"region",set2_all_vars,"covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"region",set2_all_vars,"non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"region",set3_all_vars,"all_Female_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"region",set3_all_vars,"covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"region",set3_all_vars,"non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"region",set4_all_vars,"all_Female_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"region",set4_all_vars,"covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"region",set4_all_vars,"non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Female_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

#----------------------------------------------#
# Subgroup analyses by ethnicity (EXPLORATION) #
#----------------------------------------------#
#outcomes = c('Angina')
#outcomes = c('s1')

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - White
df_White <- df [df$sex == 'White']
df_covid_White <- df_covid [df_covid$sex == 'White']
df_non_covid_White <- df_non_covid [df_non_covid$sex == 'White']

print('White')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_white_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_white_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_white_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_white_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_white_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_white_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_white_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_white_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_white_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_white_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_white_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_white_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
}


# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Asian
df_Asian <- df [df$ethnicity == 'Asian']
df_covid_Asian <- df_covid [df_covid$ethnicity == 'Asian']
df_non_covid_Asian <- df_non_covid [df_non_covid$ethnicity == 'Asian']

print('Asian')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_SA_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_SA_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_SA_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_SA_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_SA, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_SA_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Black
df_Black <- df [df$ethnicity == 'Black']
df_covid_Black <- df_covid [df_covid$ethnicity == 'Black']
df_non_covid_Black <- df_non_covid [df_non_covid$ethnicity == 'Black']

print('Black')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Black_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Black_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Black_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Black_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Black_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Mixed
df_Mixed <- df [df$ethnicity == 'Mixed']
df_covid_Mixed <- df_covid [df_covid$ethnicity == 'Mixed']
df_non_covid_Mixed <- df_non_covid [df_non_covid$ethnicity == 'Mixed']

print('Mixed')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Mixed_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Other
df_Other <- df [df$ethnicity == 'Other']
df_covid_Other <- df_covid [df_covid$ethnicity == 'Other']
df_non_covid_Other <- df_non_covid [df_non_covid$ethnicity == 'Other']

print('Other')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Other_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Other_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Other_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Other_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Other_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
}

#----------------------------------#
# Subgroup analyses by deprivation #
#----------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 5 Least deprived - run in August 2025
df_Deprivation5 <- df [df$deprivation == '9_10_least_deprived']
df_covid_Deprivation5 <- df_covid [df_covid$deprivation == '9_10_least_deprived']
df_non_covid_Deprivation5 <- df_non_covid [df_non_covid$deprivation == '9_10_least_deprived']

print('Deprivation5 - Least deprived')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 4 - run in August 2025
df_Deprivation4 <- df [df$deprivation == '7_8']
df_covid_Deprivation4 <- df_covid [df_covid$deprivation == '7_8']
df_non_covid_Deprivation4 <- df_non_covid [df_non_covid$deprivation == '7_8']

print('Deprivation4')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 3 - run in August 2025
df_Deprivation3 <- df [df$deprivation == '5_6']
df_covid_Deprivation3 <- df_covid [df_covid$deprivation == '5_6']
df_non_covid_Deprivation3 <- df_non_covid [df_non_covid$deprivation == '5_6']

print('Deprivation3')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 2 - run in August 2025
df_Deprivation2 <- df [df$deprivation == '3_4']
df_covid_Deprivation2 <- df_covid [df_covid$deprivation == '3_4']
df_non_covid_Deprivation2 <- df_non_covid [df_non_covid$deprivation == '3_4']

print('Deprivation2')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 1 Most deprived - run in August 2025
df_Deprivation1 <- df [df$deprivation == '1_2_most_deprived']
df_covid_Deprivation1 <- df_covid [df_covid$deprivation == '1_2_most_deprived']
df_non_covid_Deprivation1 <- df_non_covid [df_non_covid$deprivation == '1_2_most_deprived']

print('Deprivation1 - Most deprived')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

#-----------------------------------#
# Subgroup analyses by BMI category #
#-----------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Underweight - run in August 2025
df_Underweight <- df [df$bmi_cat == 'Underweight']
df_covid_Underweight <- df_covid [df_covid$bmi_cat == 'Underweight']
df_non_covid_Underweight <- df_non_covid [df_non_covid$bmi_cat == 'Underweight']

print('Underweight')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Healthy - run in August 2025
df_Healthy <- df [df$bmi_cat == 'Healthy']
df_covid_Healthy <- df_covid [df_covid$bmi_cat == 'Healthy']
df_non_covid_Healthy <- df_non_covid [df_non_covid$bmi_cat == 'Healthy']

print('Healthy')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Overweight - run in August 2025
df_Overweight <- df [df$bmi_cat == 'Overweight']
df_covid_Overweight <- df_covid [df_covid$bmi_cat == 'Overweight']
df_non_covid_Overweight <- df_non_covid [df_non_covid$bmi_cat == 'Overweight']

print('Overweight')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Overweightexpo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Overweight - run in August 2025
df_Obese <- df [df$bmi_cat == 'Obese']
df_covid_Obese <- df_covid [df_covid$bmi_cat == 'Obese']
df_non_covid_Obese <- df_non_covid [df_non_covid$bmi_cat == 'Obese']

print('Obese')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_Obese_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_Obese_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_Obese_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_Obese_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_Obese_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_Obese_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Missing BMI - run in August 2025
df_MissingBMI <- df [df$bmi_cat == 'Missing']
df_covid_MissingBMI <- df_covid [df_covid$bmi_cat == 'Missing']
df_non_covid_MissingBMI <- df_non_covid [df_non_covid$bmi_cat == 'Missing']

print('Missing BMI')
for (i in outcomes){
  print(i)
  
  # Adjusted for set 1
  cox_analysis_nocvdhist_agecont(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"all_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set1")
  cox_analysis_nocvdhist_agecont(df_non_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set1_ind_vars_cat,"sex, region",set1_all_vars,"non_covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set1")
  gc()
  
  # Adjusted for set 2
  cox_analysis_nocvdhist_agecont(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat",ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"all_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set2")
  cox_analysis_nocvdhist_agecont(df_non_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set2_ind_vars_cat,"sex, region",set2_all_vars,"non_covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set2")
  gc()
  
  # Adjusted for set 3
  cox_analysis_nocvdhist_agecont(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"all_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set3")
  cox_analysis_nocvdhist_agecont(df_non_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set3_ind_vars_cat,"sex, region",set3_all_vars,"non_covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set3")
  gc()
  
  # Adjusted for set 4
  cox_analysis_nocvdhist_agecont(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"all_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set4")
  cox_analysis_nocvdhist_agecont(df_non_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set4_ind_vars_cat,"sex, region",set4_all_vars,"non_covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set4")
  gc()
  
  # Adjusted for set 5
  cox_analysis_nocvdhist_agecont(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"all_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set5")
  cox_analysis_nocvdhist_agecont(df_non_covid_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_10y_age_cat", ind_vars_cont, set5_ind_vars_cat,"sex, region",set5_all_vars,"non_covid_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set5")
  gc()
}

#--------------------------------------------#
# 3. Combine results for the 10 CVD outcomes #
#--------------------------------------------#
combine_cox_results = function(output_name){
  df_results = data.frame()
  for (i in outcomes){
    #print(i)
    df_tmp = read.csv(paste0(output_dir, "HRs_",i,"_",output_name,".csv"))
    df_tmp <- df_tmp %>% mutate(
      N_event_rounded5 = ifelse(N_event_rounded5 <= 10, 10, N_event_rounded5),
      N_rounded5 = ifelse(N_rounded5 <= 10, 10, N_rounded5)
    )
    df_tmp$outcome <- paste0(i)
    df_results <-rbind(df_results,df_tmp)
  }
  #print(df_results)
  write.csv(df_results, paste0(output_dir, "HRs_",output_name,".csv"), row.names = FALSE)
  return(df_results)
}


## T2DM age cat ##
df_all_T2DMagecat_adj_set1 <- combine_cox_results("all_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_covid_T2DMagecat_adj_set1 <- combine_cox_results("covid_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_non_covid_T2DMagecat_adj_set1 <- combine_cox_results("non_covid_expo_DDSC_T2DM_noCVDhist_adj_set1")

df_all_T2DMagecat_adj_set2 <- combine_cox_results("all_expo_DDSC_T2DM_noCVDhist_adj_set2")
df_covid_T2DMagecat_adj_set2 <- combine_cox_results("covid_expo_DDSC_T2DM_noCVDhist_adj_set2")
df_non_covid_T2DMagecat_adj_set2 <- combine_cox_results("non_covid_expo_DDSC_T2DM_noCVDhist_adj_set2")

df_all_T2DMagecat_adj_set3 <- combine_cox_results("all_expo_DDSC_T2DM_noCVDhist_adj_set3")
df_covid_T2DMagecat_adj_set3 <- combine_cox_results("covid_expo_DDSC_T2DM_noCVDhist_adj_set3")
df_non_covid_T2DMagecat_adj_set3 <- combine_cox_results("non_covid_expo_DDSC_T2DM_noCVDhist_adj_set3")

df_all_T2DMagecat_adj_set4 <- combine_cox_results("all_expo_DDSC_T2DM_noCVDhist_adj_set4")
df_covid_T2DMagecat_adj_set4 <- combine_cox_results("covid_expo_DDSC_T2DM_noCVDhist_adj_set4")
df_non_covid_T2DMagecat_adj_set4 <- combine_cox_results("non_covid_expo_DDSC_T2DM_noCVDhist_adj_set4")

df_all_T2DMagecat_adj_set5 <- combine_cox_results("all_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_covid_T2DMagecat_adj_set5 <- combine_cox_results("covid_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_non_covid_T2DMagecat_adj_set5 <- combine_cox_results("non_covid_expo_DDSC_T2DM_noCVDhist_adj_set5")

df_all_T2DMagecat_adj_set6 <- combine_cox_results("all_expo_DDSC_T2DM_noCVDhist_adj_set6")

## T2DM age cat (outcome with 1st position only in HES APC) ##
df_all_T2DMagecat_pos1_adj_set1 <- combine_cox_results("all_expo_DDSC_T2DM_pos1_adj_set1")
df_covid_T2DMagecat_pos1_adj_set1 <- combine_cox_results("covid_expo_DDSC_T2DM_pos1_adj_set1")
df_non_covid_T2DMagecat_pos1_adj_set1 <- combine_cox_results("non_covid_expo_DDSC_T2DM_pos1_adj_set1")

df_all_T2DMagecat_pos1_adj_set2 <- combine_cox_results("all_expo_DDSC_T2DM_pos1_adj_set2")
df_covid_T2DMagecat_pos1_adj_set2 <- combine_cox_results("covid_expo_DDSC_T2DM_pos1_adj_set2")
df_non_covid_T2DMagecat_pos1_adj_set2 <- combine_cox_results("non_covid_expo_DDSC_T2DM_pos1_adj_set2")

df_all_T2DMagecat_pos1_adj_set3 <- combine_cox_results("all_expo_DDSC_T2DM_pos1_adj_set3")
df_covid_T2DMagecat_pos1_adj_set3 <- combine_cox_results("covid_expo_DDSC_T2DM_pos1_adj_set3")
df_non_covid_T2DMagecat_pos1_adj_set3 <- combine_cox_results("non_covid_expo_DDSC_T2DM_pos1_adj_set3")

df_all_T2DMagecat_pos1_adj_set4 <- combine_cox_results("all_expo_DDSC_T2DM_pos1_adj_set4")
df_covid_T2DMagecat_pos1_adj_set4 <- combine_cox_results("covid_expo_DDSC_T2DM_pos1_adj_set4")
df_non_covid_T2DMagecat_pos1_adj_set4 <- combine_cox_results("non_covid_expo_DDSC_T2DM_pos1_adj_set4")

df_all_T2DMagecat_pos1_adj_set5 <- combine_cox_results("all_expo_DDSC_T2DM_pos1_adj_set5")
df_covid_T2DMagecat_pos1_adj_set5 <- combine_cox_results("covid_expo_DDSC_T2DM_pos1_adj_set5")
df_non_covid_T2DMagecat_pos1_adj_set5 <- combine_cox_results("non_covid_expo_DDSC_T2DM_pos1_adj_set5")

## T1DM age cat (with 10 years age band for 0 to 50 and then 50+) ##
df_all_T1DMagecat_adj_set1 <- combine_cox_results("all_expo_DDSC_T1DM_noCVDhist_adj_set1")
df_covid_T1DMagecat_adj_set1 <- combine_cox_results("covid_expo_DDSC_T1DM_noCVDhist_adj_set1")
df_non_covid_T1DMagecat_adj_set1 <- combine_cox_results("non_covid_expo_DDSC_T1DM_noCVDhist_adj_set1")

df_all_T1DMagecat_adj_set2 <- combine_cox_results("all_expo_DDSC_T1DM_noCVDhist_adj_set2")
df_covid_T1DMagecat_adj_set2 <- combine_cox_results("covid_expo_DDSC_T1DM_noCVDhist_adj_set2")
df_non_covid_T1DMagecat_adj_set2 <- combine_cox_results("non_covid_expo_DDSC_T1DM_noCVDhist_adj_set2")

df_all_T1DMagecat_adj_set3 <- combine_cox_results("all_expo_DDSC_T1DM_noCVDhist_adj_set3")
df_covid_T1DMagecat_adj_set3 <- combine_cox_results("covid_expo_DDSC_T1DM_noCVDhist_adj_set3")
df_non_covid_T1DMagecat_adj_set3 <- combine_cox_results("non_covid_expo_DDSC_T1DM_noCVDhist_adj_set3")

df_all_T1DMagecat_adj_set4 <- combine_cox_results("all_expo_DDSC_T1DM_noCVDhist_adj_set4")
df_covid_T1DMagecat_adj_set4 <- combine_cox_results("covid_expo_DDSC_T1DM_noCVDhist_adj_set4")
df_non_covid_T1DMagecat_adj_set4 <- combine_cox_results("non_covid_expo_DDSC_T1DM_noCVDhist_adj_set4")

df_all_T1DMagecat_adj_set5 <- combine_cox_results("all_expo_DDSC_T1DM_noCVDhist_adj_set5")
df_covid_T1DMagecat_adj_set5 <- combine_cox_results("covid_expo_DDSC_T1DM_noCVDhist_adj_set5")
df_non_covid_T1DMagecat_adj_set5 <- combine_cox_results("non_covid_expo_DDSC_T1DM_noCVDhist_adj_set5")


## T2DM duration cat ##
df_all_T2DMdurationcat_adj_set1 <- combine_cox_results("all_expo_T2DM_duration_noCVDhist_adj_set1")
df_covid_T2DMdurationcat_adj_set1 <- combine_cox_results("covid_expo_T2DM_duration_noCVDhist_adj_set1")
df_non_covid_T2DMdurationcat_adj_set1 <- combine_cox_results("non_covid_expo_T2DM_duration_noCVDhist_adj_set1")

df_all_T2DMdurationcat_adj_set2 <- combine_cox_results("all_expo_T2DM_duration_noCVDhist_adj_set2")
df_covid_T2DMdurationcat_adj_set2 <- combine_cox_results("covid_expo_T2DM_duration_noCVDhist_adj_set2")
df_non_covid_T2DMdurationcat_adj_set2 <- combine_cox_results("non_covid_expo_T2DM_duration_noCVDhist_adj_set2")

df_all_T2DMdurationcat_adj_set3 <- combine_cox_results("all_expo_T2DM_duration_noCVDhist_adj_set3")
df_covid_T2DMdurationcat_adj_set3 <- combine_cox_results("covid_expo_T2DM_duration_noCVDhist_adj_set3")
df_non_covid_T2DMdurationcat_adj_set3 <- combine_cox_results("non_covid_expo_T2DM_duration_noCVDhist_adj_set3")

df_all_T2DMdurationcat_adj_set4 <- combine_cox_results("all_expo_T2DM_duration_noCVDhist_adj_set4")
df_covid_T2DMdurationcat_adj_set4 <- combine_cox_results("covid_expo_T2DM_duration_noCVDhist_adj_set4")
df_non_covid_T2DMdurationcat_adj_set4 <- combine_cox_results("non_covid_expo_T2DM_duration_noCVDhist_adj_set4")

df_all_T2DMdurationcat_adj_set5 <- combine_cox_results("all_expo_T2DM_duration_noCVDhist_adj_set5")
df_covid_T2DMdurationcat_adj_set5 <- combine_cox_results("covid_expo_T2DM_duration_noCVDhist_adj_set5")
df_non_covid_T2DMdurationcat_adj_set5 <- combine_cox_results("non_covid_expo_T2DM_duration_noCVDhist_adj_set5")


# T2DM age cat - Males ##
df_Male_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Male_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_Male_covid_T2DMagecat_adj_set1 <- combine_cox_results("covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_Male_non_covid_T2DMagecat_adj_set1 <- combine_cox_results("non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set1")

df_Male_all_T2DMagecat_adj_set2 <- combine_cox_results("all_Male_expo_DDSC_T2DM_noCVDhist_adj_set2")
df_Male_covid_T2DMagecat_adj_set2 <- combine_cox_results("covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set2")
df_Male_non_covid_T2DMagecat_adj_set2 <- combine_cox_results("non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set2")

df_Male_all_T2DMagecat_adj_set3 <- combine_cox_results("all_Male_expo_DDSC_T2DM_noCVDhist_adj_set3")
df_Male_covid_T2DMagecat_adj_set3 <- combine_cox_results("covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set3")
df_Male_non_covid_T2DMagecat_adj_set3 <- combine_cox_results("non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set3")

df_Male_all_T2DMagecat_adj_set4 <- combine_cox_results("all_Male_expo_DDSC_T2DM_noCVDhist_adj_set4")
df_Male_covid_T2DMagecat_adj_set4 <- combine_cox_results("covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set4")
df_Male_non_covid_T2DMagecat_adj_set4 <- combine_cox_results("non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set4")

df_Male_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Male_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_Male_covid_T2DMagecat_adj_set5 <- combine_cox_results("covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_Male_non_covid_T2DMagecat_adj_set5 <- combine_cox_results("non_covid_Male_expo_DDSC_T2DM_noCVDhist_adj_set5")


## T2DM cat - Females ##
df_Female_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Female_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_Female_covid_T2DMagecat_adj_set1 <- combine_cox_results("covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_Female_non_covid_T2DMagecat_adj_set1 <- combine_cox_results("non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set1")

df_Female_all_T2DMagecat_adj_set2 <- combine_cox_results("all_Female_expo_DDSC_T2DM_noCVDhist_adj_set2")
df_Female_covid_T2DMagecat_adj_set2 <- combine_cox_results("covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set2")
df_Female_non_covid_T2DMagecat_adj_set2 <- combine_cox_results("non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set2")

df_Female_all_T2DMagecat_adj_set3 <- combine_cox_results("all_Female_expo_DDSC_T2DM_noCVDhist_adj_set3")
df_Female_covid_T2DMagecat_adj_set3 <- combine_cox_results("covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set3")
df_Female_non_covid_T2DMagecat_adj_set3 <- combine_cox_results("non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set3")

df_Female_all_T2DMagecat_adj_set4 <- combine_cox_results("all_Female_expo_DDSC_T2DM_noCVDhist_adj_set4")
df_Female_covid_T2DMagecat_adj_set4 <- combine_cox_results("covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set4")
df_Female_non_covid_T2DMagecat_adj_set4 <- combine_cox_results("non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set4")

df_Female_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Female_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_Female_covid_T2DMagecat_adj_set5 <- combine_cox_results("covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_Female_non_covid_T2DMagecat_adj_set5 <- combine_cox_results("non_covid_Female_expo_DDSC_T2DM_noCVDhist_adj_set5")

## Combined T2DMcat/sex
df_all_T2DMagecat_sex_adj_set1 <- combine_cox_results("all_expo_T2DM_age_sex_noCVDhist_adj_set1")
df_all_T2DMagecat_sex_adj_set5 <- combine_cox_results("all_expo_T2DM_age_sex_noCVDhist_adj_set5")

## Combined T2DMcat/sex - with sex strata
df_all_T2DMagecat_sex_stratasex_adj_set1 <- combine_cox_results("all_expo_T2DM_age_sex_noCVDhist_stratasex_adj_set1")
df_all_T2DMagecat_sex_stratasex_adj_set5 <- combine_cox_results("all_expo_T2DM_age_sex_noCVDhist_stratasex_adj_set5")


## T2DM cat - Deprivation ##
df_deprivation5_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_deprivation4_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_deprivation3_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_deprivation2_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_deprivation1_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set1")

df_deprivation5_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Deprivation5_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_deprivation4_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Deprivation4_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_deprivation3_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Deprivation3_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_deprivation2_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Deprivation2_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_deprivation1_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Deprivation1_expo_DDSC_T2DM_noCVDhist_adj_set5")

## T2DM cat - BMI category ##
df_Underweight_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_Healthy_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_Overweight_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_Obese_all_T2DMagecat_adj_set1 <- combine_cox_results("all_Obese_expo_DDSC_T2DM_noCVDhist_adj_set1")
df_MissingBMI_all_T2DMagecat_adj_set1 <- combine_cox_results("all_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set1")

df_Underweight_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Underweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_Healthy_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Healthy_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_Overweight_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Overweight_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_Obese_all_T2DMagecat_adj_set5 <- combine_cox_results("all_Obese_expo_DDSC_T2DM_noCVDhist_adj_set5")
df_MissingBMI_all_T2DMagecat_adj_set5 <- combine_cox_results("all_MissingBMI_expo_DDSC_T2DM_noCVDhist_adj_set5")
