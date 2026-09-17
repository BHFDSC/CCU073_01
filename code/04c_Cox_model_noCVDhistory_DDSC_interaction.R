######################################################################################
## TITLE: Cox models
##
## Description: Cox model of CVD outcomes with diabetes (type 1 and type 2)
##              stratified by sex and region
##              adjusted for
##              baseline_age, ethnicity,  deprivation, 
##              bmi category and smoking, number of COVID-19 vaccines
##              exposure: diabetes age (continuous) * diabetes history (yes/no)
##              + possible interaction with COVID-19
##
## For those without CVD history
##
## Author: Genevieve Cezard
######################################################################################
library(data.table)
library(survival)
library(tidyverse)
library(broom)

output_dir = "~output/"

#-----------#
# Load data #
#-----------#
# With diabetes type 2 duration and HbA1c new categories
df <- readRDS('~data/cohort_2022_20250509_duration_forCox.rds')
#df_pos1 <- readRDS('~data/cohort_2022_20250630_pos1_duration_forCox.rds')

#------------------------------------------------#
# Cox models with stratification and interaction #
#------------------------------------------------#
cox_analysis_nocvdhist_agecont_numexpo_interaction = function(df, i, inclusion_DMcat, exposure_var, ind_vars_cont, ind_vars_cat, interaction_term, strata_var, all_vars, output_name){
  
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
  
  # Count Number of events and Number of people for each category of covariates
  df_temp <- tmp_df %>%
    #  group_by_at(exposure_var) %>%
    summarise(
      total_n = n(),
      events = sum(out_flag == "1"),
      pyears = sum(time)/365.25,
      category=""
    )%>%
    ungroup() %>%
    #rename(category = exposure_var) %>%
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
  print(df_temp)

  # Set formula, run model and get results
  formula = as.formula(paste(c(paste0("Surv(tstart, tstop, out_flag) ~ ",exposure_var), ind_vars_cont, ind_vars_cat, interaction_term, paste0("strata(",strata_var,")")), collapse = " + "))
  #formula = as.formula("Surv(tstart, tstop, out_flag) ~ expo_T2DM_hist + expo_T2DM_hist * expo_T2DM_age10c50 + expo_T2DM_age10c50 + baseline_age + ethnicity + deprivation + strata(sex, region)
  print(formula)

  dmerge = tmerge(tmp_df,tmp_df,id=PERSON_ID,tstop=time)
  model <- coxph(formula,dmerge)
  print(model)
  
  print("Print covariance matrix - with cov()")
  cov = vcov(model)
  #print(cov)
  write.csv(cov, paste0(output_dir, "HRs_cov_matrix_",i,"_",output_name,".csv"), row.names = TRUE)

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
  #print(hr_results)
  
  # Format results with Ns for categorical variables
  df_temp3 <- left_join(df_temp,hr_results, by="term") %>% mutate(
    #HR = ifelse(variable %in% c("region"), NA, ifelse( is.na(HR), 1, HR)),
    #Results = ifelse(variable %in% c("region"), NA, ifelse(is.na(Results),"Reference",Results))
    HR = ifelse(variable %in% c("sex","region"), NA, ifelse( is.na(HR), 1, HR)),
    Results = ifelse(variable %in% c("sex","region"), NA, ifelse(is.na(Results),"Reference",Results))
    #HR = ifelse(variable %in% c("sex","region","base_age_cat"), NA, ifelse( is.na(HR), 1, HR)),
    #Results = ifelse(variable %in% c("sex","region","base_age_cat"), NA, ifelse(is.na(Results),"Reference",Results))
  )
  
  # Get results for continuous variables and interaction tern
  df_temp4 <- subset(right_join(df_temp,hr_results, by="term"), is.na(variable))
  #print(df_temp4)
  
  # Get general N
  df_n_model <- data.frame(term = "N_model", variable="",category="", 
                           N_event_rounded5 = round(sum(tmp_df$out_flag)/5.0)*5,
                           N_rounded5 = round(model$n/5.0)*5,
                           PY = round(as.numeric(gsub(" days","",sum(tmp_df$time)/365.25)), 0),
                           HR = NA, ucl = NA, lcl = NA, p_value = NA, Results = NA
  )
  
  #Combine
  df_results <- rbind(df_n_model,df_temp3,df_temp4)
  
  write.csv(df_results, paste0(output_dir, "HRs_",i,"_",output_name,".csv"), row.names = FALSE)
}

########
# T2DM #
########
# Define sets of adjustment without stratification variables (e.g. sex and region):
ind_vars_cont = c("baseline_age")

set1_ind_vars_cat_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid")
set2_ind_vars_cat_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "ethnicity", "deprivation")
set3_ind_vars_cat_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "ethnicity", "deprivation", "Smoking_status", "bmi_cat")
set4_ind_vars_cat_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "ethnicity", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num")
set5_ind_vars_cat_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "ethnicity", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num","hba1c_cat")

set1_ind_vars_cat_DMhist_BMI  = c("expo_T2DM_DDSC_hist", "bmi_cat")
set2_ind_vars_cat_DMhist_BMI  = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid")
set3_ind_vars_cat_DMhist_BMI  = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "ethnicity", "deprivation")
set4_ind_vars_cat_DMhist_BMI  = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "ethnicity", "deprivation", "Smoking_status")
set5_ind_vars_cat_DMhist_BMI  = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "ethnicity", "deprivation", "Smoking_status", "vaccination_num")
set6_ind_vars_cat_DMhist_BMI  = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "ethnicity", "deprivation", "Smoking_status", "vaccination_num","hba1c_cat")

set1_ind_vars_cat_DMhist_ethnicity  = c("expo_T2DM_DDSC_hist", "ethnicity")
set2_ind_vars_cat_DMhist_ethnicity  = c("expo_T2DM_DDSC_hist", "ethnicity", "covid")
set3_ind_vars_cat_DMhist_ethnicity  = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "deprivation")
set4_ind_vars_cat_DMhist_ethnicity  = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "deprivation", "Smoking_status", "bmi_cat")
set5_ind_vars_cat_DMhist_ethnicity  = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num")
set6_ind_vars_cat_DMhist_ethnicity  = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num","hba1c_cat")

set1_ind_vars_cat_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex")
set2_ind_vars_cat_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid")
set3_ind_vars_cat_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "deprivation", "ethnicity")
set4_ind_vars_cat_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "deprivation", "ethnicity", "Smoking_status", "bmi_cat")
set5_ind_vars_cat_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "deprivation", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num")
set6_ind_vars_cat_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "deprivation", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num","hba1c_cat")

set1_ind_vars_cat_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation")
set2_ind_vars_cat_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid")
set3_ind_vars_cat_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "ethnicity")
set4_ind_vars_cat_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "ethnicity", "Smoking_status", "bmi_cat")
set5_ind_vars_cat_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num")
set6_ind_vars_cat_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num","hba1c_cat")


# Define sets of variables to include in output file (i.e. adjustment and stratification variables):
set1_all_vars_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "sex", "region")
set2_all_vars_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation")
set3_all_vars_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "bmi_cat")
set4_all_vars_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "bmi_cat", "vaccination_num")
set5_all_vars_DMhist_covid = c("expo_T2DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "bmi_cat", "vaccination_num","hba1c_cat")

set1_all_vars_DMhist_BMI = c("expo_T2DM_DDSC_hist", "bmi_cat", "sex", "region")
set2_all_vars_DMhist_BMI = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "sex", "region")
set3_all_vars_DMhist_BMI = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "sex", "region", "ethnicity", "deprivation")
set4_all_vars_DMhist_BMI = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status")
set5_all_vars_DMhist_BMI = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "vaccination_num")
set6_all_vars_DMhist_BMI = c("expo_T2DM_DDSC_hist", "bmi_cat", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "vaccination_num","hba1c_cat")

set1_all_vars_DMhist_ethnicity = c("expo_T2DM_DDSC_hist", "ethnicity", "sex", "region")
set2_all_vars_DMhist_ethnicity = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "sex", "region")
set3_all_vars_DMhist_ethnicity = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "sex", "region", "deprivation")
set4_all_vars_DMhist_ethnicity = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "sex", "region", "deprivation","Smoking_status", "bmi_cat")
set5_all_vars_DMhist_ethnicity = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "sex", "region", "deprivation","Smoking_status", "bmi_cat","vaccination_num")
set6_all_vars_DMhist_ethnicity = c("expo_T2DM_DDSC_hist", "ethnicity", "covid", "sex", "region", "deprivation","Smoking_status", "bmi_cat","vaccination_num","hba1c_cat")

set1_all_vars_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "region")
set2_all_vars_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "region")
set3_all_vars_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "region", "deprivation", "ethnicity")
set4_all_vars_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "region", "deprivation", "ethnicity", "Smoking_status", "bmi_cat")
set5_all_vars_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "region", "deprivation", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num")
set6_all_vars_DMhist_sex  = c("expo_T2DM_DDSC_hist", "sex", "covid", "region", "deprivation", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num","hba1c_cat")

set1_all_vars_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "sex", "region")
set2_all_vars_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "sex", "region")
set3_all_vars_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "sex", "region", "ethnicity")
set4_all_vars_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "sex", "region", "ethnicity", "Smoking_status", "bmi_cat")
set5_all_vars_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "sex", "region", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num")
set6_all_vars_DMhist_deprivation  = c("expo_T2DM_DDSC_hist", "deprivation", "covid", "sex", "region", "ethnicity", "Smoking_status", "bmi_cat", "vaccination_num","hba1c_cat")


########
# T1DM #
########
# Define sets of adjustment without stratification variables (e.g. sex and region):
ind_vars_cont = c("baseline_age")

set1_ind_vars_cat_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid")
set2_ind_vars_cat_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "ethnicity", "deprivation")
set3_ind_vars_cat_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "ethnicity", "deprivation", "Smoking_status", "bmi_cat")
set4_ind_vars_cat_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "ethnicity", "deprivation", "Smoking_status", "bmi_cat","vaccination_num")
set5_ind_vars_cat_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "ethnicity", "deprivation", "Smoking_status", "bmi_cat", "vaccination_num",'hba1c_cat')

# Define sets of variables to include in output file (i.e. adjustment and stratification variables):
set1_all_vars_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "sex", "region")
set2_all_vars_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation")
set3_all_vars_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "bmi_cat")
set4_all_vars_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "bmi_cat","vaccination_num")
set5_all_vars_T1DMhist_covid = c("expo_T1DM_DDSC_hist", "covid", "sex", "region", "ethnicity", "deprivation","Smoking_status", "bmi_cat", "vaccination_num",'hba1c_cat')


# For specific outcome
outcomes = c('Angina','CHD','ARR','HF','stroke','PVD','DVT','s1','s2','s3')

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - run/rerun in May 2025
for (i in outcomes){
  print(i)
  # Interaction with Covid - rerun in May 2025
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
for (i in outcomes){
  print(i)
  # Interaction with BMI - run in May/June 2025
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_BMI,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*bmi_cat","sex, region",set1_all_vars_DMhist_BMI,"all_DDSC_DMage_DMhist_BMI_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_BMI,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*bmi_cat","sex, region",set2_all_vars_DMhist_BMI,"all_DDSC_DMage_DMhist_BMI_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_BMI,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*bmi_cat","sex, region",set3_all_vars_DMhist_BMI,"all_DDSC_DMage_DMhist_BMI_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_BMI,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*bmi_cat","sex, region",set4_all_vars_DMhist_BMI,"all_DDSC_DMage_DMhist_BMI_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_BMI,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*bmi_cat","sex, region",set5_all_vars_DMhist_BMI,"all_DDSC_DMage_DMhist_BMI_adj_set5")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set6_ind_vars_cat_DMhist_BMI,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*bmi_cat","sex, region",set6_all_vars_DMhist_BMI,"all_DDSC_DMage_DMhist_BMI_adj_set6")
  gc()
}
for (i in outcomes){
  print(i)
  # Interaction with Ethnicity - run in May/June 2025
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_ethnicity,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*ethnicity","sex, region",set1_all_vars_DMhist_ethnicity,"all_DDSC_DMage_DMhist_ethnicity_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_ethnicity,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*ethnicity","sex, region",set2_all_vars_DMhist_ethnicity,"all_DDSC_DMage_DMhist_ethnicity_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_ethnicity,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*ethnicity","sex, region",set3_all_vars_DMhist_ethnicity,"all_DDSC_DMage_DMhist_ethnicity_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_ethnicity,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*ethnicity","sex, region",set4_all_vars_DMhist_ethnicity,"all_DDSC_DMage_DMhist_ethnicity_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_ethnicity,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*ethnicity","sex, region",set5_all_vars_DMhist_ethnicity,"all_DDSC_DMage_DMhist_ethnicity_adj_set5")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set6_ind_vars_cat_DMhist_ethnicity,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*ethnicity","sex, region",set6_all_vars_DMhist_ethnicity,"all_DDSC_DMage_DMhist_ethnicity_adj_set6")
  gc()
}
for (i in outcomes){
  print(i)
  # Interaction with Deprivation - run in July 2025
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_deprivation,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*deprivation","sex, region",set1_all_vars_DMhist_deprivation,"all_DDSC_DMage_DMhist_deprivation_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_deprivation,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*deprivation","sex, region",set2_all_vars_DMhist_deprivation,"all_DDSC_DMage_DMhist_deprivation_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_deprivation,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*deprivation","sex, region",set3_all_vars_DMhist_deprivation,"all_DDSC_DMage_DMhist_deprivation_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_deprivation,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*deprivation","sex, region",set4_all_vars_DMhist_deprivation,"all_DDSC_DMage_DMhist_deprivation_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_deprivation,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*deprivation","sex, region",set5_all_vars_DMhist_deprivation,"all_DDSC_DMage_DMhist_deprivation_adj_set5")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set6_ind_vars_cat_DMhist_deprivation,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*deprivation","sex, region",set6_all_vars_DMhist_deprivation,"all_DDSC_DMage_DMhist_deprivation_adj_set6")
  gc()
}
for (i in outcomes){
  print(i)
  # Interaction with Sex - run in July 2025 -  check stratification formatting code in function before running -stratification by region only?
  #cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_sex,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*sex","region",set1_all_vars_DMhist_sex,"all_DDSC_DMage_DMhist_sex_adj_set1")
  #cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_sex,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*sex","region",set2_all_vars_DMhist_sex,"all_DDSC_DMage_DMhist_sex_adj_set2")
  #cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_sex,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*sex","region",set3_all_vars_DMhist_sex,"all_DDSC_DMage_DMhist_sex_adj_set3")
  #cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_sex,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*sex","region",set4_all_vars_DMhist_sex,"all_DDSC_DMage_DMhist_sex_adj_set4")
  #cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_sex,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*sex","region",set5_all_vars_DMhist_sex,"all_DDSC_DMage_DMhist_sex_adj_set5")
  #cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set6_ind_vars_cat_DMhist_sex,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*sex","region",set6_all_vars_DMhist_sex,"all_DDSC_DMage_DMhist_sex_adj_set6")
  #gc()
}

# Models with T2DM stratified by sex and region - for a different definition of the outcomes (1st position only in HES APC) - rerun in June 2025
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_pos1_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_pos1_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_pos1_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_pos1_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_pos1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_DDSC_DMage_DMhist_covid_pos1_adj_set5")
  gc()
}

# Models with T1DM stratified by sex and region - rerun in June 2025
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_T1DMhist_covid,"expo_T1DM_DDSC_hist*expo_T1DM_DDSC_age10c50*covid","sex, region",set1_all_vars_T1DMhist_covid,"all_DDSC_T1DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_T1DMhist_covid,"expo_T1DM_DDSC_hist*expo_T1DM_DDSC_age10c50*covid","sex, region",set2_all_vars_T1DMhist_covid,"all_DDSC_T1DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_T1DMhist_covid,"expo_T1DM_DDSC_hist*expo_T1DM_DDSC_age10c50*covid","sex, region",set3_all_vars_T1DMhist_covid,"all_DDSC_T1DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_T1DMhist_covid,"expo_T1DM_DDSC_hist*expo_T1DM_DDSC_age10c50*covid","sex, region",set4_all_vars_T1DMhist_covid,"all_DDSC_T1DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T1DM","No diabetes"), "expo_T1DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_T1DMhist_covid,"expo_T1DM_DDSC_hist*expo_T1DM_DDSC_age10c50*covid","sex, region",set5_all_vars_T1DMhist_covid,"all_DDSC_T1DMage_DMhist_covid_adj_set5")
  gc()
}

#----------------------------------------------#
# Models with duration of diabetes as exposure #
#----------------------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - run in June 2025
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "duration_T2DM_10y", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*duration_T2DM_10y*covid","sex, region",set1_all_vars_DMhist_covid,"all_T2DMduration_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "duration_T2DM_10y", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*duration_T2DM_10y*covid","sex, region",set2_all_vars_DMhist_covid,"all_T2DMduration_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "duration_T2DM_10y", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*duration_T2DM_10y*covid","sex, region",set3_all_vars_DMhist_covid,"all_T2DMduration_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "duration_T2DM_10y", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*duration_T2DM_10y*covid","sex, region",set4_all_vars_DMhist_covid,"all_T2DMduration_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "duration_T2DM_10y", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*duration_T2DM_10y*covid","sex, region",set5_all_vars_DMhist_covid,"all_T2DMduration_DMhist_covid_adj_set5")
  gc()
}  

#--------------------------#
# Subgroup analyses by sex #
#--------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Male - run in May 2025
df_Male <- df [df$sex == 'Male']
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Male_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Male_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Male_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Male_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Male, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Male_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Female - run in May 2025
df_Female <- df [df$sex == 'Female']
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Female_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Female_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Female_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Female_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Female, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Female_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

#--------------------------------#
# Subgroup analyses by ethnicity #
#--------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - White - run in May 2025
df_White <- df [df$ethnicity == 'White']
print('White')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_White_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_White_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_White_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_White_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_White, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_White_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Asian - run in May 2025
df_Asian <- df [df$ethnicity == 'Asian']
print('Asian')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Asian, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Asian_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Asian, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Asian_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Asian, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Asian_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Asian, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Asian_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Asian, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Asian_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Black - run in May 2025
df_Black <- df [df$ethnicity == 'Black']
print('Black')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Black_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Black_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Black_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Black_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Black, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Black_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Mixed - run in May 2025
df_Mixed <- df [df$ethnicity == 'Mixed']
print('Mixed')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Mixed_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Mixed_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Mixed_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Mixed_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Mixed, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Mixed_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Other - run in May 2025
df_Other <- df [df$ethnicity == 'Other']
print('Other')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Other_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Other_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Other_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Other_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Other, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Other_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

#----------------------------------#
# Subgroup analyses by deprivation #
#----------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 5 Least deprived - run in June 2025
df_Deprivation5 <- df [df$deprivation == '9_10_least_deprived']
print('Deprivation5 - Least deprived')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation5, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 4 - run in June 2025
df_Deprivation4 <- df [df$deprivation == '7_8']
print('Deprivation4')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation4, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 3 - run in June 2025
df_Deprivation3 <- df [df$deprivation == '5_6']
print('Deprivation3')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation3, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 2 - run in June 2025
df_Deprivation2 <- df [df$deprivation == '3_4']
print('Deprivation2')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation2, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Quintile 1 - run in June 2025
df_Deprivation1 <- df [df$deprivation == '1_2_most_deprived']
print('Deprivation1 - Most deprived')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Deprivation1, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

#-----------------------------------#
# Subgroup analyses by BMI category #
#-----------------------------------#
# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Underweight - run in November 2025
df_Underweight <- df [df$bmi_cat == 'Underweight']
print('Underweight')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Underweight_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Underweight_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Underweight_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Underweight_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Underweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Underweight_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Healthy weight - run in November 2025
df_Healthy <- df [df$bmi_cat == 'Healthy']
print('Healthy')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Healthy_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Healthy_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Healthy_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Healthy_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Healthy, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Healthy_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Overweight -run in November 2025
df_Overweight <- df [df$bmi_cat == 'Overweight']
print('Overweight')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Overweight_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Overweight_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Overweight_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Overweight_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Overweight, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Overweight_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Overweight - run in November 2025
df_Obese <- df [df$bmi_cat == 'Obese']
print('Obese')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_Obese_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_Obese_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_Obese_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_Obese_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_Obese, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_Obese_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

# Models with T2DM stratified by sex and region - main definition of the outcomes (all positions in HES APC) - Missing BMI - run in November 2025
df_MissingBMI <- df [df$bmi_cat == 'Missing']
print('Missing BMI')
for (i in outcomes){
  print(i)
  # Interaction with Covid
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set1_all_vars_DMhist_covid,"all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set2_all_vars_DMhist_covid,"all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set3_all_vars_DMhist_covid,"all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set4_all_vars_DMhist_covid,"all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df_MissingBMI, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region",set5_all_vars_DMhist_covid,"all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set5")
  gc()
}

#-------------------------------------------------------------#
# Sensitivity analyses with additional age cat stratification #
#-------------------------------------------------------------#
set1_all_vars_DMhist_covid_agecatstrata = c("expo_T2DM_DDSC_hist", "covid", "sex", "region","base_age_cat")
set2_all_vars_DMhist_covid_agecatstrata = c("expo_T2DM_DDSC_hist", "covid", "sex", "region","base_age_cat", "ethnicity", "deprivation")
set3_all_vars_DMhist_covid_agecatstrata = c("expo_T2DM_DDSC_hist", "covid", "sex", "region","base_age_cat", "ethnicity", "deprivation","Smoking_status", "bmi_cat")
set4_all_vars_DMhist_covid_agecatstrata = c("expo_T2DM_DDSC_hist", "covid", "sex", "region","base_age_cat", "ethnicity", "deprivation","Smoking_status", "bmi_cat", "vaccination_num")
set5_all_vars_DMhist_covid_agecatstrata = c("expo_T2DM_DDSC_hist", "covid", "sex", "region","base_age_cat", "ethnicity", "deprivation","Smoking_status", "bmi_cat", "vaccination_num", "hba1c_cat")

for (i in outcomes){
  print(i)
  # Interaction with Covid - run in June 2025
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set1_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region, base_age_cat",set1_all_vars_DMhist_covid_agecatstrata,"all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set1")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set2_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region, base_age_cat",set2_all_vars_DMhist_covid_agecatstrata,"all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set2")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set3_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region, base_age_cat",set3_all_vars_DMhist_covid_agecatstrata,"all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set3")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set4_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region, base_age_cat",set4_all_vars_DMhist_covid_agecatstrata,"all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set4")
  cox_analysis_nocvdhist_agecont_numexpo_interaction(df, i, c("T2DM","No diabetes"), "expo_T2DM_DDSC_age10c50", ind_vars_cont, set5_ind_vars_cat_DMhist_covid,"expo_T2DM_DDSC_hist*expo_T2DM_DDSC_age10c50*covid","sex, region, base_age_cat",set5_all_vars_DMhist_covid_agecatstrata,"all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set5")
  gc()
}


#--------------------------------------------#
# 3. Combine results for the 10 CVD outcomes #
#--------------------------------------------#
# Combine Cox HRs results
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

df_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_adj_set1")
df_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_adj_set2")
df_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_adj_set3")
df_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_adj_set4")
df_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_adj_set5")

df_all_adj_set1_int_T2DMage_BMI <- combine_cox_results("all_DDSC_DMage_DMhist_BMI_adj_set1")
df_all_adj_set2_int_T2DMage_BMI <- combine_cox_results("all_DDSC_DMage_DMhist_BMI_adj_set2")
df_all_adj_set3_int_T2DMage_BMI <- combine_cox_results("all_DDSC_DMage_DMhist_BMI_adj_set3")
df_all_adj_set4_int_T2DMage_BMI <- combine_cox_results("all_DDSC_DMage_DMhist_BMI_adj_set4")
df_all_adj_set5_int_T2DMage_BMI <- combine_cox_results("all_DDSC_DMage_DMhist_BMI_adj_set5")
df_all_adj_set6_int_T2DMage_BMI <- combine_cox_results("all_DDSC_DMage_DMhist_BMI_adj_set6")

df_all_adj_set1_int_T2DMage_ethnicity <- combine_cox_results("all_DDSC_DMage_DMhist_ethnicity_adj_set1")
df_all_adj_set2_int_T2DMage_ethnicity <- combine_cox_results("all_DDSC_DMage_DMhist_ethnicity_adj_set2")
df_all_adj_set3_int_T2DMage_ethnicity <- combine_cox_results("all_DDSC_DMage_DMhist_ethnicity_adj_set3")
df_all_adj_set4_int_T2DMage_ethnicity <- combine_cox_results("all_DDSC_DMage_DMhist_ethnicity_adj_set4")
df_all_adj_set5_int_T2DMage_ethnicity <- combine_cox_results("all_DDSC_DMage_DMhist_ethnicity_adj_set5")
df_all_adj_set6_int_T2DMage_ethnicity <- combine_cox_results("all_DDSC_DMage_DMhist_ethnicity_adj_set6")

df_all_adj_set1_int_T2DMage_deprivation <- combine_cox_results("all_DDSC_DMage_DMhist_deprivation_adj_set1")
df_all_adj_set2_int_T2DMage_deprivation <- combine_cox_results("all_DDSC_DMage_DMhist_deprivation_adj_set2")
df_all_adj_set3_int_T2DMage_deprivation <- combine_cox_results("all_DDSC_DMage_DMhist_deprivation_adj_set3")
df_all_adj_set4_int_T2DMage_deprivation <- combine_cox_results("all_DDSC_DMage_DMhist_deprivation_adj_set4")
df_all_adj_set5_int_T2DMage_deprivation <- combine_cox_results("all_DDSC_DMage_DMhist_deprivation_adj_set5")
df_all_adj_set6_int_T2DMage_deprivation <- combine_cox_results("all_DDSC_DMage_DMhist_deprivation_adj_set6")

df_all_adj_set1_int_T2DMage_sex <- combine_cox_results("all_DDSC_DMage_DMhist_sex_adj_set1")
df_all_adj_set2_int_T2DMage_sex <- combine_cox_results("all_DDSC_DMage_DMhist_sex_adj_set2")
df_all_adj_set3_int_T2DMage_sex <- combine_cox_results("all_DDSC_DMage_DMhist_sex_adj_set3")
df_all_adj_set4_int_T2DMage_sex <- combine_cox_results("all_DDSC_DMage_DMhist_sex_adj_set4")
df_all_adj_set5_int_T2DMage_sex <- combine_cox_results("all_DDSC_DMage_DMhist_sex_adj_set5")
df_all_adj_set6_int_T2DMage_sex <- combine_cox_results("all_DDSC_DMage_DMhist_sex_adj_set6")

df_all_adj_set1_int_T2DMage_covid_pos1 <- combine_cox_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set1")
df_all_adj_set2_int_T2DMage_covid_pos1 <- combine_cox_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set2")
df_all_adj_set3_int_T2DMage_covid_pos1 <- combine_cox_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set3")
df_all_adj_set4_int_T2DMage_covid_pos1 <- combine_cox_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set4")
df_all_adj_set5_int_T2DMage_covid_pos1 <- combine_cox_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set5")

df_all_adj_set1_int_T1DMage_covid <- combine_cox_results("all_DDSC_T1DMage_DMhist_covid_adj_set1")
df_all_adj_set2_int_T1DMage_covid <- combine_cox_results("all_DDSC_T1DMage_DMhist_covid_adj_set2")
df_all_adj_set3_int_T1DMage_covid <- combine_cox_results("all_DDSC_T1DMage_DMhist_covid_adj_set3")
df_all_adj_set4_int_T1DMage_covid <- combine_cox_results("all_DDSC_T1DMage_DMhist_covid_adj_set4")
df_all_adj_set5_int_T1DMage_covid <- combine_cox_results("all_DDSC_T1DMage_DMhist_covid_adj_set5")

#+exploration by sex and 2-term interaction
df_Male_adj_set1_int_T2DMage <- combine_cox_results("all_Male_DDSC_DMage_DMhist_adj_set1")
df_Male_adj_set5_int_T2DMage <- combine_cox_results("all_Male_DDSC_DMage_DMhist_adj_set5")
df_Female_adj_set1_int_T2DMage <- combine_cox_results("all_Female_DDSC_DMage_DMhist_adj_set1")
df_Female_adj_set5_int_T2DMage <- combine_cox_results("all_Female_DDSC_DMage_DMhist_adj_set5")

# Analyses using duration of diabetes as exposure
df_all_adj_set1_int_T2DMduration_covid <- combine_cox_results("all_T2DMduration_DMhist_covid_adj_set1")
df_all_adj_set2_int_T2DMduration_covid <- combine_cox_results("all_T2DMduration_DMhist_covid_adj_set2")
df_all_adj_set3_int_T2DMduration_covid <- combine_cox_results("all_T2DMduration_DMhist_covid_adj_set3")
df_all_adj_set4_int_T2DMduration_covid <- combine_cox_results("all_T2DMduration_DMhist_covid_adj_set4")
df_all_adj_set5_int_T2DMduration_covid <- combine_cox_results("all_T2DMduration_DMhist_covid_adj_set5")

# By sex
df_Male_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Male_DDSC_DMage_DMhist_covid_adj_set1")
df_Male_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Male_DDSC_DMage_DMhist_covid_adj_set2")
df_Male_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Male_DDSC_DMage_DMhist_covid_adj_set3")
df_Male_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Male_DDSC_DMage_DMhist_covid_adj_set4")
df_Male_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Male_DDSC_DMage_DMhist_covid_adj_set5")

df_Female_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Female_DDSC_DMage_DMhist_covid_adj_set1")
df_Female_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Female_DDSC_DMage_DMhist_covid_adj_set2")
df_Female_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Female_DDSC_DMage_DMhist_covid_adj_set3")
df_Female_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Female_DDSC_DMage_DMhist_covid_adj_set4")
df_Female_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Female_DDSC_DMage_DMhist_covid_adj_set5")

# By ethnicity
df_White_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_White_DDSC_DMage_DMhist_covid_adj_set1")
df_White_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_White_DDSC_DMage_DMhist_covid_adj_set2")
df_White_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_White_DDSC_DMage_DMhist_covid_adj_set3")
df_White_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_White_DDSC_DMage_DMhist_covid_adj_set4")
df_White_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_White_DDSC_DMage_DMhist_covid_adj_set5")

df_Asian_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set1")
df_Asian_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set2")
df_Asian_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set3")
df_Asian_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set4")
df_Asian_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set5")

df_Black_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Black_DDSC_DMage_DMhist_covid_adj_set1")
df_Black_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Black_DDSC_DMage_DMhist_covid_adj_set2")
df_Black_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Black_DDSC_DMage_DMhist_covid_adj_set3")
df_Black_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Black_DDSC_DMage_DMhist_covid_adj_set4")
df_Black_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Black_DDSC_DMage_DMhist_covid_adj_set5")

df_Mixed_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set1")
df_Mixed_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set2")
df_Mixed_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set3")
df_Mixed_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set4")
df_Mixed_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set5")

df_Other_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Other_DDSC_DMage_DMhist_covid_adj_set1")
df_Other_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Other_DDSC_DMage_DMhist_covid_adj_set2")
df_Other_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Other_DDSC_DMage_DMhist_covid_adj_set3")
df_Other_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Other_DDSC_DMage_DMhist_covid_adj_set4")
df_Other_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Other_DDSC_DMage_DMhist_covid_adj_set5")

# By deprivation
df_Deprivation5_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set1")
df_Deprivation5_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set2")
df_Deprivation5_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set3")
df_Deprivation5_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set4")
df_Deprivation5_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set5")

df_Deprivation4_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set1")
df_Deprivation4_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set2")
df_Deprivation4_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set3")
df_Deprivation4_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set4")
df_Deprivation4_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set5")

df_Deprivation3_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set1")
df_Deprivation3_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set2")
df_Deprivation3_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set3")
df_Deprivation3_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set4")
df_Deprivation3_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set5")

df_Deprivation2_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set1")
df_Deprivation2_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set2")
df_Deprivation2_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set3")
df_Deprivation2_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set4")
df_Deprivation2_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set5")

df_Deprivation1_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set1")
df_Deprivation1_all_adj_set2_int_T2DMage_covid <- combine_cox_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set2")
df_Deprivation1_all_adj_set3_int_T2DMage_covid <- combine_cox_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set3")
df_Deprivation1_all_adj_set4_int_T2DMage_covid <- combine_cox_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set4")
df_Deprivation1_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set5")

# By BMI category
df_Underweight_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Underweight_DDSC_DMage_DMhist_covid_adj_set1")
df_Underweight_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Underweight_DDSC_DMage_DMhist_covid_adj_set5")

df_Healthy_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Healthy_DDSC_DMage_DMhist_covid_adj_set1")
df_Healthy_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Healthy_DDSC_DMage_DMhist_covid_adj_set5")

df_Overweight_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Overweight_DDSC_DMage_DMhist_covid_adj_set1")
df_Overweight_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Overweight_DDSC_DMage_DMhist_covid_adj_set5")

df_Obese_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_Obese_DDSC_DMage_DMhist_covid_adj_set1")
df_Obese_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_Obese_DDSC_DMage_DMhist_covid_adj_set5")

df_MissingBMI_all_adj_set1_int_T2DMage_covid <- combine_cox_results("all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set1")
df_MissingBMI_all_adj_set5_int_T2DMage_covid <- combine_cox_results("all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set5")


# For age cat stratification
df_all_agecatstrata_adj_set1_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set1")
df_all_agecatstrata_adj_set2_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set2")
df_all_agecatstrata_adj_set3_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set3")
df_all_agecatstrata_adj_set4_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set4")
df_all_agecatstrata_adj_set5_int_T2DMage_covid <- combine_cox_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set5")

# Combine Cox Covariance matrix results
combine_cox_cov_results = function(output_name){
  df_results = data.frame()
  for (i in outcomes){
    #print(i)
    df_tmp = read.csv(paste0(output_dir, "HRs_cov_matrix_",i,"_",output_name,".csv"))
    df_tmp$outcome <- paste0(i)
    df_results <-rbind(df_results,df_tmp)
  }
  #print(df_results)
  write.csv(df_results, paste0(output_dir, "HRs_cov_matrix_",output_name,".csv"), row.names = FALSE)
  return(df_results)
}

df_covmatrix_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_adj_set1_int_T2DMage_BMI <- combine_cox_cov_results("all_DDSC_DMage_DMhist_BMI_adj_set1")
df_covmatrix_adj_set2_int_T2DMage_BMI <- combine_cox_cov_results("all_DDSC_DMage_DMhist_BMI_adj_set2")
df_covmatrix_adj_set3_int_T2DMage_BMI <- combine_cox_cov_results("all_DDSC_DMage_DMhist_BMI_adj_set3")
df_covmatrix_adj_set4_int_T2DMage_BMI <- combine_cox_cov_results("all_DDSC_DMage_DMhist_BMI_adj_set4")
df_covmatrix_adj_set5_int_T2DMage_BMI <- combine_cox_cov_results("all_DDSC_DMage_DMhist_BMI_adj_set5")
df_covmatrix_adj_set6_int_T2DMage_BMI <- combine_cox_cov_results("all_DDSC_DMage_DMhist_BMI_adj_set6")

df_covmatrix_adj_set1_int_T2DMage_ethnicity <- combine_cox_cov_results("all_DDSC_DMage_DMhist_ethnicity_adj_set1")
df_covmatrix_adj_set2_int_T2DMage_ethnicity <- combine_cox_cov_results("all_DDSC_DMage_DMhist_ethnicity_adj_set2")
df_covmatrix_adj_set3_int_T2DMage_ethnicity <- combine_cox_cov_results("all_DDSC_DMage_DMhist_ethnicity_adj_set3")
df_covmatrix_adj_set4_int_T2DMage_ethnicity <- combine_cox_cov_results("all_DDSC_DMage_DMhist_ethnicity_adj_set4")
df_covmatrix_adj_set5_int_T2DMage_ethnicity <- combine_cox_cov_results("all_DDSC_DMage_DMhist_ethnicity_adj_set5")
df_covmatrix_adj_set6_int_T2DMage_ethnicity <- combine_cox_cov_results("all_DDSC_DMage_DMhist_ethnicity_adj_set6")

df_covmatrix_adj_set1_int_T2DMage_deprivation <- combine_cox_cov_results("all_DDSC_DMage_DMhist_deprivation_adj_set1")
df_covmatrix_adj_set2_int_T2DMage_deprivation <- combine_cox_cov_results("all_DDSC_DMage_DMhist_deprivation_adj_set2")
df_covmatrix_adj_set3_int_T2DMage_deprivation <- combine_cox_cov_results("all_DDSC_DMage_DMhist_deprivation_adj_set3")
df_covmatrix_adj_set4_int_T2DMage_deprivation <- combine_cox_cov_results("all_DDSC_DMage_DMhist_deprivation_adj_set4")
df_covmatrix_adj_set5_int_T2DMage_deprivation <- combine_cox_cov_results("all_DDSC_DMage_DMhist_deprivation_adj_set5")
df_covmatrix_adj_set6_int_T2DMage_deprivation <- combine_cox_cov_results("all_DDSC_DMage_DMhist_deprivation_adj_set6")

df_covmatrix_adj_set1_int_T2DMage_sex <- combine_cox_cov_results("all_DDSC_DMage_DMhist_sex_adj_set1")
df_covmatrix_adj_set2_int_T2DMage_sex <- combine_cox_cov_results("all_DDSC_DMage_DMhist_sex_adj_set2")
df_covmatrix_adj_set3_int_T2DMage_sex <- combine_cox_cov_results("all_DDSC_DMage_DMhist_sex_adj_set3")
df_covmatrix_adj_set4_int_T2DMage_sex <- combine_cox_cov_results("all_DDSC_DMage_DMhist_sex_adj_set4")
df_covmatrix_adj_set5_int_T2DMage_sex <- combine_cox_cov_results("all_DDSC_DMage_DMhist_sex_adj_set5")
df_covmatrix_adj_set6_int_T2DMage_sex <- combine_cox_cov_results("all_DDSC_DMage_DMhist_sex_adj_set6")

df_covmatrix_adj_set1_int_T2DMage_covid_pos1 <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set1")
df_covmatrix_adj_set2_int_T2DMage_covid_pos1 <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set2")
df_covmatrix_adj_set3_int_T2DMage_covid_pos1 <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set3")
df_covmatrix_adj_set4_int_T2DMage_covid_pos1 <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set4")
df_covmatrix_adj_set5_int_T2DMage_covid_pos1 <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_pos1_adj_set5")

df_covmatrix_adj_set1_int_T1DMage_covid <- combine_cox_cov_results("all_DDSC_T1DMage_DMhist_covid_adj_set1")
df_covmatrix_adj_set2_int_T1DMage_covid <- combine_cox_cov_results("all_DDSC_T1DMage_DMhist_covid_adj_set2")
df_covmatrix_adj_set3_int_T1DMage_covid <- combine_cox_cov_results("all_DDSC_T1DMage_DMhist_covid_adj_set3")
df_covmatrix_adj_set4_int_T1DMage_covid <- combine_cox_cov_results("all_DDSC_T1DMage_DMhist_covid_adj_set4")
df_covmatrix_adj_set5_int_T1DMage_covid <- combine_cox_cov_results("all_DDSC_T1DMage_DMhist_covid_adj_set5")

# Analyses using duration of diabetes as exposure
df_covmatrix_all_adj_set1_int_T2DMduration_covid <- combine_cox_cov_results("all_T2DMduration_DMhist_covid_adj_set1")
df_covmatrix_all_adj_set2_int_T2DMduration_covid <- combine_cox_cov_results("all_T2DMduration_DMhist_covid_adj_set2")
df_covmatrix_all_adj_set3_int_T2DMduration_covid <- combine_cox_cov_results("all_T2DMduration_DMhist_covid_adj_set3")
df_covmatrix_all_adj_set4_int_T2DMduration_covid <- combine_cox_cov_results("all_T2DMduration_DMhist_covid_adj_set4")
df_covmatrix_all_adj_set5_int_T2DMduration_covid <- combine_cox_cov_results("all_T2DMduration_DMhist_covid_adj_set5")

# By sex
df_covmatrix_Male_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Male_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Male_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Male_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Male_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Male_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Male_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Male_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Male_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Male_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Female_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Female_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Female_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Female_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Female_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Female_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Female_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Female_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Female_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Female_DDSC_DMage_DMhist_covid_adj_set5")

# By ethnicity
df_covmatrix_White_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_White_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_White_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_White_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_White_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_White_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_White_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_White_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_White_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_White_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Asian_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Asian_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Asian_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Asian_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Asian_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Asian_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Black_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Black_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Black_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Black_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Black_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Black_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Black_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Black_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Black_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Black_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Mixed_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Mixed_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Mixed_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Mixed_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Mixed_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Mixed_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Other_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Other_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Other_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Other_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Other_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Other_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Other_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Other_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Other_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Other_DDSC_DMage_DMhist_covid_adj_set5")

# By deprivation
df_covmatrix_Deprivation5_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Deprivation5_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Deprivation5_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Deprivation5_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Deprivation5_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation5_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Deprivation4_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Deprivation4_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Deprivation4_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Deprivation4_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Deprivation4_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation4_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Deprivation3_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Deprivation3_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Deprivation3_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Deprivation3_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Deprivation3_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation3_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Deprivation2_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Deprivation2_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Deprivation2_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Deprivation2_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Deprivation2_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation2_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Deprivation1_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Deprivation1_all_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set2")
df_covmatrix_Deprivation1_all_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set3")
df_covmatrix_Deprivation1_all_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set4")
df_covmatrix_Deprivation1_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Deprivation1_DDSC_DMage_DMhist_covid_adj_set5")


# By BMI category
df_covmatrix_Underweight_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Underweight_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Underweight_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Underweight_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Healthy_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Healthy_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Healthy_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Healthy_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Overweight_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Overweight_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Overweight_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Overweight_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_Obese_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_Obese_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_Obese_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_Obese_DDSC_DMage_DMhist_covid_adj_set5")

df_covmatrix_MissingBMI_all_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set1")
df_covmatrix_MissingBMI_all_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_MissingBMI_DDSC_DMage_DMhist_covid_adj_set5")


# For age cat stratification
df_covmatrix_all_agecatstrata_adj_set1_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set1")
df_covmatrix_all_agecatstrata_adj_set2_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set2")
df_covmatrix_all_agecatstrata_adj_set3_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set3")
df_covmatrix_all_agecatstrata_adj_set4_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set4")
df_covmatrix_all_agecatstrata_adj_set5_int_T2DMage_covid <- combine_cox_cov_results("all_DDSC_DMage_DMhist_covid_agecatstrata_adj_set5")


#--------------------------------------------------------------------------------------------------------------#
# Test code to merge all models of a specific type with all 10 outcomes and all sets of adjustment into 1 file # 
#--------------------------------------------------------------------------------------------------------------#
combine_cox_results_v2 = function(output_name){
  df_results = data.frame()
  for (i in outcomes){
    for (j in set){
      df_tmp = read.csv(paste0(output_dir, "HRs_",i,"_",output_name,"_adj_",j,".csv"))
      df_tmp <- df_tmp %>% mutate(
        N_event_rounded5 = ifelse(N_event_rounded5 <= 10, 10, N_event_rounded5),
        N_rounded5 = ifelse(N_rounded5 <= 10, 10, N_rounded5)
      )
      df_tmp$outcome <- paste0(i)
      df_tmp$adjustment <- paste0(j)
      df_results <-rbind(df_results,df_tmp)
      #}
    }
  }
  #print(df_results)
  write.csv(df_results, paste0(output_dir, "HRs_",output_name,".csv"), row.names = FALSE)
  return(df_results)
}
set = c("set1","set2","set3","set4")
#set = c("set0","set1","set2","set3")
test <- combine_cox_results_v2("all_DDSC_DMage_DMhist_covid") # => ~ 1430 rows so too long

set = c("set1","set2","set5")
df_all_int_T2DMage_BMI <- combine_cox_results_v2("all_DDSC_DMage_DMhist_BMI") # => ~ 1250 rows so too long

