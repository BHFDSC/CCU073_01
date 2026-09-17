######################################################################################
## TITLE: Create table 1 for CCU073 project
##
## Author: Genevieve Cezard, Wen Shi
######################################################################################
library(data.table)
library(dplyr)

df = readRDS('~data/cohort_2022_20250509_duration.rds')

# Exclude people with prevalent s1
df = df[hist_s1_flag==0]
# Create covid and non_covid cohort datasets
df_covid <- df [df$covid == 'Yes']
df_non_covid <- df [df$covid == 'No']

#-------------------------------------------------#
# 1. Get descriptive statistics with real numbers #
#-------------------------------------------------#
table_1 <- function(df){
  n_total = df[,.N]
  tmpc = data.table()
  for (i in c('sex','base_age_cat','ethnicity','region','deprivation',
              'bmi_cat','egfr_calc_cat','glucose_cat','hba1c_cat','Smoking_status','vaccination_num',
              'diabetes_status_DDSC')){
    levs  = levels(df[,get(i)])
    for (j in levs){
      n = df[j,.N,on=i]
      percentage = round(n/n_total*100,1)
      tmp = data.table(Characteristic=i,Group=j,N=paste0(n,' (',percentage,'%)'))  
      tmpc = rbind(tmpc,tmp)
    }
  }
  for (i in c('age_T1DM_DDSC_10y_cat','age_T2DM_DDSC_10y_cat')){
    levs  = levels(df[,get(i)])
    for (j in levs){
      n = df[j,.N,on=i]
      tmp = data.table(Characteristic=i,Group=j,N=n)  
      tmpc = rbind(tmpc,tmp) 
      }
  }
  
  for (i in c('baseline_age','age_T1DM_DDSC','age_T2DM_DDSC','hba1c_value')){
    quants = df[!is.na(get(i)),quantile(get(i))]
    q1 = round(quants[2],1)
    q2 = round(quants[3],1)
    q3 = round(quants[4],1)
    n_na = df[is.na(get(i)),.N]
    percentage = round(n_na/n_total*100,1)
    tmp = data.table(Characteristic=rep(i,2),Group=c('Median(1st,3rd)','Missing'),
                     N=c(paste0(q2,' (',q1,',',q3,')'),paste0(n_na,' (',percentage,'%)')))
    tmpc = rbind(tmpc,tmp)
  }
  
  for (i in c('bmi','tchol','hdl','Glucose_level','sbp','Creatinine','egfr_calculated','dbp')){
    quants = df[!is.na(get(paste0('cov_',i,'_value'))),quantile(get(paste0('cov_',i,'_value')))]
    q1 = round(quants[2],1)
    q2 = round(quants[3],1)
    q3 = round(quants[4],1)
    n_na = df[is.na(get(paste0('cov_',i,'_date'))),.N]
    percentage = round(n_na/n_total*100,1)
    tmp = data.table(Characteristic=rep(i,2),Group=c('Median(1st,3rd)','Missing'),
                     N=c(paste0(q2,' (',q1,',',q3,')'),paste0(n_na,' (',percentage,'%)')))
    tmpc = rbind(tmpc,tmp)
  }
  return(tmpc)
}
tmpc <- table_1(df)

tmpc_covid <- table_1(df_covid)
names(tmpc_covid)[names(tmpc_covid) == 'N'] <- 'N_Covid_people'

tmpc_non_covid <- table_1(df_non_covid)
names(tmpc_non_covid)[names(tmpc_non_covid) == 'N'] <- 'N_non_Covid_people'

# merge() doesn't keep the characteristics in the intended order
#tmpc_1 <- merge(tmpc,tmpc_covid, by = c("Characteristic", "Group"))
#tmpc_final <- merge(tmpc_1,tmpc_non_covid, by = c("Characteristic", "Group"))

# left_join() keeps characteristics presented in the expected order
tmpc_1 <- left_join(tmpc,tmpc_covid, by = c("Characteristic", "Group"))
tmpc_final <- left_join(tmpc_1,tmpc_non_covid, by = c("Characteristic", "Group"))
  
fwrite(tmpc_final,'~output/table1.csv')


#----------------------------------------------------#
# 2. Get descriptive statistics with rounded numbers #
#----------------------------------------------------#
# Rounded numbers to the nearest 5
table_1_rounded <- function(df){
  n_total = df[,.N]
  n_total_rounded5 = round(n_total/5.0)*5
  tmpc = data.table(Characteristic="Total",Group="",N=n_total_rounded5)
  for (i in c('sex','base_age_cat','ethnicity','deprivation','region',
              'Smoking_status','bmi_cat','hba1c_cat',
              'egfr_calc_cat','glucose_cat','vaccination_num',
              'diabetes_status_DDSC','age_T1DM_DDSC_10y_cat','age_T2DM_DDSC_10y_cat'
              )){
    levs  = levels(df[,get(i)])
    for (j in levs){
      n = df[j,.N,on=i]
      n_rounded5 = round(n/5.0)*5
      
      percentage = round(n_rounded5/n_total_rounded5*100,1)
      
      tmp = data.table(Characteristic=i,Group=j,N=paste0(n_rounded5,' (',percentage,'%)'))  
      tmpc = rbind(tmpc,tmp)
    }
  }
  
  for (i in c('baseline_age','age_T1DM_DDSC','age_T2DM_DDSC','hba1c_value')){
    quants = df[!is.na(get(i)),quantile(get(i))]
    q1 = round(quants[2],1)
    q2 = round(quants[3],1)
    q3 = round(quants[4],1)
    n_na = df[is.na(get(i)),.N]
    n_na_rounded5 = round(n_na/5.0)*5
    tmp = data.table(Characteristic=rep(i,2),Group=c('Median(1st,3rd)','Missing'),
                     N=c(paste0(q2,' (',q1,',',q3,')'),n_na_rounded5))
    tmpc = rbind(tmpc,tmp)
  }
  for (i in c('bmi','tchol','hdl','Glucose_level','sbp','Creatinine','egfr_calculated','dbp')){
    quants = df[!is.na(get(paste0('cov_',i,'_value'))),quantile(get(paste0('cov_',i,'_value')))]
    q1 = round(quants[2],1)
    q2 = round(quants[3],1)
    q3 = round(quants[4],1)
    n_na = df[is.na(get(paste0('cov_',i,'_date'))),.N]
    n_na_rounded5 = round(n_na/5.0)*5
    tmp = data.table(Characteristic=rep(i,2),Group=c('Median(1st,3rd)','Missing'),
                     N=c(paste0(q2,' (',q1,',',q3,')'),n_na_rounded5))
    tmpc = rbind(tmpc,tmp)
  }
  return(tmpc)
}
tmpc_rounded <- table_1_rounded(df)

tmpc_covid_rounded <- table_1_rounded(df_covid)
names(tmpc_covid_rounded)[names(tmpc_covid_rounded) == 'N'] <- 'N_Covid_people'

tmpc_non_covid_rounded <- table_1_rounded(df_non_covid)
names(tmpc_non_covid_rounded)[names(tmpc_non_covid_rounded) == 'N'] <- 'N_non_Covid_people'

tmpc_1_rounded <- left_join(tmpc_rounded,tmpc_covid_rounded, by = c("Characteristic", "Group"))
tmpc_final_rounded <- left_join(tmpc_1_rounded,tmpc_non_covid_rounded, by = c("Characteristic", "Group"))

# fwrite(tmpc_final_rounded,'~output/table1_rounded.csv')
tmpc_final_rounded_export=tmpc_final_rounded[!seq(48,54,2)]
fwrite(tmpc_final_rounded_export,'~output/table1_rounded_exclS1.csv')


#--------------------------------------#
# Get age at diagnosis by ethnic group #
#--------------------------------------#
table_1_by_var_rounded <- function(df,var){
  
  #tmpc <- data.table()
  tmpc = data.table(Group=var, N="", Value.age_T1DM_DDSC = 'Median (Q1,Q3)', Value.age_T2DM_DDSC = 'Median(Q1,Q3)')

  n_total = df[,.N]
  n_total_rounded5 = round(n_total/5.0)*5

  levs  = levels(df[,get(var)])
  for (j in levs){
    n = df[j,.N,on=var]
    n_rounded5 = round(n/5.0)*5
    tmp_rounded = data.table(Group =j,N=n_rounded5)
    
    tmp <- data.table()
    for (i in c('age_T1DM_DDSC','age_T2DM_DDSC')){
      quants = df[get(var) == j &!is.na(get(i)),quantile(get(i))]
      q1 = round(quants[2],1)
      q2 = round(quants[3],1)
      q3 = round(quants[4],1)
      n_na = df[is.na(get(i)),.N]
      n_na_rounded5 = round(n_na/5.0)*5
      tmp2 = data.table(
        Characteristic=i,
        Group=j, 
        N = n_rounded5,
        Value = c(paste0(q2,' (',q1,',',q3,')'))
        )
      tmp = rbind(tmp,tmp2)
    }
    #print(tmp)
    
    tmp_wide <- reshape(
      tmp[, .(Group, Characteristic, Value)],
      idvar = "Group",
      timevar = "Characteristic",
      v.names = "Value",
      direction = "wide"
    )
    #print(tmp_wide)
    tmp_wide2 = inner_join(tmp_rounded, tmp_wide, by='Group')
    #print(tmp_wide2)
    
    tmpc = rbind(tmpc,tmp_wide2)
  }

  return(tmpc)
}
tmpc_ethnicity_rounded <- table_1_by_var_rounded(df, 'ethnicity')
tmpc_sex_rounded <- table_1_by_var_rounded(df, 'sex')

tmp_by_to_export = rbind(tmpc_ethnicity_rounded,tmpc_sex_rounded)

fwrite(tmpc_ethnicity_rounded,'~output/table1_DMage_by_var.csv')