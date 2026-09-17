######################################################################################
## TITLE: Generate age db specific hr and calculate corresponding incidence rates 
##
## Author: Wen Shi
## 12-05-2026
#######################################################################################
library(data.table)
library(MASS)
library(survival)

df = readRDS('~data/cohort_2022_20250509_duration.rds')

s1_nonfatal = grep('^nonfatal_.+_date$',colnames(df),value=TRUE)
s1_nonfatal = s1_nonfatal[!s1_nonfatal %in% c("nonfatal_Hypertensive_disease_date",
                                              "nonfatal_OtherHD_date",
                                              "nonfatal_Arter_date",
                                              "nonfatal_OtherCVD_date",
                                              "nonfatal_Sudden_death_date")]
s3_nonfatal = s1_nonfatal[!s1_nonfatal %in% c("nonfatal_Angina_date",
                                              "nonfatal_ARR_date",
                                              "nonfatal_PVD_date",
                                              "nonfatal_DVT_date")]
s1_fatal = grep('^fatal_.+_date$',colnames(df),value=TRUE)

s3_fatal = s1_fatal[!s1_fatal %in% c("fatal_Angina_date",
                                     "fatal_ARR_date",
                                     "fatal_PVD_date",
                                     "fatal_DVT_date")]
df[,c("nonfatal_s3_date","fatal_s3_date"):=.(as.IDate(do.call(pmin,c(mget(s3_nonfatal),na.rm=TRUE))),
                                             as.IDate(do.call(pmin,c(mget(s3_fatal),na.rm=TRUE))))]


df_db2 = df[diabetes_status_DDSC %in% c('No diabetes','T2DM')]
df_db2[diabetes_status_DDSC == 'No diabetes',age_T2DM_DDSC_2:=0]
df_db2[,baseline_age_sq := baseline_age^2]
df_hist0 = df_db2[hist_s3_flag==0]  
df_hist0[,out_flag:=0]
rm(list=c('df','df_db2'))
gc()

transitions = c('ep1_cvd3_nf', 'ep1_cvd3_f', 'ep1_ncvd3_f', 'epf_cvd3_f', 'epf_ncvd3_f')


get_hr = function(fit, b=b, V=V, dbage, baseage,
                  samplesize = 1, sample = FALSE) {
  # Specify comparison
  new1 = data.frame(expo_T2DM_DDSC_hist= factor('Yes', levels=c('No','Yes')), age_T2DM_DDSC_2 = dbage, baseline_age = baseage, 
                    baseline_age_sq = baseage^2)
  new0 = data.frame(expo_T2DM_DDSC_hist= factor('No', levels=c('No','Yes')), age_T2DM_DDSC_2 = 0, baseline_age = baseage, 
                    baseline_age_sq = baseage^2)
  
  # Contrast
  X1 = model.matrix(fit, new1)
  X0 = model.matrix(fit, new0)
  L = as.numeric(X1 - X0)
  
  if (sample) {
  # Simulate hr
  beta_sim = mvrnorm(samplesize, mu = b, Sigma = V)
  loghr_sim = as.numeric(beta_sim %*% L)
  hr_est = exp(loghr_sim)
  return(data.table(dbage = dbage, age = baseage, hr = hr_est, id = seq_len(samplesize)))
  } else {
    loghr = as.numeric(b %*% L)
    hr_est =  exp(loghr)
    se = as.numeric(sqrt(t(L) %*% V %*% L))
    lower = exp(loghr - 1.96*se)
    upper = exp(loghr + 1.96*se)    
    return(data.table(dbage = dbage, age = baseage, hr = hr_est, lhr = lower, uhr = upper))
  }
}

model_transition = function(data, trans, samplesize = 1, sample = FALSE) {
  data = copy(data)
  if (trans == "ep1_cvd3_nf") {
    # transition to non fatal s3
    data[pmin(DOD,data_end,nonfatal_s3_date,na.rm=TRUE)==nonfatal_s3_date & 
               (nonfatal_s3_date != fatal_s3_date|is.na(fatal_s3_date)),out_flag:=1]
    data[out_flag==1,time:=as.integer(difftime(nonfatal_s3_date,baseline_date,units = 'days'))]
    data[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,na.rm = TRUE),baseline_date,units = 'days'))]
    } else if (trans == "ep1_cvd3_f") {
    # transition to fatal s3
      data[pmin(nonfatal_s3_date,data_end,fatal_s3_date,na.rm=TRUE)==fatal_s3_date,
               out_flag:=1]
      data[out_flag==1,time:=as.integer(difftime(fatal_s3_date,baseline_date,units = 'days'))]
      data[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,nonfatal_s3_date,na.rm = TRUE),baseline_date,units = 'days'))]
    } else if (trans == "ep1_ncvd3_f") {
    # transition to non s3 fatal
      data[pmin(DOD,data_end,nonfatal_s3_date,na.rm=TRUE)==DOD & is.na(fatal_s3_date),
               out_flag:=1]
      data[out_flag==1,time:=as.integer(difftime(DOD,baseline_date,units = 'days'))]
      data[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,nonfatal_s3_date,na.rm = TRUE),baseline_date,units = 'days'))]
    } else if (trans == "epf_cvd3_f") {
    # transition non fatal s3 to fatal s3
      data = data[!is.na(nonfatal_s3_date) & nonfatal_s3_date<=data_end]
      data[,baseline_date:=nonfatal_s3_date]
      data[pmin(data_end,fatal_s3_date,na.rm=TRUE)==fatal_s3_date,
              out_flag:=1]
      data[out_flag==1,time:=as.integer(difftime(fatal_s3_date,baseline_date,units = 'days'))]
      data[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,na.rm = TRUE),baseline_date,units = 'days'))]
    } else {
    # transition non fatal s3 to non s3 fatal
      data = data[!is.na(nonfatal_s3_date) & nonfatal_s3_date<=data_end]
      data[,baseline_date:=nonfatal_s3_date]
      data[!is.na(DOD)&is.na(fatal_s3_date)&DOD <= data_end,
              out_flag:=1]
      data[out_flag==1,time:=as.integer(difftime(DOD,baseline_date,units = 'days'))]
      data[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,na.rm = TRUE),baseline_date,units = 'days'))]
  }
  data[,time:=time+1]
  
  fit = coxph(Surv(time,out_flag)~expo_T2DM_DDSC_hist+age_T2DM_DDSC_2+baseline_age+baseline_age_sq
                       +expo_T2DM_DDSC_hist:baseline_age+expo_T2DM_DDSC_hist:baseline_age_sq+age_T2DM_DDSC_2:baseline_age
                       +age_T2DM_DDSC_2:baseline_age_sq
                       ,data=data)
  # Coefs and variance
  b = coef(fit)
  V = vcov(fit)
  
 
  modelled_hr = data.table()
  for (i in seq(30,70,10)) {
    if(i==30) {
      for (j in seq(40,95,5)) {
        modelled_hr = rbind(modelled_hr, get_hr(fit, b, V, i, j, samplesize, sample))}
      } else {
        for (j in seq(i,95,5)) {
          modelled_hr = rbind(modelled_hr, get_hr(fit, b, V, i, j, samplesize, sample))}
        }
  } 
  modelled_hr[, transition:=trans]
 
  modelled_hr  
}  

# get hr by sex
# set.seed(123)
hrl = rbindlist(lapply(c('Male', 'Female'), 
         function(x) {tmp = rbindlist(lapply(transitions, model_transition, data = df_hist0[sex==x]))
          tmp[, sex := x]
          tmp
          }
         )
)

fwrite(hrl,'~output/hr_db_5y_sex.csv')

