######################################################################################
## TITLE: Calculate sex specific IR for each transition type
##
## Author: Wen Shi
## 12-05-2026
#######################################################################################  

library(data.table)

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



ir_agegp = function(df = df, agevar, sex = TRUE, covid = FALSE) {
  df_hist0 = df[hist_s3_flag==0]
  
  # Transition from well to non fatal cvd s3
  df_hist0[,out_flag:=0]
  
  df_hist0[pmin(DOD,data_end,nonfatal_s3_date,na.rm=TRUE)==nonfatal_s3_date & 
             (nonfatal_s3_date != fatal_s3_date|is.na(fatal_s3_date)),out_flag:=1]
  df_hist0[out_flag==1,time:=as.integer(difftime(nonfatal_s3_date,baseline_date,units = 'days'))]
  df_hist0[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,na.rm = TRUE),baseline_date,units = 'days'))]
  df_hist0[,time:=time+1]
  
  if (sex & covid) {
  cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                         by=.(sex, covid, age = get(agevar))][order(sex, covid, age)]
  } else if (sex & !covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, age = get(agevar))][order(sex, age)]
  } else if (!sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(covid, age = get(agevar))][order(covid, age)]
  } else {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(age = get(agevar))][order(age)]
  }
  cvd_rate_10_t2[,c('transition'):='ep1_cvd3_nf']
  tmpc = cvd_rate_10_t2
  print("well to nonfatal cvd complete")
  
  # Transition from well to fatal cvd s3
  df_hist0[,out_flag:=0]
  
  df_hist0[pmin(nonfatal_s3_date,data_end,fatal_s3_date,na.rm=TRUE)==fatal_s3_date,
           out_flag:=1]
  df_hist0[out_flag==1,time:=as.integer(difftime(fatal_s3_date,baseline_date,units = 'days'))]
  df_hist0[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,nonfatal_s3_date,na.rm = TRUE),baseline_date,units = 'days'))]
  df_hist0[,time:=time+1]
  
  if (sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, covid, age = get(agevar))][order(sex, covid, age)]
  } else if (sex & !covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, age = get(agevar))][order(sex, age)]
  } else if (!sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(covid, age = get(agevar))][order(covid, age)]
  } else {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(age = get(agevar))][order(age)]
  }
  cvd_rate_10_t2[,c('transition'):='ep1_cvd3_f']
  tmpc = rbind(tmpc,cvd_rate_10_t2)
  print("well to fatal cvd complete")
  
  # Transition from well to fatal non-cvd
  df_hist0[,out_flag:=0]
  
  df_hist0[pmin(DOD,data_end,nonfatal_s3_date,na.rm=TRUE)==DOD & is.na(fatal_s3_date),
           out_flag:=1]
  df_hist0[out_flag==1,time:=as.integer(difftime(DOD,baseline_date,units = 'days'))]
  df_hist0[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,nonfatal_s3_date,na.rm = TRUE),baseline_date,units = 'days'))]
  df_hist0[,time:=time+1]
  
  if (sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, covid, age = get(agevar))][order(sex, covid, age)]
  } else if (sex & !covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, age = get(agevar))][order(sex, age)]
  } else if (!sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(covid, age = get(agevar))][order(covid, age)]
  } else {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(age = get(agevar))][order(age)]
  }
  cvd_rate_10_t2[,c('transition'):='ep1_ncvd3_f']
  tmpc = rbind(tmpc,cvd_rate_10_t2)
  print("well to fatal noncvd complete")
  
  # Transition from nonfatal cvd to fatal cvd
  df_hist0 = df_hist0[!is.na(nonfatal_s3_date) & nonfatal_s3_date<=data_end]
  df_hist0[,baseline_date:=nonfatal_s3_date]
  df_hist0[,out_flag:=0]
  df_hist0[pmin(data_end,fatal_s3_date,na.rm=TRUE)==fatal_s3_date,
          out_flag:=1]
  df_hist0[out_flag==1,time:=as.integer(difftime(fatal_s3_date,baseline_date,units = 'days'))]
  df_hist0[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,na.rm = TRUE),baseline_date,units = 'days'))]
  df_hist0[,time:=time+1]
  
  if (sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, covid, age = get(agevar))][order(sex, covid, age)]
  } else if (sex & !covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, age = get(agevar))][order(sex, age)]
  } else if (!sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(covid, age = get(agevar))][order(covid, age)]
  } else {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(age = get(agevar))][order(age)]
  }
  cvd_rate_10_t2[,c('transition'):='epf_cvd3_f']
  tmpc = rbind(tmpc,cvd_rate_10_t2)
  print("nonfatal to fatal including hist s3 complete")
  
  # Transition from first nonfatal cvd s3 to fatal non-cvd 
  df_hist0[,out_flag:=0]
  
  df_hist0[!is.na(DOD)&is.na(fatal_s3_date)&DOD <= data_end,
          out_flag:=1]
  df_hist0[out_flag==1,time:=as.integer(difftime(DOD,baseline_date,units = 'days'))]
  df_hist0[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,na.rm = TRUE),baseline_date,units = 'days'))]
  df_hist0[,time:=time+1]
  
  if (sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, covid, age = get(agevar))][order(sex, covid, age)]
  } else if (sex & !covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(sex, age = get(agevar))][order(sex, age)]
  } else if (!sex & covid) {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(covid, age = get(agevar))][order(covid, age)]
  } else {
    cvd_rate_10_t2 = df_hist0[expo_T2DM_DDSC_10y_age_cat_2=='No diabetes'][,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                                                                           by=.(age = get(agevar))][order(age)]
  }
  cvd_rate_10_t2[,c('transition'):='epf_ncvd3_f']
  tmpc = rbind(tmpc,cvd_rate_10_t2)
  print("nonfatal to fatal non-cvd including hist s3 complete")
  
  tmpc = tmpc[,c('n_event','py'):=.(round(n_event/5)*5,round(py/5)*5)]
  tmpc = tmpc[n_event<10,n_event:=NA]
  tmpc = tmpc[py<10,py:=NA]
  tmpc
}


df[,baseage_5y := cut(baseline_age, 
                      breaks = c(seq(40,95,5),122), 
                      labels = c(paste0(seq(40,90,5),'_',seq(44,94,5)),'95_over'), 
                      right = FALSE)]
ir_nodb = ir_agegp(df = df, agevar = 'baseage_5y', sex = TRUE, covid = FALSE)

fwrite(ir_nodb,'~output/ir_nodb_5y_sex.csv')