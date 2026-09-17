######################################################################################
## TITLE: Create (supplementary) table 2 for CCU073 project
##
## Author: Wen Shi
######################################################################################
library(data.table)
df = readRDS('data/cohort_2022_20241213_v2.rds')
# df = readRDS('data/cohort_2022_20241213_pos1_v2.rds')

tmpc = data.table()
tmpc2 = data.table()
outcomes = c('Angina','CHD','ARR','HF','stroke','PVD','DVT','s1','s2','s3')

for (i in outcomes){
  df[,out_flag:=0]
  
  df[pmin(DOD,data_end,get(paste0('out_',i,'_date')),na.rm=TRUE)==get(paste0('out_',i,'_date')),out_flag:=1]
  df[out_flag==1,time:=as.integer(difftime(get(paste0('out_',i,'_date')),baseline_date,units = 'days'))]
  df[out_flag==0,time:=as.integer(difftime(pmin(DOD,data_end,na.rm = TRUE),baseline_date,units = 'days'))]
  df[,time:=time+1]
  cvd_rate_10_t1 = df[,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                by=c('sex','covid','hist_s1_flag','base_age_cat','age_T1DM_DDSC_10y_cat')][order(sex,covid,hist_s1_flag,base_age_cat,age_T1DM_DDSC_10y_cat)]
  colnames(cvd_rate_10_t1)[which(colnames(cvd_rate_10_t1)=='age_T1DM_DDSC_10y_cat')]='diab_age_cat'
  cvd_rate_10_t1[,c('outcome'):=i]
  cvd_rate_10_t1[,c('diab_type'):='T1DM']
  tmpc = rbind(tmpc,cvd_rate_10_t1)
  cvd_rate_10_t2 = df[,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),
                by=c('sex','covid','hist_s1_flag','base_age_cat','age_T2DM_DDSC_10y_cat')][order(sex,covid,hist_s1_flag,base_age_cat,age_T2DM_DDSC_10y_cat)]
  colnames(cvd_rate_10_t2)[which(colnames(cvd_rate_10_t2)=='age_T2DM_DDSC_10y_cat')]='diab_age_cat'
  cvd_rate_10_t2[,c('outcome'):=i]
  cvd_rate_10_t2[,c('diab_type'):='T2DM']
  tmpc = rbind(tmpc,cvd_rate_10_t2)
  
  df2 = df[hist_s1_flag==0]
  cvd_rate2 = df2[,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000),by='covid']
  
  cvd_rate2[,c('outcome'):=i]
  tmpc2 = rbind(tmpc2,cvd_rate2)
  cvd_rate2_total =  df2[,.(n_event=sum(out_flag),py=sum(time)/365.25,rate=sum(out_flag)/sum(time)*365.25*100000)]
  cvd_rate2_total[,c('covid','outcome'):=.('wholepop',i)]
  tmpc2 = rbind(tmpc2, cvd_rate2_total)
}

tmpc3 = dcast(tmpc2,outcome~covid,value.var = 'n_event')[c('Angina','ARR','CHD','HF','stroke','PVD','DVT','s1','s2','s3'),c('outcome','wholepop','Yes','No')]
tmpc3[,c('wholepop','Yes','No'):=.(round(wholepop/5)*5,round(Yes/5)*5,round(No/5)*5)]

fwrite(tmpc3,'~output/table2_allpos_exclS1_v2.csv')
# fwrite(tmpc3,'~output/table2_pos1_exclS1_v2.csv')

tmpc = tmpc[,c('n_event','py'):=.(round(n_event/5)*5,round(py/5)*5)]
tmpc = tmpc[n_event<10,n_event:=NA]
tmpc = tmpc[py<10,py:=NA]
fwrite(tmpc,'~output/table2_allpos_detail_v2_s123update.csv')

