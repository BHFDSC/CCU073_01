######################################################################################
## TITLE: Get extra descriptive statistics for maps by diabetes characteristics
##
## Author: Wen Shi
######################################################################################
install.packages("R.utils", repos = "https://packages.sde.digital.nhs.uk/repository/cran-mirror/", dependencies = TRUE)
library(R.utils)
library(data.table)
df = readRDS('~data/cohort_2022_20250509_duration.rds')
df = df[hist_s1_flag==0]
la = fread('~data/ccu073_01_year2022_la_20250619.csv.gz')
df[,.N]
tmpc = la[df,on='PERSON_ID']
tmpc[,.N]
tmpc[,length(unique(PERSON_ID))]
tmpc[is.na(la),.N]
tmpc[la%in%c('Waveney','Suffolk Coastal'),la:='East Suffolk']
tmpc[la %in%c('St Edmundsbury','Forest Heath'),la:='West Suffolk']
tmpc[la %in%c('Northampton','Daventry','South Northamptonshire'),la:='West Northamptonshire']
tmpc[la %in%c('Kettering','Wellingborough','Corby','East Northamptonshire'),la:='North Northamptonshire']
tmpc[la %in%c('Bournemouth','Poole','Christchurch'),la:='Bournemouth, Christchurch and Poole']
tmpc[la %in%c('East Dorset','West Dorset','North Dorset','Weymouth and Portland','Purbeck'),la:='Dorset']
tmpc[la %in%c('Wycombe','Chiltern','South Bucks','Aylesbury Vale'),la:='Buckinghamshire']
tmpc[la %in%c('Taunton Deane','West Somerset'),la:='Somerset West and Taunton']

totn = tmpc[,.N,by=.(la)]
dbn  =tmpc[diabetes_status_DDSC=='T2DM',.(n=.N),by=la]
tmp = dbn[totn,on='la']
tmp[,p:=n/N]
dbage = tmpc[diabetes_status_DDSC=='T2DM',.(m=median(age_T2DM_DDSC_2)),by=la]
tmp = dbage[tmp,on='la']
tmp[,c('n','N'):=.(round(n/5)*5,round(N/5)*(5))]
tmp[n<10,n:=NA]
tmp[N<10,N:=NA]
setnames(tmp,c('m','n'),c('db2_age','n_db2'))
fwrite(tmp,'~output/db2_la_rmCVD.csv')
