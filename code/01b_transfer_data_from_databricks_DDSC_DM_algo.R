library(DBI)

con <- dbConnect(
odbc::odbc(),
dsn = 'databricks',
HTTPPath = '',
PWD = rstudioapi::askForPassword('Please enter Databricks PAT')
)

setwd("data")

# Transfer the DDSC dataset with diabetes type and date at diagnosis
df <- DBI::dbGetQuery(con, paste0("SELECT * FROM .ccu073_01_ddsc_cohort_out_2024_07_25"))
data.table::fwrite(df,paste0("ccu073_01_ddsc_diabetes_type_2022_01_",gsub("-","",Sys.Date()),".csv.gz"))


# Transfer the DDSC dataset which contains further variables
df <- DBI::dbGetQuery(con, paste0("SELECT * FROM .ccu073_01_ddsc_cohort_2024_07_25"))
data.table::fwrite(df,paste0("ccu073_01_ddsc_diabetes_all_var_2022_01_",gsub("-","",Sys.Date()),".csv.gz"))


