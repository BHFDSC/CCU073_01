library(DBI)

con <- dbConnect(
odbc::odbc(),
dsn = 'databricks',
HTTPPath = '',
PWD = rstudioapi::askForPassword('Please enter Databricks PAT')
)

setwd("data")

# Transfer the datasets for analysis
for (i in c(2020,2022)) {
  df <- DBI::dbGetQuery(con, paste0("SELECT * FROM .ccu073_01_combine_year",i))
  data.table::fwrite(df,paste0("ccu073_01_year",i,"_",gsub("-","",Sys.Date()),".csv.gz"))
}

 


