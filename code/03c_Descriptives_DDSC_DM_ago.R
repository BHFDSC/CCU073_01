######################################################################################
## TITLE: Get number of people for each step of the DDSC algorithm
##
## Author: Genevieve Cezard
######################################################################################
library(dplyr)

#-------------------------------------------------------------#
# 1. Get PERSON_ID for those included in the study population #
#-------------------------------------------------------------#
df = readRDS('~data/cohort_2022_20241213_v2.rds')
df = df[hist_s1_flag==0] # Exclude people with prevalent cvd (s1 definition)

#-----------------------------------------#
# 2. Get info on the DDSC algorithm steps #
#-----------------------------------------#
df_DDSC_DM_type <- data.table::fread("data/ccu073_01_ddsc_diabetes_type_2022_01_20241202.csv.gz", data.table = FALSE)

#----------#
# 3. Merge #
#----------#
df_select <- df %>% select(PERSON_ID, sex, ethnicity)
df <- left_join(df_select, df_DDSC_DM_type, by = "PERSON_ID")

#------------------------#
# 4. Get step counts out #
#------------------------#
# Using SQL in R    
library(odbc)
library(DBI) 

con <- dbConnect(drv = RSQLite::SQLite(), dbname = ":memory:")
dbWriteTable(conn = con, name = "df", value = df)

step_1 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_1 AS Category, 
                                          COUNT(PERSON_ID) AS step_1
                          FROM df
                          GROUP BY df.step_1"
)
step_2 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_2 AS Category, 
                                          COUNT(PERSON_ID) AS step_2
                          FROM df
                          GROUP BY df.step_2"
)
step_3 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_3 AS Category, 
                                          COUNT(PERSON_ID) AS step_3
                          FROM df
                          GROUP BY df.step_3"
)
step_3_1 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_3_1 AS Category, 
                                          COUNT(PERSON_ID) AS step_3_1
                          FROM df
                          GROUP BY df.step_3_1"
)
step_4 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_4 AS Category, 
                                          COUNT(PERSON_ID) AS step_4
                          FROM df
                          GROUP BY df.step_4"
)
step_5 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_5 AS Category, 
                                          COUNT(PERSON_ID) AS step_5
                          FROM df
                          GROUP BY df.step_5"
)
step_6 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_6 AS Category, 
                                          COUNT(PERSON_ID) AS step_6
                          FROM df
                          GROUP BY df.step_6"
)
step_6_1 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_6_1 AS Category, 
                                          COUNT(PERSON_ID) AS step_6_1
                          FROM df
                          GROUP BY df.step_6_1"
)
step_7 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_7 AS Category, 
                                          COUNT(PERSON_ID) AS step_7
                          FROM df
                          GROUP BY df.step_7"
)
step_8 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_8 AS Category, 
                                          COUNT(PERSON_ID) AS step_8
                          FROM df
                          GROUP BY df.step_8"
)
step_9 <- DBI::dbGetQuery(conn = con, statement = "
                          SELECT DISTINCT df.step_9 AS Category, 
                                          COUNT(PERSON_ID) AS step_9
                          FROM df
                          GROUP BY df.step_9"
)

df_list <- list(step_1, step_2, step_3, step_3_1, step_4, step_5, step_6, step_6_1, step_7, step_8, step_9)      

#merge all data frames together
N_steps <- df_list %>% reduce(full_join, by='Category')
N_steps_rounded <- N_steps %>% mutate(across(where(is.numeric), ~ round(./5.0)*5, 0))

fwrite(N_steps_rounded,'~output/DDSC_DM_algo_steps_counts.csv') # for current cohort excluding those with prior CVD
#fwrite(N_steps_rounded,'~output/DDSC_DM_algo_steps_counts_previousCVD.csv') # for previous cohort including those with prior CVD

# Check resulting diabetes classification
table(df$out_diabetes)
