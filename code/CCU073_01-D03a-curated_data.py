# Databricks notebook source
# MAGIC %md # CCU073_01-D03a-curated_data
# MAGIC
# MAGIC **Description** This notebook creates the curated tables for hospitalisation and death data.
# MAGIC
# MAGIC **Authors** Wen Shi, Genevieve Cezard
# MAGIC
# MAGIC **Created on** 2023.10.30
# MAGIC
# MAGIC **Last updated on** 2023.12.19
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC hes apc, death
# MAGIC
# MAGIC **Data output** \
# MAGIC CCU073_01_cur_hes_apc_all_years_long\
# MAGIC CCU073_01_cur_deaths_sing\
# MAGIC CCU073_01_cur_deaths_long
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard, Tom Bolton, Wen shi
# MAGIC
# MAGIC **Reviewed** 2024.10.25
# MAGIC
# MAGIC **Acknowledgements** Based on previous work by Tom Bolton (John Nolan, Elena Raffetti) for CCU018_01 and the earlier CCU002 sub-projects.
# MAGIC
# MAGIC **Notes**

# COMMAND ----------

spark.sql('CLEAR CACHE')

# COMMAND ----------

# DBTITLE 1,Libraries
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Window

from functools import reduce

import databricks.koalas as ks
import pandas as pd
import numpy as np

import re
import io
import datetime

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import dates as mdates
import seaborn as sns

print("Matplotlib version: ", matplotlib.__version__)
print("Seaborn version: ", sns.__version__)
_datetimenow = datetime.datetime.now() # .strftime("%Y%m%d")
print(f"_datetimenow:  {_datetimenow}")

# COMMAND ----------

# DBTITLE 1,Functions
# MAGIC %run "/Shared/SHDS/common/functions"

# COMMAND ----------

# MAGIC %md # 0 Parameters

# COMMAND ----------

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

# MAGIC %md # 1 Data

# COMMAND ----------

hes_apc = extract_batch_from_archive(parameters_df_datasets, 'hes_apc')
deaths  = extract_batch_from_archive(parameters_df_datasets, 'deaths')

# COMMAND ----------

# MAGIC %md # 2 HES_APC

# COMMAND ----------

_hes_apc = hes_apc\
   .withColumnRenamed('PERSON_ID_DEID', 'PERSON_ID')\
   .select(['PERSON_ID', 'EPIKEY', 'EPISTART'] + [col for col in list(hes_apc.columns) if re.match(r'^DIAG_(3|4)_\d\d$', col)])\
   .orderBy('PERSON_ID', 'EPIKEY')  
# display(_hes_apc)

# COMMAND ----------

# check
# count_var(hes_apc, 'PERSON_ID_DEID'); print()
# count_var(hes_apc, 'EPIKEY'); print()

# COMMAND ----------

# MAGIC %md ## 2.1 Diag

# COMMAND ----------

# MAGIC %md ### 2.1.1 Create

# COMMAND ----------

# check null EPISTART and potential of using ADMIDATE to supplement
# tmp1 = hes_apc\
#   .select('EPISTART', 'ADMIDATE')\
#   .withColumn('_EPISTART', f.when(f.col('EPISTART').isNotNull(), 1))\
#   .withColumn('_ADMIDATE', f.when(f.col('ADMIDATE').isNotNull(), 1))
# tmpt = tab(tmp1, '_EPISTART', '_ADMIDATE', var2_unstyled=1); print()
# display(hes_apc.where(f.col('EPISTART').isNull()))

# little to be gained from using ADMIDATE and no other date variables 

 

# COMMAND ----------

# reshape twice, tidy, and remove records with missing code
hes_apc_long = reshape_wide_to_long_multi(_hes_apc, i=['PERSON_ID', 'EPIKEY', 'EPISTART'], j='POSITION', stubnames=['DIAG_4_', 'DIAG_3_'])
hes_apc_long = reshape_wide_to_long_multi(hes_apc_long, i=['PERSON_ID', 'EPIKEY', 'EPISTART', 'POSITION'], j='DIAG_DIGITS', stubnames=['DIAG_'])\
  .withColumnRenamed('POSITION', 'DIAG_POSITION')\
  .withColumn('DIAG_POSITION', f.regexp_replace('DIAG_POSITION', r'^[0]', ''))\
  .withColumn('DIAG_DIGITS', f.regexp_replace('DIAG_DIGITS', r'[_]', ''))\
  .withColumn('DIAG_', f.regexp_replace('DIAG_', r'X$', ''))\
  .withColumn('DIAG_', f.regexp_replace('DIAG_', r'[.,\-\s]', ''))\
  .withColumnRenamed('DIAG_', 'CODE')\
  .where((f.col('CODE').isNotNull()) & (f.col('CODE') != ''))\
  .orderBy(['PERSON_ID', 'EPIKEY', 'DIAG_DIGITS', 'DIAG_POSITION'])

# COMMAND ----------

# MAGIC %md ### 2.1.2 Check

# COMMAND ----------

# check
# count_var(hes_apc_long, 'PERSON_ID'); print()
# count_var(hes_apc_long, 'EPIKEY'); print()

# check removal of trailing X
# tmpt = hes_apc_long\
#   .where(f.col('CODE').rlike('X'))\
#   .withColumn('flag', f.when(f.col('CODE').rlike('^X.*'), 1).otherwise(0))
# tmpt = tab(tmpt, 'flag'); print()

# COMMAND ----------

# display(hes_apc_long)

# COMMAND ----------

# MAGIC %md ### 2.1.3 Save

# COMMAND ----------

outName = f'{proj}_cur_hes_apc_all_years_long'.lower()  
hes_apc_long.write.mode('overwrite').saveAsTable(f'{dsa}.{outName}')

# COMMAND ----------

# MAGIC %md ## 2.2 Oper

# COMMAND ----------

# MAGIC %md ### 2.2.1 Create

# COMMAND ----------

# # TODO use OPDATE (c.f., EPISTART), which can be found in the OTR files

# oper_out = hes_apc\
#   .withColumnRenamed('PERSON_ID_DEID', 'PERSON_ID')\
#   .select(['PERSON_ID', 'EPIKEY', 'EPISTART'] + [col for col in list(hes_apc.columns) if re.match(r'^OPERTN_(3|4)_\d\d$', col)])\
#   .orderBy('PERSON_ID', 'EPIKEY')
  
# # reshape twice, tidy, and remove records with missing code
# oper_out_long = reshape_wide_to_long_multi(oper_out, i=['PERSON_ID', 'EPIKEY', 'EPISTART'], j='DIAG_POSITION', stubnames=['OPERTN_4_', 'OPERTN_3_'])
# oper_out_long = reshape_wide_to_long_multi(oper_out_long, i=['PERSON_ID', 'EPIKEY', 'EPISTART', 'DIAG_POSITION'], j='DIAG_DIGITS', stubnames=['OPERTN_'])\
#   .withColumn('DIAG_POSITION', f.regexp_replace('DIAG_POSITION', r'^[0]', ''))\
#   .withColumn('DIAG_DIGITS', f.regexp_replace('DIAG_DIGITS', r'[_]', ''))\
#   .withColumn('OPERTN_', f.regexp_replace('OPERTN_', r'[.,\-\s]', ''))\
#   .withColumnRenamed('OPERTN_', 'CODE')\
#   .where((f.col('CODE').isNotNull()) & (f.col('CODE') != ''))

# COMMAND ----------

# MAGIC %md ### 2.2.2 Check

# COMMAND ----------

# # check
# count_var(oper_out_long, 'PERSON_ID'); print()
# count_var(oper_out_long, 'EPIKEY'); print()

# tmpt = tab(oper_out_long, 'DIAG_DIGITS'); print()
# tmpt = tab(oper_out_long, 'DIAG_POSITION'); print() 
# tmpt = tab(oper_out_long, 'CODE'); print() 
# # TODO - add valid OPCS-4 code checker...

# COMMAND ----------

# display(oper_out_long)

# COMMAND ----------

# MAGIC %md ### 2.2.3 Save

# COMMAND ----------

# outName = f'{proj}_cur_hes_apc_all_years_archive_oper_long'.lower()  
# oper_out_long.write.mode('overwrite').saveAsTable(f'{dbc}.{outName}')
# spark.sql(f'ALTER TABLE {dbc}.{outName} OWNER TO {dbc}')

# COMMAND ----------

# MAGIC %md # 3 Deaths

# COMMAND ----------

# MAGIC %md ## 3.1 Create

# COMMAND ----------

# setting recommended by the data wranglers to avoid the following error message:
# SparkUpgradeException: [INCONSISTENT_BEHAVIOR_CROSS_VERSION.PARSE_DATETIME_BY_NEW_PARSER] 
# You may get a different result due to the upgrading to Spark >= 3.0:
# Caused by: DateTimeParseException: Text '2009029' could not be parsed at index 6
spark.sql("set spark.sql.legacy.timeParserPolicy=LEGACY")

# check
count_var(deaths, 'DEC_CONF_NHS_NUMBER_CLEAN_DEID')
assert dict(deaths.dtypes)['REG_DATE'] == 'string'
assert dict(deaths.dtypes)['REG_DATE_OF_DEATH'] == 'string'

# define window for the purpose of creating a row number below as per the skinny patient table
_win = Window\
  .partitionBy('PERSON_ID')\
  .orderBy(f.desc('REG_DATE'), f.desc('REG_DATE_OF_DEATH'), f.desc('S_UNDERLYING_COD_ICD10'))

# rename ID
# remove records with missing IDs
# reformat dates
# reduce to a single row per individual as per the skinny patient table
# select columns required
# rename column ahead of reshape below
# sort by ID
deaths_out = deaths\
  .withColumnRenamed('DEC_CONF_NHS_NUMBER_CLEAN_DEID', 'PERSON_ID')\
  .where(f.col('PERSON_ID').isNotNull())\
  .withColumn('REG_DATE', f.to_date(f.col('REG_DATE'), 'yyyyMMdd'))\
  .withColumn('REG_DATE_OF_DEATH', f.to_date(f.col('REG_DATE_OF_DEATH'), 'yyyyMMdd'))\
  .withColumn('_rownum', f.row_number().over(_win))\
  .where(f.col('_rownum') == 1)\
  .select(['PERSON_ID', 'REG_DATE', 'REG_DATE_OF_DEATH', 'S_UNDERLYING_COD_ICD10'] + [col for col in list(deaths.columns) if re.match(r'^S_COD_CODE_\d(\d)*$', col)])\
  .withColumnRenamed('S_UNDERLYING_COD_ICD10', 'S_COD_CODE_UNDERLYING')\
  .orderBy('PERSON_ID')


# single row deaths 
deaths_out_sing = deaths_out

# remove records with missing DOD
deaths_out = deaths_out\
  .where(f.col('REG_DATE_OF_DEATH').isNotNull())\
  .drop('REG_DATE')


# reshape
# add 1 to diagnosis position to start at 1 (c.f., 0) - will avoid confusion with HES long, which start at 1
# rename 
# remove records with missing cause of death
deaths_out_long = reshape_wide_to_long(deaths_out, i=['PERSON_ID', 'REG_DATE_OF_DEATH'], j='DIAG_POSITION', stubname='S_COD_CODE_')\
  .withColumn('DIAG_POSITION', f.when(f.col('DIAG_POSITION') != 'UNDERLYING', f.concat(f.lit('SECONDARY_'), f.col('DIAG_POSITION'))).otherwise(f.col('DIAG_POSITION')))\
  .withColumnRenamed('S_COD_CODE_', 'CODE4')\
  .where(f.col('CODE4').isNotNull())\
  .withColumnRenamed('REG_DATE_OF_DEATH', 'DATE')\
  .withColumn('CODE3', f.substring(f.col('CODE4'), 1, 3))
deaths_out_long = reshape_wide_to_long(deaths_out_long, i=['PERSON_ID', 'DATE', 'DIAG_POSITION'], j='DIAG_DIGITS', stubname='CODE')\
  .withColumn('CODE', f.regexp_replace('CODE', r'[.,\-\s]', ''))
  
# check
# count_var(deaths_out_long, 'PERSON_ID')  
# tmpt = tab(deaths_out_long, 'DIAG_POSITION', 'DIAG_DIGITS', var2_unstyled=1) 
# tmpt = tab(deaths_out_long, 'CODE')  
# TODO - add valid ICD-10 code checker...

# COMMAND ----------

# MAGIC %md ## 3.2 Check

# COMMAND ----------

# display(deaths_out_sing)

# COMMAND ----------

# display(deaths_out_long)

# COMMAND ----------

# MAGIC %md ## 3.3 Save

# COMMAND ----------

outName = f'{proj}_cur_deaths_sing'.lower()
deaths_out_sing.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')

# COMMAND ----------

outName = f'{proj}_cur_deaths_long'.lower()
deaths_out_long.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
