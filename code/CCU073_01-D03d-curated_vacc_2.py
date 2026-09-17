# Databricks notebook source
# MAGIC %md # CCU073_01-D03d-curated_vacc_2
# MAGIC
# MAGIC **Description** Curated the vaccination dataset and reshape it from long to wide format for CCU073 
# MAGIC
# MAGIC **Author(s)** Tom Bolton, Genevieve Cezard
# MAGIC
# MAGIC **Created on** 2024.04.19 (based on CCU051-D03c-curated_vacc_reshape)
# MAGIC
# MAGIC **Last updated on** 2024.04.19
# MAGIC
# MAGIC **Date last run** 2024.04.19
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC `ccu073_01_cur_vacc`\
# MAGIC `ccu073_01_cur_vacc_qa`
# MAGIC
# MAGIC **Data Output**
# MAGIC - **`ccu073_01_cur_vacc_reshaped`**
# MAGIC
# MAGIC **Reviewers** Tom Bolton, Wen Shi
# MAGIC
# MAGIC **Reviewed** 2024.10.25
# MAGIC
# MAGIC **Acknowledgements** 

# COMMAND ----------

spark.sql('CLEAR CACHE')

# COMMAND ----------

# DBTITLE 1,Libraries
import pyspark.sql.functions as f
from pyspark.sql.types import *
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

spark.sql(f'REFRESH TABLE {dsa}.{proj}_cur_vacc')
vacc = spark.table(f'{dsa}.{proj}_cur_vacc')

spark.sql(f'REFRESH TABLE {dsa}.{proj}_cur_vacc_qa')
vacc_qa = spark.table(f'{dsa}.{proj}_cur_vacc_qa')

# COMMAND ----------

# MAGIC %md # 2 Check

# COMMAND ----------

# display(vacc)

# COMMAND ----------

# display(vacc_qa)

# COMMAND ----------

# count_var(vacc, 'PERSON_ID'); print()
# count_varlist(vacc, ['PERSON_ID', 'DATE'])
# count_varlist(vacc, ['PERSON_ID', 'DATE', 'DATE_rownum'])

# count_var(vacc_qa, 'PERSON_ID'); print()

# COMMAND ----------

# tmpt = tab(vacc, 'dose_sequence'); print()

# COMMAND ----------

# MAGIC %md # 3 Prepare

# COMMAND ----------

# MAGIC %md ## 3.1 Set unknown procedure to Booster

# COMMAND ----------

# Unknown PROCEDURE may relate to dose 3
# <20 records with Unknown PROCEDURE
# set to booster

# check 
# tmpt = tab(vacc, 'PROCEDURE'); print()
count_unknown = vacc.where(f.col('PROCEDURE') == 'Unknown').count()
assert (count_unknown > 0) & (count_unknown < 50) 

# exclude Unknown PROCEDURE
_vacc = (
  vacc
  .withColumn('PROCEDURE', f.when(f.col('PROCEDURE') == 'Unknown', 'B').otherwise(f.col('PROCEDURE')))
)

# check
# tmpt = tab(_vacc, 'PROCEDURE'); print()

# COMMAND ----------

# MAGIC %md # 4 Exclusions

# COMMAND ----------

tmp0 = _vacc
tmpc = count_var(tmp0, 'PERSON_ID', ret=1, df_desc='original', indx=0); print()

# COMMAND ----------

# MAGIC %md ## 4.1 Exclude records before 20201208 (QA4)

# COMMAND ----------

# filter out records before 20201208
tmp0 = (
  tmp0
  .withColumn('qa_4', f.when(f.col('DATE') < f.to_date(f.lit('2020-12-08')), 1).otherwise(0))   
)

# check 
# tmpt = tab(tmp0, 'qa_4'); print()
# tmpt = tabstat(tmp0, 'DATE', byvar='qa_4', date=1); print()

# filter out and tidy
tmp1 = (
  tmp0
  .where(f.col('qa_4') == 0)
  .drop('qa_4')
  .drop('dose_sequence', 'PERSON_ID_rownum', 'PERSON_ID_rownummax')
)

# check
tmpt = count_var(tmp1, 'PERSON_ID', ret=1, df_desc='post exclusion of records before 20201208', indx=1); print()
tmpc = tmpc.unionByName(tmpt)

# check
# print(tmp1.orderBy('PERSON_ID', 'DATE', 'DATE_rownum').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 4.2 Exclude duplicate PERSON_ID and DATE (QA1)

# COMMAND ----------

# QA1
# Keep the record with the latest PROCEDURE (assuming different prodecures) (i.e. if dose 1 and dose 2, keep record for dose 2)
# If same PROCEDURE but different PRODUCT, keep the record at the latest time or keep the first record/randomly if no time specified.

# filter out duplicate PERSON_ID and DATE

_win_PERSON_ID_DATE_ord = Window\
  .partitionBy('PERSON_ID', 'DATE')\
  .orderBy(f.desc('PROCEDURE'), f.desc('TIME'), 'DATE_rownum')
_win_PERSON_ID_DATE = Window\
  .partitionBy('PERSON_ID', 'DATE')
_win_PERSON_ID_ord = Window\
  .partitionBy('PERSON_ID')\
  .orderBy('DATE', f.desc('PROCEDURE'), f.desc('TIME'), 'DATE_rownum')
_win_PERSON_ID_egen = Window\
  .partitionBy('PERSON_ID')\
  .orderBy('DATE', f.desc('PROCEDURE'), f.desc('TIME'), 'DATE_rownum')\
  .rowsBetween(Window.unboundedPreceding, Window.unboundedFollowing)
tmp1 = (
  tmp1
  .withColumn('_PERSON_ID_DATE_rownum', f.row_number().over(_win_PERSON_ID_DATE_ord))
  .withColumn('_PERSON_ID_DATE_rownummax', f.count('PERSON_ID').over(_win_PERSON_ID_DATE))
  .withColumn('_PERSON_ID_rownum', f.row_number().over(_win_PERSON_ID_ord))
  .withColumn('_PERSON_ID_egen_rownummax', f.max(f.col('_PERSON_ID_DATE_rownummax')).over(_win_PERSON_ID_egen))
)

# check
# tmpt = tab(tmp1.where(f.col('_PERSON_ID_DATE_rownum') == 1), '_PERSON_ID_DATE_rownummax'); print()
# tmpt = tab(tmp1.where(f.col('_PERSON_ID_rownum') == 1), '_PERSON_ID_egen_rownummax'); print()
# count_varlist(tmp1, ['PERSON_ID', 'DATE'])

# filter out and tidy
tmp2 = (
  tmp1
  .where(f.col('_PERSON_ID_DATE_rownum') == 1)
  .drop('_PERSON_ID_DATE_rownum', '_PERSON_ID_DATE_rownummax', '_PERSON_ID_rownum', '_PERSON_ID_egen_rownummax')
  .drop('DATE_rownum', 'TIME')
)

# check
# count_varlist(tmp2, ['PERSON_ID', 'DATE'])

# check
tmpt = count_var(tmp2, 'PERSON_ID', ret=1, df_desc='post exclusion of duplicate PERSON_ID and DATE', indx=2); print()
tmpc = tmpc.unionByName(tmpt)

# check
# print(tmp2.orderBy('PERSON_ID', 'DATE').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 4.3 Exclude multiple dose 1 and 2 PROCEDURE (QA2)

# COMMAND ----------

# QA2
# keep earliest DATE of dose 1 and earliest DATE of dose 2 (as per CCU002)?

# note
# - the earliest could be before 20201208
# - the latest record might be more appropriate with subsequent records (e.g., consistent PRODUCT)... 



# check
# tmpt = tab(tmp2, 'PROCEDURE'); print()

# add row number for PERSON_ID and PROCEDURE
_win_PERSON_ID_PROCEDURE_ord = Window\
  .partitionBy('PERSON_ID', 'PROCEDURE')\
  .orderBy('DATE')
tmp2 = (
  tmp2
  .withColumn('PROCEDURE_rownum', f.row_number().over(_win_PERSON_ID_PROCEDURE_ord))
  .withColumn('PROCEDURE_rownum', f.when(f.col('PROCEDURE').isin('1', '2'), f.col('PROCEDURE_rownum')).otherwise(f.lit(1)))
)

# check
# count_varlist(tmp2, ['PERSON_ID', 'PROCEDURE'])
# tmpt = tab(tmp2, 'PROCEDURE_rownum', 'PROCEDURE'); print()

# filter
tmp3 = (
  tmp2
  .where(f.col('PROCEDURE_rownum') == 1)
)

# check
# count_varlist(tmp3, ['PERSON_ID', 'PROCEDURE'])
# tmpt = tab(tmp3, 'PROCEDURE', 'PROCEDURE_rownum'); print()

# tidy
tmp3 = tmp3.drop('PROCEDURE_rownum')

# check
tmpt = count_var(tmp3, 'PERSON_ID', ret=1, df_desc='post exclusion of multiple dose 1 and 2 PROCEDURE ', indx=3); print()
tmpc = tmpc.unionByName(tmpt)

# check
# print(tmp3.orderBy('PERSON_ID', 'DATE').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 4.4 Exclude dose less than 20 days after previous dose (QA6)

# COMMAND ----------

# Note: Could consider an alternative (more conservative) iterative approach
#       For example, say an individual has 3 vaccinations with inter-vaccination intervals of 14 and 14 (e.g., 20210301, 20210315, 20210329)
#       Then, currently the below would exclude vaccinations 2 and 3. An iterative approach would only drop vaccination 2, with the interval between vaccination 1 and 3 being 28 days. 
#       However, the below shows that there are only 35 individuals with more than one inter-vaccination interval < 20 days, therefore it is not a priority to code an iterative approach.

# filter out dose less than 20 days after previous dose
_win_PERSON_ID_ord = Window\
  .partitionBy('PERSON_ID')\
  .orderBy('DATE')
_win_PERSON_ID_egen = Window\
  .partitionBy('PERSON_ID')\
  .orderBy('DATE')\
  .rowsBetween(Window.unboundedPreceding, Window.unboundedFollowing)
tmp3 = (
  tmp3
  .withColumn('_tmp_DATE_diff', f.datediff(f.col('DATE'), f.lag(f.col('DATE'), 1).over(_win_PERSON_ID_ord)))  
  .withColumn('qa_6', f.when(f.col('_tmp_DATE_diff') < 20, 1).otherwise(0))   
  .withColumn('_PERSON_ID_rownum', f.row_number().over(_win_PERSON_ID_ord))
  .withColumn('_PERSON_ID_egen_sum_qa_6', f.sum(f.col('qa_6')).over(_win_PERSON_ID_egen))  
)

# check 
# tmpt = tab(tmp3, 'qa_6'); print()
# tmpt = tabstat(tmp3, '_tmp_DATE_diff', byvar='qa_6'); print()
# tmpt = tab(tmp3.where(f.col('_PERSON_ID_rownum') == 1), '_PERSON_ID_egen_sum_qa_6'); print()

# filter out and tidy
tmp4 = (
  tmp3
  .where(f.col('qa_6') == 0)
  .drop('_tmp_DATE_diff', 'qa_6', '_PERSON_ID_rownum', '_PERSON_ID_egen_sum_qa_6')
)

# check
tmpt = count_var(tmp4, 'PERSON_ID', ret=1, df_desc='post exclusion of dose less than 20 days after previous dose', indx=4); print()
tmpc = tmpc.unionByName(tmpt)

# check
# print(tmp4.orderBy('PERSON_ID', 'DATE').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 4.5 Exclude dose sequence with very small numbers

# COMMAND ----------

# MAGIC %md ### 4.5.1 Recreate dose sequence (QA7)

# COMMAND ----------

print('---------------------------------------------------------------------------------')
print('recreate dose sequence')
print('---------------------------------------------------------------------------------')
# currently PROCEDURE takes values '1', '2', and 'B'
# not possible to distinguish between boosters (i.e., booster 1, booster 2, etc.) by SNOMED codes, simply 'B' 
# so we number boosters according to DATE and DATE_rownum, assuming a dose sequence start number of 3 after first and second dose
_win_cumulative = Window\
  .partitionBy('PERSON_ID')\
  .orderBy('DATE')\
  .rowsBetween(Window.unboundedPreceding, Window.currentRow)

tmp4 = (
  tmp4
  .withColumn('_booster', f.when(f.col('PROCEDURE') == 'B', f.lit(1)).otherwise(0))
  .withColumn('_booster_cumulative', f.sum(f.col('_booster')).over(_win_cumulative))
  .withColumn('_booster_cumulative_plus2', f.col('_booster_cumulative') + f.lit(2))  
  .withColumn('dose_sequence', 
              f.when(f.col('PROCEDURE').isin(['1','2']), f.col('PROCEDURE').cast(t.IntegerType()))
              .when(f.col('PROCEDURE') == 'B', f.col('_booster_cumulative_plus2'))
              .otherwise(999999))
  .drop('_booster', '_booster_cumulative', '_booster_cumulative_plus2')
)

# check
# count_varlist(tmp4, ['PERSON_ID', 'dose_sequence'])
# tmpt = tab(tmp4, 'dose_sequence'); print()
# tmpt = tab(tmp4, 'dose_sequence', 'PROCEDURE'); print()


print('---------------------------------------------------------------------------------')
print('recreate row number')
print('---------------------------------------------------------------------------------')
# add row number that we will compare _dose_seq against
_win_PERSON_ID_ord = Window\
  .partitionBy('PERSON_ID')\
  .orderBy('DATE')
_win_PERSON_ID = Window\
  .partitionBy('PERSON_ID')

tmp4 = (
  tmp4
  .withColumn('PERSON_ID_rownum', f.row_number().over(_win_PERSON_ID_ord))
  .withColumn('PERSON_ID_rownummax', f.count('PERSON_ID').over(_win_PERSON_ID)))

# check
# count_varlist(tmp4, ['PERSON_ID', 'dose_sequence'])
# tmpt = tab(tmp4, 'PERSON_ID_rownum'); print()
# tmpt = tab(tmp4.where(f.col('PERSON_ID_rownum') == 1), 'PERSON_ID_rownummax'); print()


print('---------------------------------------------------------------------------------')
print('compare dose_sequence and PERSON_ID_rownum')
print('---------------------------------------------------------------------------------')
# tmpt = tab(tmp4, 'PERSON_ID_rownum', 'dose_sequence'); print()

# check
# print(tmp4.orderBy('PERSON_ID', 'DATE').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ### 4.5.2 Exclude dose sequence with very big numbers

# COMMAND ----------

# check
# tmpt = tab(tmp4, 'dose_sequence'); print()
# assert tmp4.where(f.col('dose_sequence') > 10).count() < 100 (assertion error)

# filter out
tmp5 = (
  tmp4
  .where(f.col('dose_sequence') <= 10)
)

# check
# tmpt = tab(tmp5, 'dose_sequence'); print()

# check
tmpt = count_var(tmp5, 'PERSON_ID', ret=1, df_desc='post exlusion of dose sequence with very small numbers', indx=5); print()
tmpc = tmpc.unionByName(tmpt)

# check
# print(tmp5.orderBy('PERSON_ID', 'DATE').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 4.6 Flow diagram

# COMMAND ----------

# check flow table
tmpp = (
  tmpc
  .orderBy('indx')
  .select('indx', 'df_desc', 'n', 'n_id', 'n_id_distinct')
  .withColumnRenamed('df_desc', 'stage')
  .toPandas()
)
tmpp['n'] = tmpp['n'].astype(int)
tmpp['n_id'] = tmpp['n_id'].astype(int)
tmpp['n_id_distinct'] = tmpp['n_id_distinct'].astype(int)
tmpp['n_diff'] = (tmpp['n'] - tmpp['n'].shift(1)).fillna(0).astype(int)
tmpp['n_id_diff'] = (tmpp['n_id'] - tmpp['n_id'].shift(1)).fillna(0).astype(int)
tmpp['n_id_distinct_diff'] = (tmpp['n_id_distinct'] - tmpp['n_id_distinct'].shift(1)).fillna(0).astype(int)
for v in [col for col in tmpp.columns if re.search("^n.*", col)]:
  tmpp.loc[:, v] = tmpp[v].map('{:,.0f}'.format)
for v in [col for col in tmpp.columns if re.search(".*_diff$", col)]:  
  tmpp.loc[tmpp['stage'] == 'original', v] = ''
# tmpp = tmpp.drop('indx', axis=1)
print(tmpp.to_string()); print()

# COMMAND ----------

# MAGIC %md # 5 Reshape

# COMMAND ----------

# MAGIC %md ## 5.1 Prepare

# COMMAND ----------

tmpr1 = tmp5

# COMMAND ----------

# Shorten product
tmpr1 = (
  tmpr1
  .withColumn('type', 
              f.when(f.col('PRODUCT') == 'AstraZeneca', 'AZ')
              .when(f.col('PRODUCT').isin(['Pfizer', 'Pfizer child']), 'PB')
              .when(f.col('PRODUCT') == 'Moderna', 'Mo')
              .otherwise(f.lit('Other'))))

# check
# tmpt = tab(tmpr1, 'PRODUCT', 'type'); print()

# tidy
tmpr2 = (
  tmpr1
  .drop('PRODUCT')
)

# check
# print(tmpr2.orderBy('PERSON_ID', 'DATE').limit(10).toPandas().to_string()); print()

# tidy
tmpr3 = (
  tmpr2
  .select('PERSON_ID', f.col('DATE').alias('date'), 'dose_sequence', 'type')
  .orderBy('PERSON_ID', 'dose_sequence')
)  

# check
# count_varlist(tmpr3, ['PERSON_ID', 'dose_sequence']); print()
# count_var(tmpr3, 'PERSON_ID'); print()
# print(tmpr3.orderBy('PERSON_ID', 'dose_sequence').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 5.2 Reshape

# COMMAND ----------

# reshape
tmpr4 = (
  tmpr3
  .groupBy('PERSON_ID')
  .pivot('dose_sequence')
  .agg(
    f.min('date').alias('date')
    , f.first('type').alias('type')
  )
)

# rename columns
dict_rename = {}
for col in tmpr4.columns:
  col_rematch = re.match(r'(\d+)\_(date|type)', col)
  if(col_rematch):
    dict_rename[col] = col_rematch.group(2) + '_' + col_rematch.group(1)
tmpr5 = rename_columns(tmpr4, dict_rename); print()

# check
# count_var(tmpr5, 'PERSON_ID'); print()
# print(len(tmpr5.columns)); print()
# print(pd.DataFrame({f'_cols': tmpr5.columns}).to_string()); print()

# COMMAND ----------

# MAGIC %md # 6 Check

# COMMAND ----------

# check
# display(tmpr5)

# COMMAND ----------

# MAGIC %md # 7 Save

# COMMAND ----------

# save name
outName = f'{proj}_cur_vacc_reshaped'.lower()

# save
tmpr5.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
