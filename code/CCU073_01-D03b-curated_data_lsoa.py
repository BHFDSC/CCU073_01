# Databricks notebook source
# MAGIC %md # CCU073_01-D03b-curated_data_lsoa
# MAGIC
# MAGIC **Description** This notebook creates the curated lsoa lookup tables.
# MAGIC
# MAGIC **Authors** Tom Bolton (Health Data Science Team, BHF Data Science Centre) + Genevieve Cezard (adapted from CCU004_01)
# MAGIC
# MAGIC **Created on** 2023.12.19
# MAGIC
# MAGIC **Last updated on** 2023.12.19
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC ONS geo listing and deprivation indices
# MAGIC
# MAGIC **Data Output**
# MAGIC - **`ccu073_01_cur_lsoa_region_lookup`** : Region lookup for LSOA codes
# MAGIC - **`ccu073_01_cur_lsoa_imd_lookup`** : IMD lookup for LSOA codes
# MAGIC
# MAGIC **Reviewers** Wen Shi
# MAGIC
# MAGIC **Reviewed** 2024.10.25
# MAGIC
# MAGIC **Acknowledgements** Based on CCU004_01.
# MAGIC
# MAGIC **Notes** SQL converted to Pyspark to provide greater transparency and flexibility
# MAGIC

# COMMAND ----------

# MAGIC %md # 0. Setup

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

# DBTITLE 1,Common Functions
# MAGIC %run "/Shared/SHDS/common/functions"

# COMMAND ----------

# MAGIC %md # 1. Parameters

# COMMAND ----------

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

# MAGIC %md # 2. Data

# COMMAND ----------

geog = spark.table(path_ref_geog)
imd = spark.table(path_ref_imd)

# COMMAND ----------

# MAGIC %md # 3. LSOA Region Lookup

# COMMAND ----------

# MAGIC %md ## 3.1 Components

# COMMAND ----------

# MAGIC %md ### 3.1.1 LSOA

# COMMAND ----------

# ENTITY_CODE == E01 relates to LSOA
lsoa = (
  geog
  .where(f.col('ENTITY_CODE') == 'E01')
  .withColumn('flag_parent', f.when(f.col('PARENT_GEOGRAPHY_CODE').rlike(r'^E02.*'), 1).otherwise(0))
)

# check
count_var(lsoa, 'GEOGRAPHY_CODE'); print()
tmpt = tab(lsoa, 'flag_parent'); print() # all E02
assert lsoa.where(f.col('flag_parent') == 0).count() == 0

# prepare
lsoa = (
  lsoa
  .select(
    f.col('GEOGRAPHY_CODE').alias('lsoa_code')
    , f.col('GEOGRAPHY_NAME').alias('lsoa_name')
    , f.col('DATE_OF_TERMINATION').alias('lsoa_dot')
    , f.col('PARENT_GEOGRAPHY_CODE').alias('msoa_code')
  )
)

# identify duplicates 
win_lsoa_code_ord = Window\
  .partitionBy('lsoa_code')\
  .orderBy('lsoa_name', 'lsoa_dot', 'msoa_code')
win_lsoa_code = Window\
  .partitionBy('lsoa_code')
lsoa = (
  lsoa
  .withColumn('rownum', f.row_number().over(win_lsoa_code_ord))
  .withColumn('rownummax', f.count(f.lit(1)).over(win_lsoa_code))
)

# check
tmpt = tab(lsoa.where(f.col('rownum') == 1), 'rownummax'); print()
print(lsoa.orderBy('lsoa_code').where(f.col('rownummax') > 1).toPandas().to_string()); print()

# tidy
lsoa = (
  lsoa
  .drop('rownum', 'rownummax')
)

# check
count_varlist(lsoa, ['lsoa_code', 'lsoa_name']); print()
count_varlist(lsoa, ['lsoa_code', 'lsoa_name', 'lsoa_dot']); print()
count_varlist(lsoa, ['lsoa_code', 'lsoa_name', 'lsoa_dot', 'msoa_code']); print()
print(lsoa.orderBy('lsoa_code', 'lsoa_name', 'lsoa_dot', 'msoa_code').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ### 3.1.2 MSOA

# COMMAND ----------

# ENTITY_CODE == E02 relates to MSOA
msoa = (
  geog
  .where(f.col('ENTITY_CODE') == 'E02')
  .withColumn('PARENT_GEOGRAPHY_CODE_3', f.substring(f.col('PARENT_GEOGRAPHY_CODE'), 1, 3))
  .withColumn('flag_parent', f.when((f.col('PARENT_GEOGRAPHY_CODE_3').isin(['E06', 'E07', 'E08', 'E09'])) | (f.col('PARENT_GEOGRAPHY_CODE_3').isNull()), 1).otherwise(0))
)  

# check
count_var(msoa, 'GEOGRAPHY_CODE'); print()
tmpt = tab(msoa, 'PARENT_GEOGRAPHY_CODE_3'); print() # all E06, E07, E08, E09
tmpt = tab(msoa, 'PARENT_GEOGRAPHY_CODE_3', 'flag_parent'); print() 
assert msoa.where(f.col('flag_parent') == 0).count() == 0

# prepare
msoa = (
  msoa
  .select(
    f.col('GEOGRAPHY_CODE').alias('msoa_code')
    , f.col('GEOGRAPHY_NAME').alias('msoa_name')
    , f.col('DATE_OF_TERMINATION').alias('msoa_dot')
    , f.col('PARENT_GEOGRAPHY_CODE').alias('la_code')
  )
)

# check
count_var(msoa, 'msoa_code'); print()
print(msoa.orderBy('msoa_code').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ### 3.1.3 LA

# COMMAND ----------

# ENTITY_CODE == E06, E07, E08, E09 relate to LA (Local Authority)
la = (
  geog
  .where(f.col('ENTITY_CODE').isin(['E06', 'E07', 'E08', 'E09']))
  .withColumn('PARENT_GEOGRAPHY_CODE_3', f.substring(f.col('PARENT_GEOGRAPHY_CODE'), 1, 3))
  .withColumn('flag_parent', f.when(f.col('PARENT_GEOGRAPHY_CODE_3').isin(['E10', 'E12']), 1).otherwise(0))
)  

# check
count_var(la, 'GEOGRAPHY_CODE'); print()
tmpt = tab(la, 'PARENT_GEOGRAPHY_CODE_3'); print() # all E10, E12
tmpt = tab(la, 'PARENT_GEOGRAPHY_CODE_3', 'flag_parent'); print() 
assert la.where(f.col('flag_parent') == 0).count() == 0

# prepare
la = (
  la
  .select(
    f.col('GEOGRAPHY_CODE').alias('la_code')
    , f.col('GEOGRAPHY_NAME').alias('la_name')
    , f.col('DATE_OF_TERMINATION').alias('la_dot')    
    , 'PARENT_GEOGRAPHY_CODE'
    , 'PARENT_GEOGRAPHY_CODE_3'
  )
  .withColumn('county_code', f.when(f.col('PARENT_GEOGRAPHY_CODE_3') == 'E10', f.col('PARENT_GEOGRAPHY_CODE')))
  .withColumn('region_code_1', f.when(f.col('PARENT_GEOGRAPHY_CODE_3') == 'E12', f.col('PARENT_GEOGRAPHY_CODE')))
  .drop('PARENT_GEOGRAPHY_CODE', 'PARENT_GEOGRAPHY_CODE_3')
)  

# check
count_var(la, 'la_code'); print()

# identify duplicates 
win_la_code_ord = Window\
  .partitionBy('la_code')\
  .orderBy('la_name', 'la_dot', 'county_code', 'region_code_1')
win_la_code = Window\
  .partitionBy('la_code')
la = (
  la
  .withColumn('rownum', f.row_number().over(win_la_code_ord))
  .withColumn('rownummax', f.count(f.lit(1)).over(win_la_code))
)
              
# check
tmpt = tab(la.where(f.col('rownum') == 1), 'rownummax'); print()
print(la.orderBy('la_code').where(f.col('rownummax') > 1).toPandas().to_string()); print()
assert la.where(f.col('rownum') == 1).where(f.col('rownummax') > 1).count() == 1
assert la.where(f.col('rownummax') > 1).count() == 2
assert la.where(f.col('rownummax') > 1).where(f.col('rownum') == 1).select('la_code').collect()[0][0] == 'E07000112'
assert la.where(f.col('rownummax') > 1).where(f.col('rownum') == 1).select('la_name').collect()[0][0] == 'Folkestone and Hythe'
assert la.where(f.col('rownummax') > 1).where(f.col('rownum') == 2).select('la_name').collect()[0][0] == 'Shepway'
assert la.where(f.col('rownummax') > 1).where(f.col('rownum') == 1).select('county_code').collect()[0][0] == 'E10000016'
assert la.where(f.col('rownummax') > 1).where(f.col('rownum') == 2).select('county_code').collect()[0][0] == 'E10000016'

# note: https://en.wikipedia.org/wiki/Folkestone_and_Hythe_District
# Folkestone and Hythe is a local government district in Kent ... The authority was renamed from Shepway in April 2018 ...

# filter
# tidy
la = (
  la
  .where(f.col('rownum') == 1)
  .drop('rownum', 'rownummax')
)

# check
count_var(la, 'la_code'); print()
print(la.where(f.col('la_code') == 'E07000112').toPandas().to_string()); print()
print(la.orderBy('la_code').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ### 3.1.4 County 

# COMMAND ----------

# ENTITY_CODE == E10 relates to County
county = (
  geog
  .where(f.col('ENTITY_CODE') == 'E10')
  .withColumn('PARENT_GEOGRAPHY_CODE_3', f.substring(f.col('PARENT_GEOGRAPHY_CODE'), 1, 3))
  .withColumn('flag_parent', f.when(f.col('PARENT_GEOGRAPHY_CODE_3').isin(['E12']), 1).otherwise(0))
)  

# check
count_var(county, 'GEOGRAPHY_CODE'); print()
tmpt = tab(county, 'PARENT_GEOGRAPHY_CODE_3'); print() # all E12
tmpt = tab(county, 'PARENT_GEOGRAPHY_CODE_3', 'flag_parent'); print() 
assert county.where(f.col('flag_parent') == 0).count() == 0

# prepare
county = (
  county
  .select(
    f.col('GEOGRAPHY_CODE').alias('county_code')
    , f.col('GEOGRAPHY_NAME').alias('county_name')
    , f.col('DATE_OF_TERMINATION').alias('county_dot')
    , f.col('PARENT_GEOGRAPHY_CODE').alias('region_code_2')
  )
)

# check
count_var(county, 'county_code'); print()
print(county.orderBy('county_code').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ### 3.1.5 Region 

# COMMAND ----------

# ENTITY_CODE == E12 relates to Region
region = (
  geog
  .where(f.col('ENTITY_CODE') == 'E12')
)  

# check
count_var(region, 'GEOGRAPHY_CODE'); print()
assert region.count() == 9

# prepare
region = (
  region
  .select(
    f.col('GEOGRAPHY_CODE').alias('region_code')
    , f.col('GEOGRAPHY_NAME').alias('region_name')
    , f.col('DATE_OF_TERMINATION').alias('region_dot')
  )
)

# check
count_var(region, 'region_code'); print()
print(region.orderBy('region_code').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 3.2 Join

# COMMAND ----------

# check
count_var(lsoa, 'lsoa_code'); print()
count_var(msoa, 'msoa_code'); print()
count_var(la, 'la_code'); print()
count_var(county, 'county_code'); print()
count_var(region, 'region_code'); print()


print('-----------------------------------------------------------------')
print('merge LSOA and MSOA')
print('-----------------------------------------------------------------')
tmpg1 = merge(lsoa, msoa, ['msoa_code'], validate='m:1', assert_results=['both'], keep_results=['both'], indicator=0); print()


print('-----------------------------------------------------------------')
print('merge above and LA')
print('-----------------------------------------------------------------')
# prepare
tmpg1 = (
  tmpg1
  .withColumn('flag_la_code_null', f.when(f.col('la_code').isNull(), 1).otherwise(0))
)

# merge
tmpg2 = merge(tmpg1, la, ['la_code'], validate='m:1', keep_results=['both', 'left_only']); print()

# check
tmpt = tab(tmpg2, '_merge', 'flag_la_code_null'); print()

# tidy
tmpg2 = (
  tmpg2
  .drop('_merge', 'flag_la_code_null')
)


print('-----------------------------------------------------------------')
print('merge above and County')
print('-----------------------------------------------------------------')
# prepare
tmpg2 = (
  tmpg2
  .withColumn('flag_county_code_null', f.when(f.col('county_code').isNull(), 1).otherwise(0))
)

# merge
tmpg3 = merge(tmpg2, county, ['county_code'], validate='m:1', keep_results=['both', 'left_only']); print()

# check
tmpt = tab(tmpg3, '_merge', 'flag_county_code_null'); print()

# tidy
tmpg3 = (
  tmpg3
  .withColumn('region_code', f.coalesce(f.col('region_code_1'), f.col('region_code_2')))
  .drop('_merge', 'flag_county_code_null', 'region_code_1', 'region_code_2')
)


print('-----------------------------------------------------------------')
print('merge above and Region')
print('-----------------------------------------------------------------')
# prepare
tmpg3 = (
  tmpg3
  .withColumn('flag_region_code_null', f.when(f.col('region_code').isNull(), 1).otherwise(0))
)

# merge
tmpg4 = merge(tmpg3, region, ['region_code'], validate='m:1', keep_results=['both', 'left_only']); print()

# check
tmpt = tab(tmpg4, '_merge', 'flag_region_code_null'); print()

# tidy
tmpg5 = (
  tmpg4
  .drop('_merge', 'flag_region_code_null')
  .select(
    'lsoa_code', 'lsoa_name', 'lsoa_dot'
    , 'msoa_code', 'msoa_name', 'msoa_dot'
    , 'la_code', 'la_name', 'la_dot'
    , 'county_code', 'county_name', 'county_dot'
    , 'region_code', 'region_name', 'region_dot'
  )
  .orderBy('lsoa_code', 'lsoa_name', 'lsoa_dot', 'msoa_code')
  .withColumn('flag_all_dot_null', 
              f.when(
                (f.col('lsoa_dot').isNull())
                & (f.col('msoa_dot').isNull()) 
                & (f.col('la_dot').isNull()) 
                & (f.col('county_dot').isNull()) 
                & (f.col('region_dot').isNull())
              , 1)
              .otherwise(0)
             )
)

# check
tmpt = tab(tmpg5, 'flag_all_dot_null'); print()

# COMMAND ----------

# check
count_varlist(tmpg5, ['lsoa_code'])

# COMMAND ----------

# check
# display(tmpg5)

# COMMAND ----------

# MAGIC %md ## 3.3 Exclude LSOA with missing region

# COMMAND ----------

# filter LSOA with missing region (due to missing LA)

# check
count_varlist(tmpg5, ['lsoa_code', 'region_name'])

# check
tmpf = (
  tmpg5
  .withColumn('la_code_null', f.when(f.col('la_code').isNull(), 1).otherwise(0))
  .withColumn('region_name_null', f.when(f.col('region_name').isNull(), 1).otherwise(0))
)
tmpt = tab(tmpf, 'la_code_null', 'region_name_null'); print()

# filter
tmpg6 = (
  tmpg5
  .where(f.col('la_code').isNotNull())
)

# check
count_varlist(tmpg6, ['lsoa_code', 'region_name'])
count_varlist(tmpg6, ['lsoa_code'])

# COMMAND ----------

tmpf = (
  tmpg5
  .where(f.col('la_code').isNull())
)
tmpt = tab(tmpf, 'msoa_dot'); print()
# all msoa_dot 2011-12-30 - old termination date

# COMMAND ----------

# check
# display(tmpf.orderBy('lsoa_code'))

# COMMAND ----------

# check
# display(tmpg5.where(f.col('lsoa_code') == 'E01000019'))

# COMMAND ----------

# MAGIC %md ## 3.4 Exclude LSOA duplicates

# COMMAND ----------

# check duplicate LSOAs for region conflict
w1 = Window\
  .partitionBy('lsoa_code')\
  .orderBy('region_name')
w2 = Window\
  .partitionBy('lsoa_code')
w3 = Window\
  .partitionBy('lsoa_code')\
  .orderBy('region_name')\
  .rowsBetween(Window.unboundedPreceding, Window.unboundedFollowing)
tmpf = (
  tmpg6
  .withColumn('rownum', f.row_number().over(w1))
  .withColumn('rownummax', f.count(f.lit(1)).over(w2))
  .withColumn('region_name_lag1', f.lag(f.col('region_name'), 1).over(w1))
  .withColumn('region_lag1_diff', f.when(f.col('region_name') != f.col('region_name_lag1'), 1).otherwise(0))
  .withColumn('region_name_conflict', f.max(f.col('region_lag1_diff')).over(w3))
  .withColumn('region_name_null', f.when(f.col('region_name').isNull(), 1).otherwise(0))
)
  
# check  
tmpt = tab(tmpf, 'region_name_null'); print()
tmpt = tab(tmpf.where(f.col('rownum') == 1), 'rownummax', 'region_name_conflict'); print()  

# COMMAND ----------

# select unique LSOA (null DATE_OF_TERMINATION i.e., current)

# check
w1 = Window\
  .partitionBy('lsoa_code')\
  .orderBy(f.asc_nulls_first('lsoa_dot')) # default
w2 = Window\
  .partitionBy('lsoa_code')
tmpg7 = (
  tmpg6
  .withColumn('rownum', f.row_number().over(w1))
  .withColumn('rownummax', f.count(f.lit(1)).over(w2))
  .withColumn('flag_lsoa_dot_null', f.when(f.col('lsoa_dot').isNull(), 1).otherwise(0))
  .withColumn('flag_lsoa_dot_not_null', f.when(f.col('lsoa_dot').isNotNull(), 1).otherwise(0))
  .withColumn('sum_lsoa_dot_null', f.sum(f.col('flag_lsoa_dot_null')).over(w2))
  .withColumn('sum_lsoa_dot_not_null', f.sum(f.col('flag_lsoa_dot_not_null')).over(w2))
)
  
# check  
tmpt = tab(tmpg7, 'flag_lsoa_dot_null', 'flag_lsoa_dot_not_null'); print()
tmpt = tab(tmpg7.where(f.col('rownum') == 1), 'rownummax'); print() 
tmpt = tab(tmpg7.where(f.col('rownum') == 1), 'sum_lsoa_dot_null', 'sum_lsoa_dot_not_null'); print() 
tmpt = tab(tmpg7.where((f.col('rownum') == 1) & (f.col('rownummax') > 1)), 'sum_lsoa_dot_null', 'sum_lsoa_dot_not_null'); print() 

# filter
tmpg8 = (
  tmpg7
  .where(f.col('rownum') == 1)
)

# check  
tmpt = tab(tmpg8, 'flag_lsoa_dot_null', 'flag_lsoa_dot_not_null'); print()

# tidy
tmpg8 = (
  tmpg8
  .drop('rownum', 'rownummax', 'flag_lsoa_dot_not_null', 'sum_lsoa_dot_null', 'sum_lsoa_dot_not_null')
)

# check
count_var(tmpg8, 'lsoa_code'); print()
tmpt = tab(tmpg8, 'flag_all_dot_null', 'flag_lsoa_dot_null'); print()

# COMMAND ----------

# MAGIC %md ## 3.5 Check

# COMMAND ----------

# display(tmpg7.where(f.col('rownummax') > 1))

# COMMAND ----------

# check
# display(tmpg8)

# COMMAND ----------

# check
count_var(tmpg8, 'lsoa_code'); print()

# COMMAND ----------

# MAGIC %md ## 3.6 Save

# COMMAND ----------

# Select variables of interest only
tmpg8 = ( 
    tmpg8.select( f.col('lsoa_code').alias('lsoa'),  f.col('lsoa_name'),  f.col('region_code'),  f.col('region_name').alias('region'))
)

# COMMAND ----------

outName = f'{proj}_cur_lsoa_region_lookup'
tmpg8.write.mode('overwrite').saveAsTable(f'{dsa}.{outName}')

# COMMAND ----------

# MAGIC
# MAGIC %md # 4. LSOA IMD Lookup

# COMMAND ----------

# MAGIC %md ## 4.1 Create

# COMMAND ----------

# check
print(imd.toPandas().head(5)); print()
count_var(imd, 'LSOA_CODE_2011'); print()
ltmpt = tab(imd, 'DECI_IMD', 'IMD_YEAR', var2_unstyled=1); print()

# tidy
ltmp1 = imd\
  .where(f.col('IMD_YEAR') == 2019)\
  .select('LSOA_CODE_2011', 'DECI_IMD')\
  .withColumnRenamed('LSOA_CODE_2011', 'LSOA')\
  .withColumn('IMD_2019_QUINTILES',
    f.when(f.col('DECI_IMD').isin([1,2]), 1)\
     .when(f.col('DECI_IMD').isin([3,4]), 2)\
     .when(f.col('DECI_IMD').isin([5,6]), 3)\
     .when(f.col('DECI_IMD').isin([7,8]), 4)\
     .when(f.col('DECI_IMD').isin([9,10]), 5)\
     .otherwise(None)\
  )\
  .withColumnRenamed('DECI_IMD', 'IMD_2019_DECILES')

# check
ltmpt = tab(ltmp1, 'IMD_2019_DECILES', 'IMD_2019_QUINTILES', var2_unstyled=1); print()
print(ltmp1.toPandas().head(5)); print()

# COMMAND ----------

# MAGIC %md ## 4.2 Save

# COMMAND ----------

outName=f'{proj}_cur_lsoa_imd_lookup'
ltmp1.write.mode('overwrite').saveAsTable(f'{dsa}.{outName}')
