# Databricks notebook source
# MAGIC %md # CCU073_01-D04b-LSOA
# MAGIC  
# MAGIC **Description** This notebook creates the covariates - region,deprivation based on LSOA using multi-sourced LSOA table updated monthly by BHF data science centre.
# MAGIC
# MAGIC
# MAGIC
# MAGIC LSOA conflicts (different `LSOA` values which have the same `REPORTING_PERIOD_END_DATE`) will be solved using mode.
# MAGIC  
# MAGIC <br>**Authors** Wen Shi
# MAGIC
# MAGIC **Created on** 2023.12.
# MAGIC
# MAGIC **Last updated on** 2024.10.25
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard
# MAGIC
# MAGIC **Reviewed** 2024.11.17
# MAGIC
# MAGIC **Data inputs** \
# MAGIC functions - libraries - parameters\
# MAGIC gdppr\
# MAGIC ccu073_01_cur_region_imd_lookup\
# MAGIC ccu073_01_cur_lsoa_imd_lookup
# MAGIC
# MAGIC **Data outputs**
# MAGIC - **`ccu073_01_out_lsoa_year2020`**
# MAGIC - **`ccu073_01_out_lsoa_year2022`** 
# MAGIC
# MAGIC **Acknowledgements** BHF datq science centre
# MAGIC
# MAGIC **Notes**
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC # 0. Setup

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
from datetime import date

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

lsoa_mult = spark.table(path_tmp_lsoa)

# COMMAND ----------

# -----------------------------------------------------------------------------
# LSOA Curated Data
# -----------------------------------------------------------------------------
lsoa_region = spark.table(path_cur_lsoa_region)
lsoa_imd    = spark.table(path_cur_lsoa_imd)



# COMMAND ----------

lsoa_region.printSchema()

# COMMAND ----------

lsoa_imd.printSchema()

# COMMAND ----------

# display(lsoa_region)

# COMMAND ----------

# display(lsoa_imd)

# COMMAND ----------

# MAGIC %md
# MAGIC #3 Apply rules to conflicted nearest LSOA

# COMMAND ----------

for i in cohort:
    tmp = (lsoa_mult
           .withColumnRenamed('person_id','PERSON_ID')
           .withColumn('lsoa_diff',f.datediff(f.col('record_date') ,f.lit(date(i,1,1))))
            .withColumn('lsoa_diff_abs',f.abs('lsoa_diff'))
            .withColumn('lsoa_record_window',f.when(f.col('lsoa_diff')<=0,1).when((f.col('lsoa_diff') >0) & (f.col('lsoa_diff') <=365),2).otherwise(3))
            .distinct()
        )
    _win =(Window
            .partitionBy('PERSON_ID')
            .orderBy('lsoa_record_window','lsoa_diff_abs')
            )
    tmp = tmp.withColumn('_rank',f.dense_rank().over(_win)).where('_rank ==1').drop('_rank')
    tmp = (tmp.pandas_api().groupby('PERSON_ID').agg('mode').reset_index()
            .to_spark().withColumn('_rownum', f.row_number().over(Window.partitionBy('PERSON_ID').orderBy('lsoa_diff_abs')))
            .where('_rownum = 1').drop('_rownum')
            
          )
    
    print('='*20,'Year', i,'='*20)
    # count_var(tmp,'PERSON_ID');print()
    # tab(tmp,'lsoa_record_window');print()
    
    outName = f'{proj}_tmp_rank_mode_lsoa_year{i}'.lower()
    tmp.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')

# COMMAND ----------

# MAGIC %md # 4. Add region & imd

# COMMAND ----------

for i in cohort:
  tmp0 = spark.table(f'{path_tmp_rank_lsoa}{i}')
  tmp0 = (tmp0
        .withColumn('LSOA_1',f.substring(f.col('lsoa'),1,1))
        )
  tmp1 = merge(tmp0, lsoa_region, ['lsoa'], validate='m:1', keep_results=['both', 'left_only'])
        

  tmp2 = (
    tmp1
    .withColumn('region',
                f.when(f.col('LSOA_1') == 'W', 'Wales')
                .when(f.col('LSOA_1') == 'S', 'Scotland')
                .otherwise(f.col('region'))
              )
    .drop('_merge')
  )


  tmp3 = merge(tmp2.withColumnRenamed('lsoa','LSOA'), lsoa_imd, ['LSOA'], validate='m:1', keep_results=['both', 'left_only'])
  assert tmp3.count()==tmp3.select('PERSON_ID').distinct().count()

  outName = f'{proj}_out_lsoa_year{i}'.lower()
  tmp3.select('PERSON_ID','region','IMD_2019_DECILES').write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')

  # check
#   tmp = spark.table(f'{dsa}.{outName}')
#   tab(tmp,'region');print()
#   tab(tmp,'IMD_2019_DECILES');print()


 

