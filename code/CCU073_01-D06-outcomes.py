# Databricks notebook source
# MAGIC %md # CCU073_01-D06-outcomes
# MAGIC
# MAGIC **Description** This notebook adds flag and extracts ealiest date for each outcome for each person.
# MAGIC
# MAGIC **Authors** Wen Shi
# MAGIC
# MAGIC **Created on** 2023.12.17
# MAGIC
# MAGIC **Last updated on** 2024.11.12
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC Codelist outcomes\
# MAGIC ccu073_01_cur_hes_apc_all_years_long\
# MAGIC ccu073_01_cur_deaths_long\
# MAGIC ccu073_01_out_cohort_year2020\
# MAGIC ccu073_01_out_cohort_year2022\
# MAGIC gdppr
# MAGIC
# MAGIC **Data outputs**
# MAGIC - **`ccu073_01_out_outcome_year2020`**
# MAGIC - **`ccu073_01_out_outcome_year2022`**
# MAGIC
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard
# MAGIC
# MAGIC **Reviewed on** 2024.01.23/ 2024.04.25/ 2024.11.17
# MAGIC
# MAGIC **Acknowledgements** 
# MAGIC
# MAGIC **Notes**

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

spark.sql(f"""REFRESH TABLE {path_out_codelist_outcome}""")
codelist     = spark.table(path_out_codelist_outcome)
spark.sql(f"""REFRESH TABLE {path_cur_hes_apc_long}""")
hes_apc_long = spark.table(path_cur_hes_apc_long)
spark.sql(f"""REFRESH TABLE {path_cur_deaths_long}""")
deaths_long  = spark.table(path_cur_deaths_long)
gdppr   = extract_batch_from_archive(parameters_df_datasets, 'gdppr')

# COMMAND ----------

gdppr = gdppr.where((f.col('CODE').isNotNull()) | (f.col('CODE') != '')).select(f.col('NHS_NUMBER_DEID').alias('PERSON_ID'),f.coalesce('DATE','RECORD_DATE').alias('DATE'),'CODE')

# COMMAND ----------

# MAGIC %md # 2 Prepare

# COMMAND ----------

inci_out = ['Angina','CHD','ARR','HF','stroke','PVD','DVT']
death_out = ['Hypertensive_disease','CHD','OtherHD','stroke','Arter','DVT','Sudden_death','OtherCVD','Angina','ARR','HF','PVD']
hist_out = inci_out+death_out
# hist_out.remove('Sudden_death')
hist_out = set(hist_out)

# COMMAND ----------

for j in cohort:
    print(f'{j:-^80}')
    spark.sql(f"""REFRESH TABLE {path_out_inc_exc_cohort}{j}""")
    tmp_cohort = spark.table(f'{path_out_inc_exc_cohort}{j}')
    k = 0
    for i in hist_out:
        outcode = codelist.filter((f.col('name')==i)&(f.col('terminology')=="ICD10")).select('code','code_type').distinct()
        outcode_sno = codelist.filter((f.col('name')==i)&(f.col('terminology')=="SNOMED")).select('code','code_type').distinct()
        

        _hes_apc = hes_apc_long.join(f.broadcast(outcode),'code')
        _hes_apc = _hes_apc.select(['PERSON_ID', 'EPISTART','DIAG_POSITION','code_type']).withColumnRenamed('EPISTART', 'DATE')

        _death = deaths_long.join(f.broadcast(outcode),'code').where(f.col('DIAG_POSITION') == 'UNDERLYING').select('PERSON_ID','DATE')

        _gdppr = gdppr.join(f.broadcast(outcode_sno),'code').select('PERSON_ID','DATE','code_type')

        tmp1 = (
        _gdppr
        .withColumn('DIAG_POSITION', f.lit(None))
        .withColumn('source', f.lit('gdppr'))
        .unionByName(_hes_apc.withColumn('source', f.lit('hes_apc')))
        .unionByName(_death.withColumn('DIAG_POSITION', f.lit(None))
                        .withColumn('code_type',f.lit(None))
                        .withColumn('source', f.lit('deaths'))))
        tmp1.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{proj}_tmp1_{i}')
        tmp1 = spark.table(f'{dsa}.{proj}_tmp1_{i}')


        tmp2 = (tmp_cohort.join(tmp1,'PERSON_ID','left')
                .withColumn('DATE',f.when((f.col('DATE')<f.col('DOB'))|(f.col('DATE')>f.col('DOD')),None).otherwise(f.col('DATE'))))
        tmp2.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{proj}_tmp2_{i}_year{j}')
        tmp2 = spark.table(f'{dsa}.{proj}_tmp2_{i}_year{j}')

        tmp2 = (tmp2
                .withColumn(f'hist_{i}_date',f.when((f.col('DATE')<datetime.date(j,1,1))&(f.col('source').isin(['hes_apc','gdppr'])),f.col('DATE')).otherwise(None))
                .withColumn(f'nonfatal_{i}_date',
                        f.when((f.col('source').isin(['hes_apc','gdppr']))&(f.col('DATE')>=datetime.date(j,1,1))&(f.col('code_type').isin(['1'])),f.col('DATE')).otherwise(None))
                .withColumn(f'fatal_{i}_date',
                        f.when((f.col('source').isin(['deaths']))&(f.col('DATE')>=datetime.date(j,1,1)),f.col('DATE')).otherwise(None))
        )   
        tmp3 = (tmp2
                .groupBy('PERSON_ID')
                .agg(f.min(f.col(f'hist_{i}_date')).alias(f'hist_{i}_date'),
                        f.min(f.col(f'nonfatal_{i}_date')).alias(f'nonfatal_{i}_date'),
                        f.min(f.col(f'fatal_{i}_date')).alias(f'fatal_{i}_date')
                        )
        )
        tmp3.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{proj}_tmp3_{i}_year{j}')
        tmp3 = spark.table(f'{dsa}.{proj}_tmp3_{i}_year{j}')
        print(f'outcome {i} finished')
        if k==0: tmpc = tmp3
        else: tmpc = tmpc.join(tmp3,'PERSON_ID')
        k = k + 1
    tmpc.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{proj}_out_outcome_year{j}')


# COMMAND ----------

# MAGIC %md
# MAGIC #3 Check

# COMMAND ----------

# out1 = spark.table(f'{dsa}.{proj}_out_outcome_v2_year2020')
# for i in inci_out:
    
#         print(f"{i}: {out1.where(f.col(f'nonfatal_{i}_date_hesgp_pos12').isNotNull()).count()}")

# COMMAND ----------

out2 = spark.table(f'{dsa}.{proj}_out_outcome_year2020')
print(f'{"HES pos all":=^80}')
for i in inci_out:
    print(f"{i}: {out2.where(f.col(f'nonfatal_{i}_date').isNotNull()).count()}")


# COMMAND ----------

out2_2 =  spark.table(f'{dsa}.{proj}_out_outcome_year2022')
print(f'{"HES pos all":=^80}')
for i in inci_out:
    print(f"{i}: {out2_2.where(f.col(f'nonfatal_{i}_date').isNotNull()).count()}")

# COMMAND ----------

# display(out2)
