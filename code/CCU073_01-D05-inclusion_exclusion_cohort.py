# Databricks notebook source
# MAGIC %md # CCU073_01-D05-inclusion_exclusion
# MAGIC  
# MAGIC **Description** This notebook creates the cohort for baseline 2020/1/1 and 2022/1/1 after applying inclusion exclusion criteria.
# MAGIC  
# MAGIC **Authors** Wen Shi
# MAGIC
# MAGIC **Created on** 2023.12.
# MAGIC
# MAGIC **Last updated on** 2024.10.28
# MAGIC  
# MAGIC **Reviewers** Genevieve Cezard
# MAGIC
# MAGIC **Reviewed** 2024.11.17
# MAGIC
# MAGIC **Data inputs** \
# MAGIC functions - libraries - parameters\
# MAGIC ccu073_01_skinny\
# MAGIC ccu073_01_quality_assurance\
# MAGIC ccu073_01_out_lsoa_year2020\
# MAGIC ccu073_01_out_lsoa_year2022
# MAGIC
# MAGIC **Data outputs**
# MAGIC - **`ccu073_01_out_cohort_year2020`**
# MAGIC - **`ccu073_01_out_cohort_year2022`** 
# MAGIC - **`ccu073_01_out_flow_inc_exc_year2020`**
# MAGIC - **`ccu073_01_out_flow_inc_exc_year2022`** 
# MAGIC
# MAGIC **Acknowledgements** Based on previous work by Tom Bolton (John Nolan, Elena Raffetti) for CCU018_01 and the earlier CCU002 sub-projects.
# MAGIC
# MAGIC **Notes**
# MAGIC Patients will be excluded if they meet ANY the following criteria:
# MAGIC
# MAGIC 1.Exclusion of patients not in GDPPR
# MAGIC <br>2.Exclusion of patients who died before baseline
# MAGIC <br>3.Exclusion of patients with region missing or outside England
# MAGIC <br>4.Exclusion of patients without deprivation
# MAGIC <br>5.Exclusion of patients who failed the quality assurance
# MAGIC <br>6.Exclusion of patients who were aged under 40 years old at baseline date

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

spark.sql(f"""REFRESH TABLE {dsa}.{proj}_skinny""")
spark.sql(f"""REFRESH TABLE {dsa}.{proj}_quality_assurance""")


skinny  = spark.table(f'{dsa}.{proj}_skinny')
qa      = spark.table(f'{dsa}.{proj}_quality_assurance')


# COMMAND ----------

# MAGIC %md # 2 Prepare

# COMMAND ----------

print('---------------------------------------------------------------------------------')
print('skinny')
print('---------------------------------------------------------------------------------')
# reduce
_skinny = skinny.select('PERSON_ID', 'DOB', 'SEX', 'ETHNIC_CAT','DOD', 'in_gdppr')


# check
# count_var(_skinny, 'PERSON_ID'); print()


print('---------------------------------------------------------------------------------')
print('quality assurance')
print('---------------------------------------------------------------------------------')
_qa = qa

# check
# count_var(_qa, 'PERSON_ID'); print()
# tmpt = tab(_qa, '_rule_total'); print()


print('---------------------------------------------------------------------------------')
print('merged')
print('---------------------------------------------------------------------------------')
_merged = (
  merge(_skinny, _qa, ['PERSON_ID'], validate='1:1', assert_results=['both'])
  .drop('_merge')); print()

# COMMAND ----------

# MAGIC %md # 3 Inclusion / exclusion

# COMMAND ----------

for i in cohort:
    spark.sql(f"""REFRESH TABLE {dsa}.{proj}_out_lsoa_year{i}""")
    cov_lsoa   = spark.table(f'{dsa}.{proj}_out_lsoa_year{i}')
    
    # merge with lsoa

    _merged_with_lsoa = _merged.join(cov_lsoa,'PERSON_ID','left')

    # add baseline date
    _merged_with_lsoa = (_merged_with_lsoa.withColumn('baseline_date', f.when(f.col('DOB') > f'{i}-01-01', f.col('DOB')).otherwise(f.to_date(f.lit(f'{i}-01-01')))))

    # temp save
    #outName = f'{proj}_tmp_inc_exc_merged_year{i}'.lower()
    #_merged_with_lsoa.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
    #_merged_with_lsoa = spark.table(f'{dsa}.{outName}')
    tmp0 = _merged_with_lsoa
    tmpc = count_var(tmp0, 'PERSON_ID', ret=1, df_desc='original', indx=0)
    # 3.1 Exclude patients not in GDPPR
    tmp1 = tmp0.where(f.col('in_gdppr') == 1)
    tmpt = count_var(tmp1, 'PERSON_ID', ret=1, df_desc='post exclusion of patients not in GDPPR', indx=1)
    tmpc = tmpc.unionByName(tmpt)
    # 3.2 Exclude patients who died before baseline
    tmp1 = (tmp1
        .withColumn('flag_DOD_lt_baseline', f.when(f.col('DOD') < f.col('baseline_date'), 1).otherwise(0)))
    tmp2 = (tmp1
        .where(f.col('flag_DOD_lt_baseline') == 0)
        .drop('flag_DOD_lt_baseline'))
    tmpt = count_var(tmp2, 'PERSON_ID', ret=1, df_desc='post exlusion of patients who died before baseline', indx=2)
    tmpc = tmpc.unionByName(tmpt)
    # 3.3 Exclude patients with region missing or outside England
    tmp3 = tmp2.where(~f.col('region').isin(['Scotland','Wales']))
    tmpt = count_var(tmp3, 'PERSON_ID', ret=1, df_desc='post exclusion of patients with region outside England', indx=3)
    tmpc = tmpc.unionByName(tmpt)
    # 3.4 Exclude patients without deprivation
    tmp4 = tmp3.where(~f.col('IMD_2019_DECILES').isNull())
    tmpt = count_var(tmp4, 'PERSON_ID', ret=1, df_desc='post exclusion of patients without deprivation', indx=4)
    tmpc = tmpc.unionByName(tmpt)
    # 3.5 Exclude patients who failed the quality assurance
    tmp4 = (tmp4
        .withColumn('flag_rule_total_gt_0', f.when(f.col('_rule_total') > 0, 1).otherwise(0)))
    tmp5 = (tmp4
        .where(f.col('flag_rule_total_gt_0') == 0)
        .drop('flag_rule_total_gt_0'))
    tmpt = count_var(tmp5, 'PERSON_ID', ret=1, df_desc='post exclusion of patients who failed the quality assurance', indx=5)
    tmpc = tmpc.unionByName(tmpt)
    # 3.6 Exclude patients aged under 40 at baseline
    tmp5 = (tmp5.withColumn('flag_age_under_40',f.when(f.datediff(f.col('baseline_date'), f.col('DOB')) < 40*365.25 ,1).otherwise(0)))
    tmp6 = tmp5.where('flag_age_under_40 = 0').drop('flag_age_under_40')
    tmpt = count_var(tmp6, 'PERSON_ID', ret=1, df_desc='post exclusion of patients who were aged under 40 at baseline', indx=6)
    tmpc = tmpc.unionByName(tmpt)

    # check flow table
    tmpp = (tmpc
    .orderBy('indx')
    .select('indx', 'df_desc', 'n', 'n_id', 'n_id_distinct')
    .withColumnRenamed('df_desc', 'stage')
    .toPandas())
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
    print(f"{'year'+str(i):-^50}")
    print(tmpp.to_string()); print()
    tmpf = tmp6.select('PERSON_ID', 'DOB', 'SEX','ETHNIC_CAT', 'DOD','region', 'IMD_2019_DECILES', 'baseline_date')
    # check final cohort table contains one row per person
    assert tmpf.count() ==tmpf.distinct().count()
    assert tmpf.where('DOD <DOB').count() ==0
    outName = f'{proj}_out_cohort_year{i}'.lower()
    tmpf.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
    outName = f'{proj}_out_flow_inc_exc_year{i}'.lower()
    tmpc.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
