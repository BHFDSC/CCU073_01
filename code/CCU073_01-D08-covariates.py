# Databricks notebook source
# MAGIC %md # CCU073_01-D08_covariates
# MAGIC
# MAGIC **Description** This notebook creates the covariates needed for CCU073_01:
# MAGIC - Smoking
# MAGIC - BMI
# MAGIC - Blood pressure (Systolic)
# MAGIC - Total cholesterol
# MAGIC - HDL-cholesterol
# MAGIC - eGFR
# MAGIC - Creatinine
# MAGIC - HbA1c
# MAGIC - Glucose level
# MAGIC
# MAGIC We select the last record available within 5 years lookback from Baseline.
# MAGIC However, smoking status is captured at the latest records ever before baseline.
# MAGIC
# MAGIC **Author(s)** Genevieve Cezard, Tom Bolton
# MAGIC
# MAGIC **Created on** 2024.02.05
# MAGIC
# MAGIC **Last updated on** 2024.03.26
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC codelist_covariates\
# MAGIC ccu073_01_out_cohort_year2020\
# MAGIC ccu073_01_out_cohort_year2022\
# MAGIC gdppr
# MAGIC
# MAGIC **Data output** \
# MAGIC CCU073_01_out_covariates
# MAGIC
# MAGIC **Reviewers** Wen Shi, Tom Bolton
# MAGIC
# MAGIC **Reviewed** 2024.10.31
# MAGIC
# MAGIC
# MAGIC **Acknowledgements** Adapted from previous work by Tom Bolton from CCU004_01 and CCU051.
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

#import the pyspaprk module
# from pyspark.sql.functions import col,lit,when

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

# MAGIC %md # 0. Parameters

# COMMAND ----------

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

# MAGIC %md # 1. Data

# COMMAND ----------

# codelists
codelist_covariates = spark.table(path_out_codelist_covariates)

# GDPPR
gdppr = extract_batch_from_archive(parameters_df_datasets, 'gdppr') 


# COMMAND ----------

# GDPPR curated with mono_id - created once only
cur_gdppr = (gdppr
            .select(f.col('NHS_NUMBER_DEID').alias('PERSON_ID'), 'DATE', 'CODE', 'VALUE1_CONDITION', 'VALUE2_CONDITION')
            .withColumn('mono_id', f.monotonically_increasing_id())
)

outName = f'{proj}_cur_gdppr'
cur_gdppr.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
cur_gdppr = spark.table(f'{dsa}.{outName}')
# display(cur_gdppr)

# COMMAND ----------

# MAGIC %md # 2. Prepare

# COMMAND ----------

# MAGIC %md ## 2.1. List of covariates - Values 'out of range' 

# COMMAND ----------

# List covariates
covariates_smoking = ['smoking_current', 'smoking_ex', 'smoking_never']
covariates = ['bmi', 'tchol', 'hdl', 'hba1c', 'Glucose_level', 'egfr', 'Creatinine', 'sbp']

# Treat 'dbp' separately as we need to use a combinaison of value 1 and value2

# COMMAND ----------

codelist_range = """
name,lower_bound,upper_bound
bmi,10,100
tchol, 1.75,20
hdl,0.2,10
hba1c,0,195
egfr,0,175
Creatinine,0,250
sbp,60,250
dbp,30,250
Glucose_level,1,100
"""

codelist_range = (spark.createDataFrame(
  pd.DataFrame(pd.read_csv(io.StringIO(codelist_range)))
  .fillna('')
  #.astype(str)
))
#display(codelist_range)

# COMMAND ----------

# MAGIC %md ## 3.2. Select covariates last date and value

# COMMAND ----------

win_row_number = Window.partitionBy('PERSON_ID', 'name').orderBy(f.desc('DATE'), 'mono_id')

for i in cohort:

    print(f"{'year'+str(i):-^50}")

    # 1. Get cohort
    spark.sql(f"""REFRESH TABLE {path_out_inc_exc_cohort}{i}""")
    cohort = spark.table(f'{path_out_inc_exc_cohort}{i}')
    tmp_cohort = (cohort.select(['PERSON_ID', 'DOB','baseline_date']).withColumn('CENSOR_DATE_START', f.add_months(f.col('baseline_date'), -12*5)))

    # 2. Select the last Smoking record of each type (current, ex, never) before baseline and create smoking status.
    tmp_gdppr = (cur_gdppr
                .select('PERSON_ID', 'DATE', 'CODE')
                .join(tmp_cohort, on=['PERSON_ID'], how='inner')
                .where((f.col('DATE') >= f.col('DOB')) & (f.col('DATE')<=f.col('baseline_date')))
                )

    tmp_covariates_smoking = tmp_cohort.select('PERSON_ID','DOB')
    for j in covariates_smoking:
        codelist_smoking = codelist_covariates.filter(f'''name like "{j}"''').where('terminology="SNOMED"').select('name','code').distinct()
        tmp_gdppr_code_smoking = (tmp_gdppr
                                  .select('PERSON_ID', 'DATE', 'CODE')
                                  .join(codelist_smoking, on=['code'], how='inner')
                                  .dropDuplicates(['PERSON_ID', 'name', 'DATE'])
                                  .groupBy('PERSON_ID').agg(f.max('DATE').alias(f'cov_{j}_date'))
                                  )
        tmp_covariates_smoking = (tmp_covariates_smoking.join(tmp_gdppr_code_smoking,'PERSON_ID','left'))

    tmp_covariates = (tmp_covariates_smoking
                    .select('PERSON_ID','cov_smoking_current_date','cov_smoking_ex_date','cov_smoking_never_date')
                    .withColumn('Smoking_status',
                                f.when((f.col('cov_smoking_current_date').isNotNull()) & (f.col('cov_smoking_ex_date').isNull()) & (f.col('cov_smoking_never_date').isNull()),'Current')
                                .when((f.col('cov_smoking_current_date')>=f.col('cov_smoking_ex_date')),'Current')
                                .when((f.col('cov_smoking_ex_date').isNull()) & (f.col('cov_smoking_current_date')>=f.col('cov_smoking_never_date')),'Current')
                                .when((f.col('cov_smoking_current_date')>=f.col('cov_smoking_ex_date')) & (f.col('cov_smoking_current_date')>=f.col('cov_smoking_never_date')),'Current')
                                .when((f.col('cov_smoking_current_date')<f.col('cov_smoking_ex_date')),'Ex')
                                .when(f.col('cov_smoking_current_date')<f.col('cov_smoking_never_date'),'Ex')
                                .when((f.col('cov_smoking_current_date').isNull()) & (f.col('cov_smoking_ex_date').isNotNull()) ,'Ex')
                                .when((f.col('cov_smoking_current_date').isNull()) & (f.col('cov_smoking_ex_date').isNull()) & (f.col('cov_smoking_never_date').isNotNull()),'Never')
                                .when((f.col('cov_smoking_current_date').isNull()) & (f.col('cov_smoking_ex_date').isNull()) & (f.col('cov_smoking_never_date').isNull()),'Missing')
                    )
    )  

    # 3. Merge gdppr with cohort and apply date range restriction for other covariates
    tmp_gdppr = (cur_gdppr
                 .select('PERSON_ID', 'DATE', 'CODE', 'VALUE1_CONDITION', 'VALUE2_CONDITION', 'mono_id')
                 .join(tmp_cohort, on=['PERSON_ID'], how='inner')
                 .where((f.col('DATE') >= f.col('CENSOR_DATE_START')) & (f.col('DATE')<=f.col('baseline_date')))
                 .where(f.col('VALUE1_CONDITION').isNotNull())
                 )

    #4. Select last date and values (in range) for each covariate
    for j in covariates: 
        covcode = codelist_covariates.filter(f'''name like "{j}"''').where('terminology="SNOMED"').select('name','code').distinct()
        tmp_gdppr_code = (tmp_gdppr
                        .join(covcode, on=['code'], how='inner')
                        .join(codelist_range.select('name', 'lower_bound', 'upper_bound'), on=['name'], how='left')
                        .where((f.col('lower_bound').isNull()) | (f.col('VALUE1_CONDITION') >= f.col('lower_bound')))
                        .where((f.col('upper_bound').isNull()) | (f.col('VALUE1_CONDITION') <= f.col('upper_bound')))
                        .dropDuplicates(['PERSON_ID', 'name', 'DATE', 'VALUE1_CONDITION'])
                        .withColumn('VALUE1_CONDITION', f.round(f.col('VALUE1_CONDITION'), 2))
                        .withColumn('rownum', f.row_number().over(win_row_number))
                        .where(f.col('rownum') == 1)
                        .select('PERSON_ID', f.col('DATE').alias(f'cov_{j}_date'), f.col('VALUE1_CONDITION').alias(f'cov_{j}_value'))
        )
        tmp_covariates = (tmp_covariates.join(tmp_gdppr_code,'PERSON_ID','left'))

    #5. Select either value 1 or value 2 for dbp and take last date and values (in range)
    codelist_dbp = codelist_covariates.filter(f'''name like "dbp"''').where('terminology="SNOMED"').select('name','code').distinct()
    tmp_gdppr_code = (tmp_gdppr
                    .join(codelist_dbp, on=['code'], how='inner')
                    .join(codelist_range.select('name', 'lower_bound', 'upper_bound'), on=['name'], how='left')
                    .withColumn('value',
                                f.when((f.col('VALUE1_CONDITION').isNotNull()) & (f.col('VALUE2_CONDITION').isNull()),f.round(f.col('VALUE1_CONDITION'), 2))
                                .when((f.col('VALUE1_CONDITION').isNotNull()) & (f.col('VALUE2_CONDITION').isNotNull()) & (f.col('VALUE1_CONDITION')<=f.col('VALUE2_CONDITION')), f.round(f.col('VALUE1_CONDITION'), 2))
                                .when((f.col('VALUE1_CONDITION').isNotNull()) & (f.col('VALUE2_CONDITION').isNotNull()) & (f.col('VALUE1_CONDITION')>f.col('VALUE2_CONDITION')), f.round(f.col('VALUE2_CONDITION'), 2))
                    )
                    .where((f.col('lower_bound').isNull()) | (f.col('value') >= f.col('lower_bound')))
                    .where((f.col('upper_bound').isNull()) | (f.col('value') <= f.col('upper_bound')))
                    .dropDuplicates(['PERSON_ID', 'name', 'DATE', 'value'])
                    .withColumn('rownum', f.row_number().over(win_row_number))
                    .where(f.col('rownum') == 1)
                    .select('PERSON_ID', f.col('DATE').alias(f'cov_dbp_date'), f.col('value').alias(f'cov_dbp_value'))
    )
    tmp_covariates = (tmp_covariates.join(tmp_gdppr_code,'PERSON_ID','left'))

    tmp_covariates.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{proj}_out_covariates_year{i}')
    tmp_covariates = spark.table(f'{dsa}.{proj}_out_covariates_year{i}')

# COMMAND ----------

# MAGIC %md # 3. Check

# COMMAND ----------

#display(tmp_covariates)

# COMMAND ----------

#tab(tmp_covariates,'Smoking_status')
