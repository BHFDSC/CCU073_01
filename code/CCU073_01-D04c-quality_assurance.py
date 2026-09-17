# Databricks notebook source
# MAGIC %md # CCU073_01-D04c-quality_assurance
# MAGIC  
# MAGIC **Description** This notebook creates the quality assurance table, which indentifies person identifiers to remove from the analyses due to conflicting data, with reference to previous work/coding by CCU002 and Spiros Denaxas.
# MAGIC  
# MAGIC **Authors** Wen Shi
# MAGIC
# MAGIC **Created on** 2023.11.17
# MAGIC
# MAGIC **Last updated on** 2024.10.28 
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC codelist_qa, skinny, death, gdppr
# MAGIC
# MAGIC **Data output** \
# MAGIC CCU073_01_quality_assurance
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard
# MAGIC
# MAGIC **Last Reviewed** 2024.11.17
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

# MAGIC %md 
# MAGIC # 0 Parameters

# COMMAND ----------

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

# MAGIC %md # 1 Data

# COMMAND ----------

codelist_qa = spark.table(path_out_codelist_quality_assurance)
skinny      = spark.table(path_out_skinny)
death_sing      = spark.table(path_cur_deaths_sing)
gdppr       = extract_batch_from_archive(parameters_df_datasets, 'gdppr')

# COMMAND ----------

# MAGIC %md # 2 Prepare

# COMMAND ----------


_skinny = skinny.select('PERSON_ID', 'DOB', 'SEX','DOD','death_flag')

_death = death_sing.select('PERSON_ID','REG_DATE')

_merge = _skinny.join(_death,'PERSON_ID','left')

_gdppr = (
  gdppr
  .select(f.col('NHS_NUMBER_DEID').alias('PERSON_ID'), 'DATE', 'RECORD_DATE', 'CODE')
  .where(f.col('PERSON_ID').isNotNull()))
  

# COMMAND ----------

# MAGIC %md # 3 Medical conditions

# COMMAND ----------

# MAGIC %md ## 3.1 Codelist

# COMMAND ----------

# check
tmpt = tab(codelist_qa, 'name', 'terminology', var2_unstyled=1); print()

# check
_list_terms = list(
  codelist_qa
  .select('terminology')
  .distinct()
  .toPandas()['terminology'])
assert set(_list_terms) <= set(['SNOMED'])

# COMMAND ----------

# MAGIC %md ## 3.2 Create

# COMMAND ----------

# prepare
_gdppr_hx = (
  _gdppr
  .withColumn('DATE', f.when(f.col('DATE').isNull(), f.col('RECORD_DATE')).otherwise(f.col('DATE')))  
  .where(f.col('DATE').isNotNull())
  .withColumn('CENSOR_DATE_START', f.lit(None))
  .withColumn('CENSOR_DATE_END', f.lit(None)))       

# dictionary - dataset, codelist, and ordering in the event of tied records
_hx_in = {
    'gdppr':   ['_gdppr_hx',   'codelist_qa', 1]
}

# run codelist match and codelist match summary functions
_hx, _hx_1st, _hx_1st_wide = codelist_match(_hx_in, _name_prefix=f'hx_')
_hx_summ_name, _hx_summ_name_code = codelist_match_summ(_hx_in, _hx)

# temp save
_hx_1st = _hx_1st.drop('CENSOR_DATE_START', 'CENSOR_DATE_END')
outName = f'{proj}_tmp_quality_assurance_hx_1st'
_hx_1st.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
_hx_1st = spark.table(f'{dsa}.{outName}')

# temp save
_hx_1st_wide = _hx_1st_wide.drop('CENSOR_DATE_START', 'CENSOR_DATE_END')
outName = f'{proj}_tmp_quality_assurance_hx_1st_wide'
_hx_1st_wide.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
_hx_1st_wide = spark.table(f'{dsa}.{outName}')

# COMMAND ----------

# MAGIC %md ## 3.3 Check

# COMMAND ----------

# tmpt = tab(_hx_1st, 'name'); print()
# tmpt = tab(_hx_1st_wide, 'hx_pregnancy_flag', 'hx_prostate_cancer_flag', var2_unstyled=1); print()



# COMMAND ----------

# MAGIC %md ## 3.4 Merge

# COMMAND ----------

_merged = merge(_merge, _hx_1st_wide.select('PERSON_ID', 'hx_pregnancy_flag', 'hx_prostate_cancer_flag'), ['PERSON_ID'], validate='1:1', assert_results=['both', 'left_only'], indicator=0); print()

# COMMAND ----------

# check
tmpt = tab(_merged, 'SEX', 'hx_pregnancy_flag', var2_unstyled=1); print()
tmpt = tab(_merged, 'SEX', 'hx_prostate_cancer_flag', var2_unstyled=1); print()

# COMMAND ----------

# MAGIC %md # 4 Rules

# COMMAND ----------

# MAGIC %md ## 4.1 Prepare

# COMMAND ----------

# ------------------------------------------------------------------------------
# preparation: rule 8 (Patients have all missing record_dates and dates)
# ------------------------------------------------------------------------------
# identify records with null date
_gdppr_null = _gdppr\
  .select('PERSON_ID', 'DATE', 'RECORD_DATE')\
  .withColumn('_null',\
    f.when(\
      ((f.col('DATE').isNotNull()) | (f.col('DATE') == '') | (f.col('DATE') == ' '))\
      | ((f.col('RECORD_DATE').isNotNull()) | (f.col('RECORD_DATE') == '') | (f.col('RECORD_DATE') == ' '))\
    , 0).otherwise(1)\
  )\

# check
tmpt = tab(_gdppr_null, '_null'); print()

# summarise per individual
_gdppr_null_summ = _gdppr_null\
  .groupBy('PERSON_ID')\
  .agg(\
    f.sum(f.when(f.col('_null') == 0, 1).otherwise(0)).alias('_n_gdppr_notnull')\
    , f.sum(f.col('_null')).alias('_n_gdppr_null')\
  )\
  .where(f.col('_n_gdppr_null') > 0)

# cache
_gdppr_null_summ.cache().count()

# check
print(_gdppr_null_summ.toPandas().to_string()); print()

# check
tmp = _gdppr_null_summ\
  .select('_n_gdppr_null')\
  .groupBy()\
  .sum()\
  .collect()[0][0]
print(tmp); print()

# merge
_merged = merge(_merged, _gdppr_null_summ, ['PERSON_ID'], validate='1:1', assert_results=['both', 'left_only'], indicator=0); 

# COMMAND ----------

# MAGIC %md ## 4.2 Create

# COMMAND ----------

# Rule 1: Year of birth is after the year of death
# Rule 2: Patient does not have mandatory fields completed (nhs_number, sex, Date of birth)
# Rule 3: Year of Birth Predates NHS Established Year or Year is over the Current Date
# Rule 4: Remove those with only null/invalid dates of death
# Rule 5: Remove those where registered date of death before the actual date of death
# Rule 6: Pregnancy/birth codes for men
# Rule 7: Prostate Cancer Codes for women
# Rule 8: Patients have all missing record_dates and dates

_qax = _merged\
  .withColumn('YOB', f.year(f.col('DOB')))\
  .withColumn('YOD', f.year(f.col('DOD')))\
  .withColumn('_rule_1', f.when(f.col('YOB') > f.col('YOD'), 1).otherwise(0))\
  .withColumn('_rule_2',\
    f.when(\
      (f.col('SEX').isNull()) | (~f.col('SEX').isin([1,2]))\
      | (f.col('DOB').isNull())\
      | (f.col('PERSON_ID').isNull()) | (f.col('PERSON_ID') == '') | (f.col('PERSON_ID') == ' ')\
    , 1).otherwise(0)\
  )\
  .withColumn('_rule_3',\
    f.when(\
      (f.col('YOB') < 1793) | (f.col('YOB') > datetime.datetime.today().year)\
    , 1).otherwise(0)\
  )\
  .withColumn('_rule_4',\
    f.when(\
      (f.col('death_flag') == 1)\
      & (\
        (f.col('DOD').isNull())\
        | (f.col('DOD') <= f.to_date(f.lit('1900-01-01')))\
        | (f.col('DOD') > f.current_date())\
      )\
    , 1).otherwise(0)\
  )\
  .withColumn('_rule_5', f.when(f.col('DOD') > f.col('REG_DATE'), 1).otherwise(0))\
  .withColumn('_rule_6', f.when((f.col('SEX') == 1) & (f.col('hx_pregnancy_flag') == 1) , 1).otherwise(0))\
  .withColumn('_rule_7', f.when((f.col('SEX') == 2) & (f.col('hx_prostate_cancer_flag') == 1) , 1).otherwise(0))\
  .withColumn('_rule_8', f.when((f.col('_n_gdppr_null') > 0) & (f.col('_n_gdppr_notnull') == 0) , 1).otherwise(0))
  

# row total and concat
_qax = _qax\
  .withColumn('_rule_total', sum([f.col(col) for col in _qax.columns if re.match('^_rule_.*$', col)]))\
  .withColumn('_rule_concat', f.concat(*[f'_rule_{i}' for i in list(range(1, 9))]))

# temp save
outName = f'{proj}_tmp_quality_assurance_qax'
_qax.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
_qax = spark.table(f'{dsa}.{outName}')

# COMMAND ----------

# MAGIC %md ## 4.3 Checks

# COMMAND ----------

# check
count_var(_qax, 'PERSON_ID'); print()
tmpt = tab(_qax, '_rule_total'); print()
tmpt = tab(_qax, '_rule_concat'); print()

# check rule frequency
for i in list(range(1, 9)):
  tmpt = tab(_qax, f'_rule_{i}'); print()


# COMMAND ----------

# MAGIC %md # 5 Save 

# COMMAND ----------

# save name
outName = f'{proj}_quality_assurance'.lower()

# save
_qax.select('PERSON_ID','_rule_total').write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
