# Databricks notebook source
# MAGIC %md # CCU073_01_D03c_curated_data_covid
# MAGIC
# MAGIC
# MAGIC **Description** This notebook creates the covid phenotypes table of CCU073_01.
# MAGIC
# MAGIC **Author(s)** Tom Bolton, Fionna Chalmers, Anna Stevenson (Health Data Science Team, BHF Data Science Centre) + Genevieve Cezard (adapted from CCU004_01)
# MAGIC
# MAGIC **Created on** 2023.12.19
# MAGIC
# MAGIC **Last updated on** 2024.01.16
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC gdppr, hes apc, sus, sgss, chess
# MAGIC
# MAGIC **Data Output**
# MAGIC - **`ccu073_01_cur_covid`** : Covid phenotypes
# MAGIC
# MAGIC **Reviewers** Wen Shi
# MAGIC
# MAGIC **Reviewed** 2024.10.25
# MAGIC
# MAGIC **Acknowledgements** Based on CCU002_07 and subsequently CCU051_D03e_curated_data_covid and CCU004_01
# MAGIC

# COMMAND ----------

spark.sql('CLEAR CACHE')
spark.conf.set('spark.sql.legacy.allowCreatingManagedTableUsingNonemptyLocation', 'true')

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

# MAGIC %md # 0. Parameters

# COMMAND ----------

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

start_date = '2020-01-01'
end_date   = tmp_archived_on
print(start_date, end_date)

# COMMAND ----------

# MAGIC %md # 1. Data

# COMMAND ----------

codelist_covid = spark.table(path_out_codelist_covid)

# COMMAND ----------

sgss     = extract_batch_from_archive(parameters_df_datasets, 'sgss')
gdppr    = extract_batch_from_archive(parameters_df_datasets, 'gdppr')
hes_apc  = extract_batch_from_archive(parameters_df_datasets, 'hes_apc')
sus      = extract_batch_from_archive(parameters_df_datasets, 'sus')
chess    = extract_batch_from_archive(parameters_df_datasets, 'chess')


# COMMAND ----------

# MAGIC %md # 2. Prepare

# COMMAND ----------

# sgss
_sgss = sgss\
  .select(['PERSON_ID_DEID', 'Reporting_Lab_ID', 'Specimen_Date'])\
  .withColumnRenamed('PERSON_ID_DEID', 'PERSON_ID')\
  .withColumnRenamed('Specimen_Date', 'DATE')\
  .where((f.col('DATE') >= start_date) & (f.col('DATE') <= end_date))\
  .dropDuplicates()

# gdppr
# omitted: 'LSOA'
_gdppr = gdppr\
  .select(['NHS_NUMBER_DEID', 'DATE', 'CODE'])\
  .withColumnRenamed('NHS_NUMBER_DEID', 'PERSON_ID')\
  .where((f.col('DATE') >= start_date) & (f.col('DATE') <= end_date))\
  .dropDuplicates()

# hes_apc
# omitted: 'DISMETH', 'DISDEST', 'DISDATE', 'SUSRECID'
_hes_apc = hes_apc\
  .select(['PERSON_ID_DEID', 'EPISTART', 'DIAG_4_01', 'DIAG_4_CONCAT', 'OPERTN_4_CONCAT', 'SUSRECID'])\
  .withColumnRenamed('PERSON_ID_DEID', 'PERSON_ID')\
  .withColumnRenamed('EPISTART', 'DATE')\
  .where(f.col('DIAG_4_CONCAT').rlike('U07(1|2)'))\
  .where((f.col('DATE') >= start_date) & (f.col('DATE') <= end_date))\
  .dropDuplicates()


# sus
_sus = sus\
  .select(['NHS_NUMBER_DEID'
    , 'EPISODE_START_DATE'
    , 'PRIMARY_PROCEDURE_DATE'
    , 'SECONDARY_PROCEDURE_DATE_1'
    , 'DISCHARGE_DESTINATION_HOSPITAL_PROVIDER_SPELL'
    , 'DISCHARGE_METHOD_HOSPITAL_PROVIDER_SPELL'
    , 'END_DATE_HOSPITAL_PROVIDER_SPELL'           
    ]\
    + [col for col in sus.columns if re.match('.*(DIAGNOSIS|PROCEDURE)_CODE.*', col)]
  )\
  .withColumnRenamed('NHS_NUMBER_DEID', 'PERSON_ID')\
  .withColumnRenamed('EPISODE_START_DATE', 'DATE')\
  .withColumn('DIAG_CONCAT', f.concat_ws(',', *[col for col in sus.columns if re.match('.*DIAGNOSIS_CODE.*', col)]))\
  .withColumn('PROCEDURE_CONCAT', f.concat_ws(',', *[col for col in sus.columns if re.match('.*PROCEDURE_CODE.*', col)]))\
  .where((f.col('DATE') >= start_date) & (f.col('DATE') <= end_date))\
  .where(\
    ((f.col('END_DATE_HOSPITAL_PROVIDER_SPELL') >= start_date) | (f.col('END_DATE_HOSPITAL_PROVIDER_SPELL').isNull()))\
    & ((f.col('END_DATE_HOSPITAL_PROVIDER_SPELL') <= end_date) | (f.col('END_DATE_HOSPITAL_PROVIDER_SPELL').isNull()))\
  )\
  .where(f.col('DATE').isNotNull())\
  .dropDuplicates()

# chess
_chess = chess\
  .select(['PERSON_ID_DEID', 'Typeofspecimen', 'Covid19', 'AdmittedToICU', 'Highflownasaloxygen', 'NoninvasiveMechanicalventilation', 'Invasivemechanicalventilation', 'RespiratorySupportECMO', 'DateAdmittedICU', 'HospitalAdmissionDate', 'InfectionSwabDate'])\
  .withColumnRenamed('PERSON_ID_DEID', 'PERSON_ID')\
  .withColumnRenamed('InfectionSwabDate', 'DATE')\
  .where(f.col('Covid19') == 'Yes')\
  .where(\
    ((f.col('DATE') >= start_date) | (f.col('DATE').isNull()))\
    & ((f.col('DATE') <= end_date) | (f.col('DATE').isNull()))\
  )\
  .where(\
    ((f.col('HospitalAdmissionDate') >= start_date) | (f.col('HospitalAdmissionDate').isNull()))\
    & ((f.col('HospitalAdmissionDate') <= end_date) | (f.col('HospitalAdmissionDate').isNull()))\
  )\
  .where(\
    ((f.col('DateAdmittedICU') >= start_date) | (f.col('DateAdmittedICU').isNull()))\
    & ((f.col('DateAdmittedICU') <= end_date) | (f.col('DateAdmittedICU').isNull()))\
  )\
  .dropDuplicates()
  


# COMMAND ----------

# MAGIC %md
# MAGIC # 3. Covid positive

# COMMAND ----------

# sgss
# note: all records are included as every record is a "positive test"
# -- TODO: wranglers please clarify whether LAB ID 840 is still the best means of identifying pillar 1 vs 2
# -- CASE WHEN REPORTING_LAB_ID = '840' THEN "pillar_2" ELSE "pillar_1" END as description,
#   .withColumn('description', f.when(f.col('Reporting_Lab_ID') == '840', 'pillar_2').otherwise('pillar_1'))\
_sgss_pos = _sgss\
  .withColumn('covid_phenotype', f.lit('01_Covid_positive_test'))\
  .withColumn('clinical_code', f.lit(''))\
  .withColumn('description', f.lit(''))\
  .withColumn('covid_status', f.lit(''))\
  .withColumn('code', f.lit(''))\
  .withColumn('source', f.lit('sgss'))\
  .select('PERSON_ID', 'DATE', 'covid_phenotype', 'clinical_code', 'description', 'covid_status', 'code', 'source')

# gdppr
# note: need to inspect and identify which are only suspected NOT confirmed!
_codelist_gdppr = codelist_covid\
  .where(f.col('name') == 'covid19')\
  .select(['code', 'term'])

_gdppr_pos = _gdppr\
  .select(['PERSON_ID', 'DATE', 'CODE'])\
  .join(f.broadcast(_codelist_gdppr), on='code', how='inner')\
  .withColumn('covid_phenotype', f.lit('01_GP_covid_diagnosis'))\
  .withColumnRenamed('CODE', 'clinical_code')\
  .withColumnRenamed('term', 'description')\
  .withColumn('covid_status', f.lit(''))\
  .withColumn('code', f.lit('SNOMED'))\
  .withColumn('source', f.lit('gdppr'))

# COMMAND ----------

# MAGIC %md # 4. Covid admission

# COMMAND ----------

# ------------------------------------------------------------------------------
# hes_apc
# ------------------------------------------------------------------------------
# any
_hes_apc_adm_any = _hes_apc\
  .where(f.col('DIAG_4_CONCAT').rlike('U07(1|2)'))\
  .withColumn('covid_phenotype', f.lit('02_Covid_admission_any_position'))\
  .withColumn('clinical_code',\
    f.when(f.col('DIAG_4_CONCAT').rlike('U071'), 'U071')\
    .when(f.col('DIAG_4_CONCAT').rlike('U072'), 'U072')\
  )\
  .withColumn('description',\
    f.when(f.col('DIAG_4_CONCAT').rlike('U071'), 'Confirmed_COVID19')\
    .when(f.col('DIAG_4_CONCAT').rlike('U072'), 'Suspected_COVID19')\
  )\
  .withColumn('covid_status',\
    f.when(f.col('DIAG_4_CONCAT').rlike('U071'), 'confirmed')\
    .when(f.col('DIAG_4_CONCAT').rlike('U072'), 'suspected')\
  )\
  .withColumn('code', f.lit('ICD10'))\
  .withColumn('source', f.lit('hes_apc'))\
  .select('PERSON_ID', 'DATE', 'covid_phenotype', 'clinical_code', 'description', 'covid_status', 'code', 'source')

# pri
_hes_apc_adm_pri = _hes_apc\
  .where(f.col('DIAG_4_01').rlike('U07(1|2)'))\
  .withColumn('covid_phenotype', f.lit('02_Covid_admission_primary_position'))\
  .withColumn('clinical_code',\
    f.when(f.col('DIAG_4_01').rlike('U071'), 'U071')\
    .when(f.col('DIAG_4_01').rlike('U072'), 'U072')\
  )\
  .withColumn('description',\
    f.when(f.col('DIAG_4_01').rlike('U071'), 'Confirmed_COVID19')\
    .when(f.col('DIAG_4_01').rlike('U072'), 'Suspected_COVID19')\
  )\
  .withColumn('covid_status',\
    f.when(f.col('DIAG_4_01').rlike('U071'), 'confirmed')\
    .when(f.col('DIAG_4_01').rlike('U072'), 'suspected')\
  )\
  .withColumn('code', f.lit('ICD10'))\
  .withColumn('source', f.lit('hes_apc'))\
  .select('PERSON_ID', 'DATE', 'covid_phenotype', 'clinical_code', 'description', 'covid_status', 'code', 'source')


# ------------------------------------------------------------------------------
# sus
# ------------------------------------------------------------------------------
# any
_sus_adm_any = _sus\
  .where(f.col('DIAG_CONCAT').rlike('U07(1|2)'))\
  .withColumn('covid_phenotype', f.lit('02_Covid_admission_any_position'))\
  .withColumn('clinical_code',\
    f.when(f.col('DIAG_CONCAT').rlike('U071'), 'U071')\
    .when(f.col('DIAG_CONCAT').rlike('U072'), 'U072')\
  )\
  .withColumn('description',\
    f.when(f.col('DIAG_CONCAT').rlike('U071'), 'Confirmed_COVID19')\
    .when(f.col('DIAG_CONCAT').rlike('U072'), 'Suspected_COVID19')\
  )\
  .withColumn('covid_status',\
    f.when(f.col('DIAG_CONCAT').rlike('U071'), 'confirmed')\
    .when(f.col('DIAG_CONCAT').rlike('U072'), 'suspected')\
  )\
  .withColumn('code', f.lit('ICD10'))\
  .withColumn('source', f.lit('sus'))\
  .select('PERSON_ID', 'DATE', 'covid_phenotype', 'clinical_code', 'description', 'covid_status', 'code', 'source')

# pri
_sus_adm_pri = _sus\
  .where(f.col('PRIMARY_DIAGNOSIS_CODE').rlike('U07(1|2)'))\
  .withColumn('covid_phenotype', f.lit('02_Covid_admission_primary_position'))\
  .withColumn('clinical_code',\
    f.when(f.col('PRIMARY_DIAGNOSIS_CODE').rlike('U071'), 'U071')\
    .when(f.col('PRIMARY_DIAGNOSIS_CODE').rlike('U072'), 'U072')\
  )\
  .withColumn('description',\
    f.when(f.col('PRIMARY_DIAGNOSIS_CODE').rlike('U071'), 'Confirmed_COVID19')\
    .when(f.col('PRIMARY_DIAGNOSIS_CODE').rlike('U072'), 'Suspected_COVID19')\
  )\
  .withColumn('covid_status',\
    f.when(f.col('PRIMARY_DIAGNOSIS_CODE').rlike('U071'), 'confirmed')\
    .when(f.col('PRIMARY_DIAGNOSIS_CODE').rlike('U072'), 'suspected')\
  )\
  .withColumn('code', f.lit('ICD10'))\
  .withColumn('source', f.lit('sus'))\
  .select('PERSON_ID', 'DATE', 'covid_phenotype', 'clinical_code', 'description', 'covid_status', 'code', 'source')


# ------------------------------------------------------------------------------
# chess
# ------------------------------------------------------------------------------
_chess_adm = _chess\
  .select(['PERSON_ID', f.col('HospitalAdmissionDate').alias('DATE')])\
  .withColumn('covid_phenotype', f.lit('02_Covid_admission_any_position'))\
  .withColumn('clinical_code', f.lit(''))\
  .withColumn('description', f.lit('HospitalAdmissionDate IS NOT null'))\
  .withColumn('covid_status', f.lit('confirmed'))\
  .withColumn('code', f.lit(''))\
  .withColumn('source', f.lit('chess'))

# COMMAND ----------

# MAGIC %md # 5. Combine

# COMMAND ----------

tmp = _sgss_pos\
  .unionByName(_gdppr_pos)\
  .unionByName(_sus_adm_any)\
  .unionByName(_sus_adm_pri)\
  .unionByName(_hes_apc_adm_any)\
  .unionByName(_hes_apc_adm_pri)\
  .unionByName(_chess_adm)

# COMMAND ----------

# MAGIC %md # 6. Check

# COMMAND ----------

# display(tmp)

# COMMAND ----------

# tmp.count()

# COMMAND ----------

# check combined
# count_var(tmp, 'PERSON_ID'); print()
# tmpt = tab(tmp, 'covid_phenotype', 'source', var2_unstyled=1); print()
# tmp1 = tmp.withColumn('source_pheno', f.concat_ws('_', f.col('source'), f.col('covid_phenotype')))
# tmpt = tabstat(tmp1, 'DATE', byvar='source_pheno', date=1); print()

# COMMAND ----------

# check individual
# count_var(_sgss_pos, 'PERSON_ID'); print()
# count_var(_gdppr_pos, 'PERSON_ID'); print()
# count_var(_sus_adm_any, 'PERSON_ID'); print()
# count_var(_sus_adm_pri, 'PERSON_ID'); print()
# count_var(_hes_apc_adm_any, 'PERSON_ID'); print()
# count_var(_hes_apc_adm_pri, 'PERSON_ID'); print()
# count_var(_chess_adm, 'PERSON_ID'); print()

# COMMAND ----------

# MAGIC %md # 7. Save

# COMMAND ----------

outName = f'{proj}_cur_covid'
tmp.write.mode('overwrite').saveAsTable(f'{dsa}.{outName}')
