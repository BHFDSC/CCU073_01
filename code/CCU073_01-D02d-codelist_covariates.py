# Databricks notebook source
# MAGIC %md # CCU073_01-D02d_codelist_covariates
# MAGIC
# MAGIC **Description** This notebook creates the covariates codelist needed for CCU073_01:
# MAGIC - Smoking
# MAGIC - BMI
# MAGIC - Blood pressure (Systolic , Diastolic)
# MAGIC - Total cholesterol
# MAGIC - HDL-cholesterol
# MAGIC - eGFR
# MAGIC - Creatinine
# MAGIC - HbA1c
# MAGIC - Glucose level
# MAGIC
# MAGIC **Author(s)** Genevieve Cezard, Tom Bolton
# MAGIC
# MAGIC **Created on** 2023.11.27
# MAGIC
# MAGIC **Last updated on** 2023.12.13
# MAGIC
# MAGIC **Data input** functions - libraries - parameters
# MAGIC
# MAGIC **Data output** \
# MAGIC CCU073_01_out_codelist_covariates
# MAGIC
# MAGIC **Reviewers** Tom Bolton, Genevieve Cezard, Wen Shi
# MAGIC
# MAGIC **Reviewed** 2024.10.25
# MAGIC
# MAGIC **Acknowledgements** Based on previous work by Tom Bolton, John Nolan, and earlier projects including CCU051, CCU002_07 and CCU004_01.
# MAGIC
# MAGIC **Notes**

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

spark.sql(f'REFRESH TABLE .gdppr_cluster_refset')
gdppr_refset = spark.table(path_ref_gdppr_refset)


# COMMAND ----------

# MAGIC %md # 3. Codelists

# COMMAND ----------

# Smoking
# BMI
# Blood pressure (Systolic, Diastolic)
# Total cholesterol
# HDL-cholesterol
# eGFR
# Creatinine
# HbA1c
# Glucose level

# COMMAND ----------

# MAGIC %md ## 3.1. Smoking

# COMMAND ----------

# Codelist checked - smoking codelist in CCU004_01 is an improved version of CCU002_01 codelist but with duplicates
# The Codelist below is reordered and de-duplicated 
tmp_smoking_status = spark.createDataFrame(
  [
    #CURRENT SMOKER
    ("134406006","Smoking reduced (finding)","Current-smoker","Unknown"),
    ("160603005","Light cigarette smoker (1-9 cigs/day) (finding)","Current-smoker","Light"),
    ("160604004","Moderate cigarette smoker (10-19 cigs/day) (finding)","Current-smoker","Moderate"),
    ("160605003","Heavy cigarette smoker (20-39 cigs/day) (finding)","Current-smoker","Heavy"),
    ("160606002","Very heavy cigarette smoker (40+ cigs/day) (finding)","Current-smoker","Heavy"),
    ("160612007","Keeps trying to stop smoking (finding)","Current-smoker","Unknown"),
    ("160613002","Admitted tobacco consumption possibly untrue (finding)","Current-smoker","Unknown"),
    ("160616005","Trying to give up smoking (finding)","Current-smoker","Unknown"),
    ("160619003","Rolls own cigarettes (finding)","Current-smoker","Unknown"),
    ("203191000000107","Wants to stop smoking (finding)","Current-smoker","Unknown"),
    ("225934006","Smokes in bed (finding)","Current-smoker","Unknown"),
    ("230056004","Cigarette consumption (observable entity)","Current-smoker","Unknown"),
    ("230057008","Cigar consumption (observable entity)","Current-smoker","Unknown"),
    ("230058003","Pipe tobacco consumption (observable entity)","Current-smoker","Unknown"),
    ("230059006","Occasional cigarette smoker (finding)","Current-smoker","Light"),
    ("230060001","Light cigarette smoker (finding)","Current-smoker","Light"),
    ("230062009","Moderate cigarette smoker (finding)","Current-smoker","Moderate"),
    ("230063004","Heavy cigarette smoker (finding)","Current-smoker","Heavy"),
    ("230064005","Very heavy cigarette smoker (finding)","Current-smoker","Heavy"),
    ("230065006","Chain smoker (finding)","Current-smoker","Heavy"),
    ("266918002","Tobacco smoking consumption (observable entity)","Current-smoker","Unknown"),
    ("266920004","Trivial cigarette smoker (less than one cigarette/day) (finding)","Current-smoker","Light"),
    ("266929003","Smoking started (finding)","Current-smoker","Unknown"),
    ("308438006","Smoking restarted (finding)","Current-smoker","Unknown"),
    ("394871007","Thinking about stopping smoking (finding)","Current-smoker","Unknown"),
    ("394872000","Ready to stop smoking (finding)","Current-smoker","Unknown"),
    ("394873005","Not interested in stopping smoking (finding)","Current-smoker","Unknown"),
    ("401159003","Reason for restarting smoking (observable entity)","Current-smoker","Unknown"),
    ("401201003","Cigarette pack-years (observable entity)","Current-smoker","Unknown"),
    ("413173009","Minutes from waking to first tobacco consumption (observable entity)","Current-smoker","Unknown"),
    ("428041000124106","Occasional tobacco smoker (finding)","Current-smoker","Light"),
    ("446172000","Failed attempt to stop smoking (finding)","Current-smoker","Unknown"),
    ("449868002","Smokes tobacco daily (finding)","Current-smoker","Unknown"),
    ("56578002","Moderate smoker (20 or less per day) (finding)","Current-smoker","Moderate"),
    ("56771006","Heavy smoker (over 20 per day) (finding)","Current-smoker","Heavy"),
    ("59978006","Cigar smoker (finding)","Current-smoker","Unknown"),
    ("65568007","Cigarette smoker (finding)","Current-smoker","Unknown"),
    ("77176002","Smoker (finding)","Current-smoker","Unknown"),
    ("82302008","Pipe smoker (finding)","Current-smoker","Unknown"),
    ("836001000000109","Waterpipe tobacco consumption (observable entity)","Current-smoker","Unknown"),
    # EX SMOKER
    ("1092031000000108","Ex-smoker amount unknown (finding)","Ex-smoker","Unknown"),
    ("1092041000000104","Ex-very heavy smoker (40+/day) (finding)","Ex-smoker","Unknown"),
    ("1092071000000105","Ex-heavy smoker (20-39/day) (finding)","Ex-smoker","Unknown"),
    ("1092091000000109","Ex-moderate smoker (10-19/day) (finding)","Ex-smoker","Unknown"),
    ("1092111000000104","Ex-light smoker (1-9/day) (finding)","Ex-smoker","Unknown"),
    ("1092131000000107","Ex-trivial smoker (<1/day) (finding)","Ex-smoker","Unknown"),
    ("160617001","Stopped smoking (finding)","Ex-smoker","Unknown"),
    ("160620009","Ex-pipe smoker (finding)","Ex-smoker","Unknown"),
    ("160621008","Ex-cigar smoker (finding)","Ex-smoker","Unknown"),
    ("160625004","Date ceased smoking (observable entity)","Ex-smoker","Unknown"),
    ("228486009","Time since stopped smoking (observable entity)","Ex-smoker","Unknown"),
    ("266921000","Ex-trivial cigarette smoker (<1/day) (finding)","Ex-smoker","Unknown"),
    ("266922007","Ex-light cigarette smoker (1-9/day) (finding)","Ex-smoker","Unknown"),
    ("266923002","Ex-moderate cigarette smoker (10-19/day) (finding)","Ex-smoker","Unknown"),
    ("266924008","Ex-heavy cigarette smoker (20-39/day) (finding)","Ex-smoker","Unknown"),
    ("266925009","Ex-very heavy cigarette smoker (40+/day) (finding)","Ex-smoker","Unknown"),
    ("266928006","Ex-cigarette smoker amount unknown (finding)","Ex-smoker","Unknown"),
    ("281018007","Ex-cigarette smoker (finding)","Ex-smoker","Unknown"),
    ("360890004","Intolerant ex-smoker (finding)","Ex-smoker","Unknown"),
    ("360900008","Aggressive ex-smoker (finding)","Ex-smoker","Unknown"),    
    ("48031000119106","Ex-smoker for more than 1 year (finding)","Ex-smoker","Unknown"),
    ("492191000000103","Ex roll-up cigarette smoker (finding)","Ex-smoker","Unknown"),
    ("53896009","Tolerant ex-smoker (finding)","Ex-smoker","Unknown"),
    ("735112005","Date ceased using moist tobacco (observable entity)","Ex-smoker","Unknown"),
    ("735128000","Ex-smoker for less than 1 year (finding)","Ex-smoker","Unknown"),
    ("8517006","Ex-smoker (finding)","Ex-smoker","Unknown"),   
    # NEVER SMOKER
    ("221000119102","Never smoked any substance (finding)","Never-smoker","NA"),
    ("266919005","Never smoked tobacco (finding)","Never-smoker","NA")
  ],
  ['code', 'term', 'smoking_status', 'severity']
)
codelist_smoking_status = (
  tmp_smoking_status
  .distinct()
  .withColumn('name',
    f.when(f.col('smoking_status') == 'Current-smoker', f.lit('smoking_current'))
     .when(f.col('smoking_status') == 'Ex-smoker', f.lit('smoking_ex'))
     .when(f.col('smoking_status') == 'Never-smoker', f.lit('smoking_never'))
     .otherwise(f.col('smoking_status'))
  )
  .withColumn('terminology', f.lit('SNOMED'))
  .select('name', 'terminology', 'code', 'term')
)  

# check
tmpt = tab(codelist_smoking_status, 'name', 'terminology'); print()
print(codelist_smoking_status.orderBy('name', 'terminology', 'code').limit(10).toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 3.2. BMI

# COMMAND ----------

# Codelist checked - same as CCU051, from CCU002_01
list_bmi = [
  '722595002'
  , '914741000000103'
  , '914731000000107'
  , '914721000000105'
  , '35425004'
  , '48499001'
  , '301331008'
  , '6497000'
  , '310252000'
  , '427090001'
  , '408512008'
  , '162864005'
  , '162863004'
  , '412768003'
  , '60621009'
  , '846931000000101'
]

tmp1_bmi = (
  gdppr_refset
  .select('ConceptId', 'ConceptId_description')
  .where(f.col('ConceptId').isin(list_bmi))
  .dropDuplicates(['ConceptId'])
  .select(f.col('ConceptId').alias('code'), f.col('ConceptId_description').alias('term'))
  .orderBy('code')
)

# check
print(tmp1_bmi.toPandas().to_string()); print()

# BMIVAL_COD
tmp2_bmi = (
  gdppr_refset
  .select('Cluster_ID', 'Cluster_Desc', 'ConceptId', 'ConceptId_description')
  .where(f.col('Cluster_ID').isin(['BMIVAL_COD']))
  .dropDuplicates(['ConceptId'])
  .select(f.col('ConceptId').alias('code'))
  .orderBy('code')
)

# check
print(tmp2_bmi.toPandas().to_string()); print()

# merge
tmp3_bmi = merge(tmp1_bmi, tmp2_bmi, ['code'], validate='1:1', assert_results=['both', 'left_only'])

# check
print(tmp3_bmi.toPandas().to_string()); print()

# prepare
codelist_bmi = (
  tmp3_bmi
  .withColumn('name', f.lit('bmi'))
  .withColumn('terminology', f.lit('SNOMED'))
  .select('name', 'terminology', 'code', 'term')
)  

# check
tmpt = tab(codelist_bmi, 'name', 'terminology'); print()
#print(codelist_bmi.orderBy('name', 'terminology', 'code').toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ## 3.3. Blood pressure 

# COMMAND ----------

tmp_bp = (
  gdppr_refset
  .where(f.col('Cluster_ID') == 'BP_COD')
  .dropDuplicates(['ConceptId'])
  .select('ConceptId', 'ConceptId_description')
  .withColumn('flag_sys', f.when(f.lower(f.col('ConceptId_description')).rlike('systolic'), 1).otherwise(0))
  .withColumn('flag_dia', f.when(f.lower(f.col('ConceptId_description')).rlike('diastolic'), 1).otherwise(0))
  .withColumn('flag', 
              f.when(f.col('flag_sys') == 1, 'sys')
              .when(f.col('flag_dia') == 1, 'dia')
              .otherwise('gen')
             )
)

# check
tmpt = tab(tmp_bp, 'flag_sys', 'flag_dia'); print()
assert tmp_bp.where((f.col('flag_sys') == 1) & (f.col('flag_dia') == 1)).count() == 0
tmpt = tab(tmp_bp, 'flag'); print()
print(tmp_bp.drop('flag_sys', 'flag_dia').orderBy('ConceptId').toPandas().to_string()); print()

# COMMAND ----------

# MAGIC %md ### 3.3.1. Systolic blood pressure 

# COMMAND ----------

# filter
tmp2_sbp = (
  tmp_bp
  .where(f.col('flag').isin(['sys', 'gen']))
  .drop('flag_sys', 'flag_dia')
)

# check
tmpt = tab(tmp2_sbp, 'flag'); print()

# prepare
codelist_sbp = (
  tmp2_sbp
  .withColumn('name', f.lit('sbp'))
  .withColumn('terminology', f.lit('SNOMED'))
  .select('name', 'terminology', f.col('ConceptId').alias('code'), f.col('ConceptId_description').alias('term'), f.col('flag').alias('flag_bp'))
)  

# check
tmpt = tab(codelist_sbp, 'name'); print()
#print(codelist_sbp.orderBy('name', 'terminology', 'code').toPandas().to_string()); print()

# COMMAND ----------

# Exclude the following SNOMED codes: 
# name	terminology	code	term
# sbp	SNOMED	"775671000000104"	Post exercise systolic blood pressure response normal (finding)
# sbp	SNOMED	"707303003"	Post exercise systolic blood pressure response abnormal (finding)
# sbp	SNOMED	"707304009"	Post exercise systolic blood pressure response normal (finding)

# flag
codelist_sbp = (
  codelist_sbp
  .withColumn('flag_exclusion', 
              f.when((f.col('name') == 'sbp') & (f.col('code').isin(['775671000000104', '707303003', '707304009'])), 1)
              .otherwise(0)
             )
)

# check
tmpt = tab(codelist_sbp, 'flag_exclusion'); print()
print(codelist_sbp.where(f.col('flag_exclusion') == 1).orderBy('name', 'terminology', 'code').toPandas().to_string()); print()

# filter
codelist_sbp = (
  codelist_sbp
  .where(f.col('flag_exclusion') == 0)
)

# check
tmpt = tab(codelist_sbp, 'flag_exclusion'); print()

# tidy
codelist_sbp = (
  codelist_sbp
  .drop('flag_exclusion')
)

# COMMAND ----------

# MAGIC %md ### 3.3.2. Diastolic blood pressure 

# COMMAND ----------

# filter
tmp3_dbp = (
  tmp_bp
  .where(f.col('flag').isin(['dia', 'gen']))
  .drop('flag_sys', 'flag_dia')
)

# check
tmpt = tab(tmp3_dbp, 'flag'); print()

# prepare
codelist_dbp = (
  tmp3_dbp
  .withColumn('name', f.lit('dbp'))
  .withColumn('terminology', f.lit('SNOMED'))
  .select('name', 'terminology', f.col('ConceptId').alias('code'), f.col('ConceptId_description').alias('term'), f.col('flag').alias('flag_bp'))
)  

# check
tmpt = tab(codelist_dbp, 'name'); print()
#print(codelist_dbp.orderBy('name', 'terminology', 'code').toPandas().to_string()); print()


# COMMAND ----------

# MAGIC %md ## 3.4. Others

# COMMAND ----------

list_cod = ['CHOL2_COD', 'HDLCCHOL_COD', 'IFCCHBAM_COD', 'EGFR_COD','CRE_COD','GLUC_COD']
tmp_meas = (
  gdppr_refset
  .select('Cluster_ID', 'Cluster_Desc', 'ConceptId', 'ConceptId_description')
  .where(f.col('Cluster_ID').isin(list_cod))
  .withColumn('name',
    f.when(f.col('Cluster_ID') == 'CHOL2_COD',    f.lit('tchol'))
     .when(f.col('Cluster_ID') == 'HDLCCHOL_COD', f.lit('hdl'))
     .when(f.col('Cluster_ID') == 'EGFR_COD',     f.lit('egfr'))
     .when(f.col('Cluster_ID') == 'CRE_COD',     f.lit('Creatinine'))
     .when(f.col('Cluster_ID') == 'IFCCHBAM_COD', f.lit('hba1c'))
     .when(f.col('Cluster_ID') == 'GLUC_COD', f.lit('Glucose_level'))
     .otherwise(f.col('Cluster_ID'))
  )
  .withColumn('terminology', f.lit('SNOMED'))
)  
  
# check  
tmpt = tab(tmp_meas, 'name', 'Cluster_ID'); print()
with pd.option_context('display.max_colwidth', None): 
  tmpt = tab(tmp_meas, 'Cluster_ID', 'Cluster_Desc', var2_wide=0); print()
tmpt = tab(tmp_meas, 'name', 'terminology'); print()
print(tmp_meas.drop('Cluster_ID', 'Cluster_Desc').orderBy('name', 'terminology', 'ConceptId').toPandas().to_string()); print()  
    
# prepare
codelist_meas = (
  tmp_meas
  .select('name', 'terminology', f.col('ConceptId').alias('code'), f.col('ConceptId_description').alias('term'))
)
 
# check
tmpt = tab(codelist_meas, 'name', 'terminology'); print()
print(codelist_meas.orderBy('name', 'terminology', 'code').toPandas().to_string()); print()  

# COMMAND ----------

# MAGIC %md # 4. Combine

# COMMAND ----------

# append (union) codelists defined above
# harmonise columns before appending
codelist_all = []
for indx, clist in enumerate([clist for clist in globals().keys() if (bool(re.match('^codelist_.*', clist))) & (clist not in ['codelist_match', 'codelist_match_summ', 'codelist_match_stages_to_run', 'codelist_match_v2_test', 'codelist_tmp', 'codelist_all'])]):
  print(f'{0 if indx<10 else ""}' + str(indx) + ' ' + clist)
  codelist_tmp = globals()[clist]
  if(indx == 0):
    codelist_all = codelist_tmp
  else:
    # pre unionByName
    for col in [col for col in codelist_tmp.columns if col not in codelist_all.columns]:
      # print('  M - adding column: ' + col)
      codelist_all = codelist_all.withColumn(col, f.lit(None))
    for col in [col for col in codelist_all.columns if col not in codelist_tmp.columns]:
      # print('  C - adding column: ' + col)
      codelist_tmp = codelist_tmp.withColumn(col, f.lit(None))
      
    # unionByName  
    codelist_all = (
      codelist_all
      .unionByName(codelist_tmp)
    )
  
# order  
codelist = (
  codelist_all
  .orderBy('name', 'terminology', 'code')
)

# COMMAND ----------

# MAGIC %md # 5. Check

# COMMAND ----------

# check 
tmpt = tab(codelist, 'name', 'terminology'); print()

# check 
tmpt = tab(codelist, 'name', 'flag_bp'); print()

# COMMAND ----------

# check
# display(codelist)

# COMMAND ----------

# MAGIC %md # 6. Save

# COMMAND ----------

outName = f'{proj}_out_codelist_covariates'
codelist.write.mode('overwrite').saveAsTable(f'{dsa}.{outName}')
