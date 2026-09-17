# Databricks notebook source
# MAGIC %md # CCU073_01-D02a-codelist_quality_assurance
# MAGIC
# MAGIC **Description** This notebook creates the codelist for the quality assurance, which includes codelists for pregnancy and prostate cancer.
# MAGIC
# MAGIC **Authors** Wen Shi
# MAGIC
# MAGIC **Created on** 2023.10.
# MAGIC
# MAGIC **Last updated on** 2023.11.03
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard, Wen Shi
# MAGIC
# MAGIC **Reviewed** 2024.10.24
# MAGIC
# MAGIC **Data input** functions - libraries - parameters
# MAGIC
# MAGIC **Data output** \
# MAGIC CCU073_01_out_codelist_quality_assurance
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

# MAGIC %md # 1 Codelists

# COMMAND ----------

# MAGIC %md ## Prostate cancer codelist

# COMMAND ----------

# prostate_cancer (SNOMED codes only)
codelist_prostate_cancer = spark.createDataFrame(
  [
    ("prostate_cancer","SNOMED","126906006","Neoplasm of prostate","",""),
    ("prostate_cancer","SNOMED","81232004","Radical cystoprostatectomy","",""),
    ("prostate_cancer","SNOMED","176106009","Radical cystoprostatourethrectomy","",""),
    ("prostate_cancer","SNOMED","176261008","Radical prostatectomy without pelvic node excision","",""),
    ("prostate_cancer","SNOMED","176262001","Radical prostatectomy with pelvic node sampling","",""),
    ("prostate_cancer","SNOMED","176263006","Radical prostatectomy with pelvic lymphadenectomy","",""),
    ("prostate_cancer","SNOMED","369775001","Gleason Score 2-4: Well differentiated","",""),
    ("prostate_cancer","SNOMED","369777009","Gleason Score 8-10: Poorly differentiated","",""),
    ("prostate_cancer","SNOMED","385377005","Gleason grade finding for prostatic cancer (finding)","",""),
    ("prostate_cancer","SNOMED","394932008","Gleason prostate grade 5-7 (medium) (finding)","",""),
    ("prostate_cancer","SNOMED","399068003","Malignant tumor of prostate (disorder)","",""),
    ("prostate_cancer","SNOMED","428262008","History of malignant neoplasm of prostate (situation)","","")
  ],
  ['name', 'terminology', 'code', 'term', 'code_type', 'RecordDate']  
)

# COMMAND ----------

# MAGIC %md ## Pregnancy & birth codelist

# COMMAND ----------

# pregnancy

# v1
tmp_pregnancy_v1 = spark.createDataFrame(
  [
    ("171057006","Pregnancy alcohol education (procedure)"),
    ("72301000119103","Asthma in pregnancy (disorder)"),
    ("10742121000119104","Asthma in mother complicating childbirth (disorder)"),
    ("10745291000119103","Malignant neoplastic disease in mother complicating childbirth (disorder)"),
    ("10749871000119100","Malignant neoplastic disease in pregnancy (disorder)"),
    ("20753005","Hypertensive heart disease complicating AND/OR reason for care during pregnancy (disorder)"),
    ("237227006","Congenital heart disease in pregnancy (disorder)"),
    ("169501005","Pregnant, diaphragm failure (finding)"),
    ("169560008","Pregnant - urine test confirms (finding)"),
    ("169561007","Pregnant - blood test confirms (finding)"),
    ("169562000","Pregnant - vaginal examination confirms (finding)"),
    ("169565003","Pregnant - planned (finding)"),
    ("169566002","Pregnant - unplanned - wanted (finding)"),
    ("413567003","Aplastic anemia associated with pregnancy (disorder)"),
    ("91948008","Asymptomatic human immunodeficiency virus infection in pregnancy (disorder)"),
    ("169488004","Contraceptive intrauterine device failure - pregnant (finding)"),
    ("169508004","Pregnant, sheath failure (finding)"),
    ("169564004","Pregnant - on abdominal palpation (finding)"),
    ("77386006","Pregnant (finding)"),
    ("10746341000119109","Acquired immune deficiency syndrome complicating childbirth (disorder)"),
    ("10759351000119103","Sickle cell anemia in mother complicating childbirth (disorder)"),
    ("10757401000119104","Pre-existing hypertensive heart and chronic kidney disease in mother complicating childbirth (disorder)"),
    ("10757481000119107","Pre-existing hypertensive heart and chronic kidney disease in mother complicating pregnancy (disorder)"),
    ("10757441000119102","Pre-existing hypertensive heart disease in mother complicating childbirth (disorder)"),
    ("10759031000119106","Pre-existing hypertensive heart disease in mother complicating pregnancy (disorder)"),
    ("1474004","Hypertensive heart AND renal disease complicating AND/OR reason for care during childbirth (disorder)"),
    ("199006004","Pre-existing hypertensive heart disease complicating pregnancy, childbirth and the puerperium (disorder)"),
    ("199007008","Pre-existing hypertensive heart and renal disease complicating pregnancy, childbirth and the puerperium (disorder)"),
    ("22966008","Hypertensive heart AND renal disease complicating AND/OR reason for care during pregnancy (disorder)"),
    ("59733002","Hypertensive heart disease complicating AND/OR reason for care during childbirth (disorder)"),
    ("171054004","Pregnancy diet education (procedure)"),
    ("106281000119103","Pre-existing diabetes mellitus in mother complicating childbirth (disorder)"),
    ("10754881000119104","Diabetes mellitus in mother complicating childbirth (disorder)"),
    ("199225007","Diabetes mellitus during pregnancy - baby delivered (disorder)"),
    ("237627000","Pregnancy and type 2 diabetes mellitus (disorder)"),
    ("609563008","Pre-existing diabetes mellitus in pregnancy (disorder)"),
    ("609566000","Pregnancy and type 1 diabetes mellitus (disorder)"),
    ("609567009","Pre-existing type 2 diabetes mellitus in pregnancy (disorder)"),
    ("199223000","Diabetes mellitus during pregnancy, childbirth and the puerperium (disorder)"),
    ("199227004","Diabetes mellitus during pregnancy - baby not yet delivered (disorder)"),
    ("609564002","Pre-existing type 1 diabetes mellitus in pregnancy (disorder)"),
    ("76751001","Diabetes mellitus in mother complicating pregnancy, childbirth AND/OR puerperium (disorder)"),
    ("526961000000105","Pregnancy advice for patients with epilepsy (procedure)"),
    ("527041000000108","Pregnancy advice for patients with epilepsy not indicated (situation)"),
    ("527131000000100","Pregnancy advice for patients with epilepsy declined (situation)"),
    ("10753491000119101","Gestational diabetes mellitus in childbirth (disorder)"),
    ("40801000119106","Gestational diabetes mellitus complicating pregnancy (disorder)"),
    ("10562009","Malignant hypertension complicating AND/OR reason for care during childbirth (disorder)"),
    ("198944004","Benign essential hypertension complicating pregnancy, childbirth and the puerperium - delivered (disorder)"),
    ("198945003","Benign essential hypertension complicating pregnancy, childbirth and the puerperium - delivered with postnatal complication (disorder)"),
    ("198946002","Benign essential hypertension complicating pregnancy, childbirth and the puerperium - not delivered (disorder)"),
    ("198949009","Renal hypertension complicating pregnancy, childbirth and the puerperium (disorder)"),
    ("198951008","Renal hypertension complicating pregnancy, childbirth and the puerperium - delivered (disorder)"),
    ("198954000","Renal hypertension complicating pregnancy, childbirth and the puerperium with postnatal complication (disorder)"),
    ("199005000","Pre-existing hypertension complicating pregnancy, childbirth and puerperium (disorder)"),
    ("23717007","Benign essential hypertension complicating AND/OR reason for care during pregnancy (disorder)"),
    ("26078007","Hypertension secondary to renal disease complicating AND/OR reason for care during childbirth (disorder)"),
    ("29259002","Malignant hypertension complicating AND/OR reason for care during pregnancy (disorder)"),
    ("65402008","Pre-existing hypertension complicating AND/OR reason for care during pregnancy (disorder)"),
    ("8218002","Chronic hypertension complicating AND/OR reason for care during childbirth (disorder)"),
    ("10752641000119102","Eclampsia with pre-existing hypertension in childbirth (disorder)"),
    ("118781000119108","Pre-existing hypertensive chronic kidney disease in mother complicating pregnancy (disorder)"),
    ("18416000","Essential hypertension complicating AND/OR reason for care during childbirth (disorder)"),
    ("198942000","Benign essential hypertension complicating pregnancy, childbirth and the puerperium (disorder)"),
    ("198947006","Benign essential hypertension complicating pregnancy, childbirth and the puerperium with postnatal complication (disorder)"),
    ("198952001","Renal hypertension complicating pregnancy, childbirth and the puerperium - delivered with postnatal complication (disorder)"),
    ("198953006","Renal hypertension complicating pregnancy, childbirth and the puerperium - not delivered (disorder)"),
    ("199008003","Pre-existing secondary hypertension complicating pregnancy, childbirth and puerperium (disorder)"),
    ("34694006","Pre-existing hypertension complicating AND/OR reason for care during childbirth (disorder)"),
    ("37618003","Chronic hypertension complicating AND/OR reason for care during pregnancy (disorder)"),
    ("48552006","Hypertension secondary to renal disease complicating AND/OR reason for care during pregnancy (disorder)"),
    ("71874008","Benign essential hypertension complicating AND/OR reason for care during childbirth (disorder)"),
    ("78808002","Essential hypertension complicating AND/OR reason for care during pregnancy (disorder)"),
    ("91923005","Acquired immunodeficiency syndrome virus infection associated with pregnancy (disorder)"),
    ("10755671000119100","Human immunodeficiency virus in mother complicating childbirth (disorder)"),
    ("721166000","Human immunodeficiency virus complicating pregnancy childbirth and the puerperium (disorder)"),
    ("449369001","Stopped smoking before pregnancy (finding)"),
    ("449345000","Smoked before confirmation of pregnancy (finding)"),
    ("449368009","Stopped smoking during pregnancy (finding)"),
    ("88144003","Removal of ectopic interstitial uterine pregnancy requiring total hysterectomy (procedure)"),
    ("240154002","Idiopathic osteoporosis in pregnancy (disorder)"),
    ("956951000000104","Pertussis vaccination in pregnancy (procedure)"),
    ("866641000000105","Pertussis vaccination in pregnancy declined (situation)"),
    ("956971000000108","Pertussis vaccination in pregnancy given by other healthcare provider (finding)"),
    ("169563005","Pregnant - on history (finding)"),
    ("10231000132102","In-vitro fertilization pregnancy (finding)"),
    ("134781000119106","High risk pregnancy due to recurrent miscarriage (finding)"),
    ("16356006","Multiple pregnancy (disorder)"),
    ("237239003","Low risk pregnancy (finding)"),
    ("276367008","Wanted pregnancy (finding)"),
    ("314204000","Early stage of pregnancy (finding)"),
    ("439311009","Intends to continue pregnancy (finding)"),
    ("713575004","Dizygotic twin pregnancy (disorder)"),
    ("80997009","Quintuplet pregnancy (disorder)"),
    ("1109951000000101","Pregnancy insufficiently advanced for reliable antenatal screening (finding)"),
    ("1109971000000105","Pregnancy too advanced for reliable antenatal screening (finding)"),
    ("237238006","Pregnancy with uncertain dates (finding)"),
    ("444661007","High risk pregnancy due to history of preterm labor (finding)"),
    ("459166009","Dichorionic diamniotic twin pregnancy (disorder)"),
    ("459167000","Monochorionic twin pregnancy (disorder)"),
    ("459168005","Monochorionic diamniotic twin pregnancy (disorder)"),
    ("459171002","Monochorionic monoamniotic twin pregnancy (disorder)"),
    ("47200007","High risk pregnancy (finding)"),
    ("60810003","Quadruplet pregnancy (disorder)"),
    ("64254006","Triplet pregnancy (disorder)"),
    ("65147003","Twin pregnancy (disorder)"),
    ("713576003","Monozygotic twin pregnancy (disorder)"),
    ("171055003","Pregnancy smoking education (procedure)"),
    ("10809101000119109","Hypothyroidism in childbirth (disorder)"),
    ("428165003","Hypothyroidism in pregnancy (disorder)")
  ],
  ['code', 'term']  
)
tmp_pregnancy = tmp_pregnancy_v1\
  .withColumn('name', f.lit('pregnancy'))\
  .withColumn('terminology', f.lit('SNOMED'))\
  .select(['name', 'terminology', 'code', 'term'])

# check
count_var(tmp_pregnancy, 'code')

# save
codelist_pregnancy = tmp_pregnancy

# COMMAND ----------

# append (union) codelists defined above
codelist = codelist_prostate_cancer\
  .select('name', 'terminology', 'code', 'term')\
  .unionByName(codelist_pregnancy)\
  .orderBy('name', 'terminology', 'code')

# COMMAND ----------

# MAGIC %md # 2 Check

# COMMAND ----------

# check 
# tmpt = tab(codelist, 'name', 'terminology', var2_unstyled=1)

# COMMAND ----------

# check
# display(codelist)

# COMMAND ----------

# MAGIC %md # 3 Save

# COMMAND ----------

# save name
outName = f'{proj}_out_codelist_quality_assurance'

# save
codelist.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
