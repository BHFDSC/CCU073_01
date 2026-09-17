# Databricks notebook source
# MAGIC %md
# MAGIC # CCU073_01-D07-vaccination
# MAGIC  
# MAGIC **Description** This notebook creates the final vaccination variables needed for this project.
# MAGIC  
# MAGIC **Author(s)** Genevieve Cezard
# MAGIC
# MAGIC **Created on** 2024.04.22
# MAGIC
# MAGIC **Last updated on** 2024.04.22
# MAGIC
# MAGIC **Date last run** 2024.10.29
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC `ccu073_01_cur_vacc_reshaped`\
# MAGIC `ccu073_01_cohort_year2022`
# MAGIC
# MAGIC **Data Output**
# MAGIC - **`ccu073_01_out_vaccination`**
# MAGIC
# MAGIC **Reviewers** Wen Shi
# MAGIC
# MAGIC **Reviewed** 2024.10.29
# MAGIC
# MAGIC **Notes** \ created 1 extra vaccination variable and renamed date and type variables

# COMMAND ----------

spark.sql('CLEAR CACHE')

# COMMAND ----------

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

# MAGIC %run "/Shared/SHDS/common/functions"

# COMMAND ----------

# MAGIC %md # 0. Parameters

# COMMAND ----------

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

# MAGIC %md # 1. Data

# COMMAND ----------

spark.sql(f"""REFRESH TABLE {dsa}.{proj}_cur_vacc_reshaped""")
vaccination = spark.table(f'{dsa}.{proj}_cur_vacc_reshaped')

spark.sql(f"""REFRESH TABLE {path_out_inc_exc_cohort}2022""")
cohort2022 = spark.table(f'{path_out_inc_exc_cohort}2022')

# COMMAND ----------

# display(vaccination)

# COMMAND ----------

# MAGIC %md # 2. Prepare

# COMMAND ----------

# Get vaccination data for people in cohort 2022 and create vaccination status as of 2022-01-01
vacc1 = (cohort2022
         .select('PERSON_ID')
         .join(vaccination,'PERSON_ID','left')
         .withColumn('vaccination_status', 
                     f.when((f.col('date_3') < '2022-01-01'), "3+")
                    .when((f.col('date_2') < '2022-01-01'), "2")
                    .when((f.col('date_1') < '2022-01-01'), "1")
                    .otherwise('0')
                    )
         )

# COMMAND ----------

# rename columns i.e. date_1 -> date_vacc_1
dict_rename = {}
for col in vacc1.columns:
  col_rematch_date = re.match(r'(date)\_(\d+)', col)
  if(col_rematch_date):
    dict_rename[col] = col_rematch_date.group(1) + '_vacc_' + col_rematch_date.group(2)
    
vacc2 = rename_columns(vacc1, dict_rename); print()

# COMMAND ----------

# MAGIC %md # 3. Check

# COMMAND ----------

# count_var(tmp, 'PERSON_ID'); print()
# tmpt = tab(tmp, 'vaccination_status'); print()

# COMMAND ----------

# MAGIC %md # 4 Save

# COMMAND ----------

# save name
outName = f'{proj}_out_vaccination'.lower()

# save
vacc2.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')
