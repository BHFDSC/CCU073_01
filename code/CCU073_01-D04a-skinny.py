# Databricks notebook source
# MAGIC %md # CCU073_01-D04a-skinny
# MAGIC  
# MAGIC **Description** This notebook creates the skinny patient table, which includes key patient characteristics, using demographics table updated monthly by BHF data science centre.
# MAGIC  
# MAGIC **Authors** Wen Shi
# MAGIC
# MAGIC **Created on** 2023.11.07
# MAGIC
# MAGIC **Last updated on** 2024.10.25 
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC
# MAGIC
# MAGIC **Data output** \
# MAGIC CCU073_01_skinny
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard
# MAGIC
# MAGIC **Reviewed** 2024.11.17
# MAGIC
# MAGIC **Acknowledgements** 
# MAGIC
# MAGIC **Notes**

# COMMAND ----------

spark.sql('CLEAR CACHE')

# COMMAND ----------

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

# MAGIC %run "/Shared/SHDS/common/functions"

# COMMAND ----------

# MAGIC %run "/Shared/SHDS/common/skinny_20221113"

# COMMAND ----------

# MAGIC %md # 0 Parameters

# COMMAND ----------

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

# MAGIC %md # 1 Data

# COMMAND ----------

tmp = spark.table(path_tmp_demo)

# COMMAND ----------

# MAGIC %md # 2 Skinny

# COMMAND ----------

# rename column names and select columns wanted
tmp = (tmp.withColumnRenamed('person_id','PERSON_ID')
          .withColumnRenamed('sex','SEX')

          .withColumnRenamed('date_of_birth','DOB')
          .withColumnRenamed('ethnicity_5_group','ETHNIC_CAT')
          .withColumnRenamed('date_of_death','DOD')
          .select('PERSON_ID','SEX','DOB','ETHNIC_CAT','death_flag','DOD','in_gdppr','gdppr_min_date')
        )


# COMMAND ----------

# change sex values
tmp = tmp.withColumn('SEX',f.when(f.col('SEX')=='M',1).when(f.col('SEX')=='F',2).otherwise(f.col('SEX')))

# COMMAND ----------

# MAGIC %md # 3 Check

# COMMAND ----------

# count_var(tmp, 'PERSON_ID')

# COMMAND ----------

# display(tmp)

# COMMAND ----------

# MAGIC %md # 4 Save

# COMMAND ----------

# save name
outName = f'{proj}_skinny'.lower()


# save
tmp.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{outName}')

