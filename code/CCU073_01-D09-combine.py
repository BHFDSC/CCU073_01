# Databricks notebook source
# MAGIC %md # CCU073_01-D09-combine
# MAGIC
# MAGIC **Description** This notebook combine all the curated data. Each person has one row.
# MAGIC
# MAGIC **Authors** Wen Shi, Genevieve Cezard
# MAGIC
# MAGIC **Created on** 2024.01.20
# MAGIC
# MAGIC **Last updated on** 2024.11.14
# MAGIC
# MAGIC **Data input** functions - libraries - parameters\
# MAGIC ccu073_01_out_cohort_year2020\
# MAGIC ccu073_01_out_cohort_year2022\
# MAGIC ccu073_01_out_outcome_year2020\
# MAGIC ccu073_01_out_outcome_year2022\
# MAGIC ccu073_01_out_covariates_year2020\
# MAGIC ccu073_01_out_covariates_year2022\
# MAGIC ccu073_01_cur_covid\
# MAGIC ccu073_01_out_vaccination
# MAGIC
# MAGIC **Data outputs**
# MAGIC - **`ccu073_01_combine_year2020`**
# MAGIC - **`ccu073_01_combine_year2022`**
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard, Wen Shi
# MAGIC
# MAGIC **Reviewed on** 2024.05.01 / 2024.11.17
# MAGIC
# MAGIC **Acknowledgements** 
# MAGIC
# MAGIC **Notes**

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

# MAGIC %run "./CCU073_01-D01-parameters"

# COMMAND ----------

# MAGIC %md # 1. Combine & add chunk id

# COMMAND ----------

for i in cohort:
    spark.sql(f"""REFRESH TABLE {path_out_inc_exc_cohort}{i}""")
    spark.sql(f"""REFRESH TABLE {path_out_outcomes}{i}""")
    spark.sql(f"""REFRESH TABLE {path_out_covariates}{i}""")

    _cohort = spark.table(f'{path_out_inc_exc_cohort}{i}')
    _outcome = spark.table(f'{path_out_outcomes}{i}')
    _covariates = spark.table(f'{path_out_covariates}{i}')
    _combine = _cohort.join(_outcome,'PERSON_ID').join(_covariates,'PERSON_ID')
    _combine = _combine.withColumn('CHUNK',f.floor(f.rand(1234)*10)+f.lit(1))
    assert _combine.select('PERSON_ID').count()==_cohort.select('PERSON_ID').count()

    _combine.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{proj}_combine_year{i}')


# COMMAND ----------

# MAGIC %md # 2. For cohort 2022 only

# COMMAND ----------

spark.sql(f"""REFRESH TABLE {dsa}.{proj}_combine_year2022""")
cohort2022 = spark.table(f'{dsa}.{proj}_combine_year2022')

spark.sql(f"""REFRESH TABLE {path_cur_covid}""")
covid = spark.table(f'{path_cur_covid}')

spark.sql(f"""REFRESH TABLE {path_out_vaccination}""")
vaccination = spark.table(f'{path_out_vaccination}')

# Add Covid Flag
tmp_covid = covid.where(f.col('DATE')<datetime.date(2022,1,1)).groupBy('PERSON_ID').agg(f.max('DATE').alias('DATE'))
tmp = (cohort2022
       .join(tmp_covid.withColumnRenamed('DATE','covid_date'),'PERSON_ID','left')
       .withColumn('covid_flag',f.when(f.col('covid_date').isNotNull(),f.lit(1)).otherwise(f.lit(0)))
      )   

# Add Vaccination status
tmp_final = (tmp
             .join(vaccination.select('PERSON_ID','vaccination_status'),'PERSON_ID','left')
            )

tmp_final.write.mode('overwrite').option('overwriteSchema', 'true').saveAsTable(f'{dsa}.{proj}_combine_year2022')

# COMMAND ----------

display(tmp_final)
