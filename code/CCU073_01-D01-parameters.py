# Databricks notebook source
# MAGIC %md
# MAGIC
# MAGIC # CCU073_01-D01-parameters
# MAGIC
# MAGIC **Description** This notebook defines a set of parameters, which is loaded in each notebook in the data curation pipeline, so that helper functions and parameters are consistently available.
# MAGIC
# MAGIC **Authors** Wen Shi, Genevieve Cezard
# MAGIC
# MAGIC **Created on** 2023.10.30
# MAGIC
# MAGIC **Last updated on** 2024.04.22/2024.10.04
# MAGIC
# MAGIC **Reviewers** Genevieve Cezard, Wen Shi
# MAGIC
# MAGIC **Last Reviewed** 2024.11.17
# MAGIC
# MAGIC **Acknowledgements** Based on previous work by Tom Bolton (John Nolan, Elena Raffetti) for CCU018_01 and the CCU002 sub-projects.
# MAGIC
# MAGIC **Notes**

# COMMAND ----------

# MAGIC %run "/Shared/SHDS/common/functions"

# COMMAND ----------

import pyspark.sql.functions as f
import pandas as pd
import re

# COMMAND ----------

# MAGIC %md ## 1. Define parameters

# COMMAND ----------

# -----------------------------------------------------------------------------
# Project
# -----------------------------------------------------------------------------
proj = 'ccu073_01'

# -----------------------------------------------------------------------------
# Databases
# -----------------------------------------------------------------------------
db = ''
dbc = f''
dsa = f''

# -----------------------------------------------------------------------------
# cohort
# -----------------------------------------------------------------------------
cohort = [2020,2022]

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# data frame of datasets
tmp_archived_on = '2024-10-01'
data = [
    ['deaths',  dbc, f'deaths_{db}_archive',            tmp_archived_on, 'DEC_CONF_NHS_NUMBER_CLEAN_DEID', 'REG_DATE_OF_DEATH']
  , ['gdppr',   dbc, f'gdppr_{db}_archive',             tmp_archived_on, 'NHS_NUMBER_DEID',                'DATE']
  , ['hes_apc', dbc, f'hes_apc_all_years_archive',      tmp_archived_on, 'PERSON_ID_DEID',                 'EPISTART'] 
  , ['hes_op',  dbc, f'hes_op_all_years_archive',       tmp_archived_on, 'PERSON_ID_DEID',                 'APPTDATE'] 
  , ['hes_ae',  dbc, f'hes_ae_all_years_archive',       tmp_archived_on, 'PERSON_ID_DEID',                 'ARRIVALDATE'] 
  , ['vacc',    dbc, f'vaccine_status_{db}_archive',    tmp_archived_on, 'PERSON_ID_DEID',                 'DATE_AND_TIME']
  , ['sgss',    dbc, f'sgss_{db}_archive',              '2024-09-02',    'PERSON_ID_DEID',                 'Specimen_Date']  
  , ['sus',     dbc, f'sus_{db}_archive',               '2022-09-30',    'NHS_NUMBER_DEID',                'EPISODE_START_DATE']     
  , ['chess',   dbc, f'chess_{db}_archive',             '2024-09-02',    'PERSON_ID_DEID',                 'InfectionSwabDate']       
  , ['pmeds',   dbc, f'primary_care_meds_{db}_archive', tmp_archived_on, 'PERSON_ID_DEID',                 'ProcessingPeriodDate']
]
parameters_df_datasets = pd.DataFrame(data, columns = ['dataset', 'database', 'table', 'archived_on', 'idVar', 'dateVar'])
print('parameters_df_datasets:\n', parameters_df_datasets.to_string())
  
# note: the below is largely listed in order of appearance within the pipeline:  

# reference tables
path_ref_bhf_phenotypes  = '.bhf_covid_uk_phenotypes_20210127'
path_ref_gdppr_refset    = '.gdppr_cluster_refset'
path_ref_geog            = '.ons_chd_geo_listings'
path_ref_imd             = '.english_indices_of_dep_v02'


# curated tables
path_cur_hes_apc_long      = f'{dsa}.{proj}_cur_hes_apc_all_years_long'
path_cur_hes_apc_oper_long = f'{dsa}.{proj}_cur_hes_apc_all_years_archive_oper_long'
path_cur_deaths_long       = f'{dsa}.{proj}_cur_deaths_long'
path_cur_deaths_sing       = f'{dsa}.{proj}_cur_deaths_sing'
path_cur_lsoa_region       = f'{dsa}.{proj}_cur_lsoa_region_lookup'
path_cur_lsoa_imd          = f'{dsa}.{proj}_cur_lsoa_imd_lookup'
path_cur_covid             = f'{dsa}.{proj}_cur_covid'
path_cur_vacc              = f'{dsa}.{proj}_cur_vacc'
path_cur_vacc_reshape      = f'{dsa}.{proj}_cur_vacc_reshaped'

# # temporary tables
path_tmp_demo                  = f"{dsa}.hds_curated_assets__demographics_{re.sub('-','_',tmp_archived_on)}"
path_tmp_lsoa                  =f"{dsa}.hds_curated_assets__lsoa_multisource_{re.sub('-','_',tmp_archived_on)}"
path_tmp_rank_lsoa             =f'{dsa}.{proj}_tmp_rank_mode_lsoa_year'
path_tmp_quality_assurance_hx_1st_wide  = f'{dsa}.{proj}_tmp_quality_assurance_hx_1st_wide'
path_tmp_quality_assurance_hx_1st       = f'{dsa}.{proj}_tmp_quality_assurance_hx_1st'
path_tmp_quality_assurance_qax          = f'{dsa}.{proj}_tmp_quality_assurance_qax'
path_tmp_cohort                         = f'{dsa}.{proj}_tmp_cohort'

# out tables
path_out_codelist_outcome           = f'{dsa}.{proj}_out_codelist_cvd_outcomes'
path_out_codelist_covariates        = f'{dsa}.{proj}_out_codelist_covariates'
path_out_codelist_covid             = f'{dsa}.{proj}_out_codelist_covid'
path_out_codelist_quality_assurance = f'{dsa}.{proj}_out_codelist_quality_assurance'

path_out_skinny                     = f'{dsa}.{proj}_skinny'
path_out_lsoa                       = f'{dsa}.{proj}_out_lsoa_year'
path_out_quality_assurance          = f'{dsa}.{proj}_quality_assurance'
path_out_inc_exc_cohort             = f'{dsa}.{proj}_out_cohort_year'
path_out_flow_inc_exc               = f'{dsa}.{proj}_out_flow_inc_exc_year'
path_out_outcomes                   = f'{dsa}.{proj}_out_outcome_year'
path_out_covariates                 = f'{dsa}.{proj}_out_covariates_year'
path_out_vaccination                = f'{dsa}.{proj}_out_vaccination'


# COMMAND ----------

# Check latest archived dates for each dataset of interest - checked on 2024-10-04

#tmp = spark.table(f'{dbc}.vaccine_status_{db}_archive')
#tmpt = tab(tmp, 'archived_on')
# -> Vaccine status dataset latest update : 2024-10-01

#tmp = spark.table(f'{dbc}.deaths_{db}_archive')
#tmpt = tab(tmp, 'archived_on')
# -> Death dataset latest update : 2024-10-01

#tmp = spark.table(f'{dbc}.gdppr_{db}_archive')
#tmpt = tab(tmp, 'archived_on')
# -> GDPPR dataset latest update : 2024-10-01

#tmp = spark.table(f'{dbc}.hes_apc_all_years_archive')
#tmpt = tab(tmp, 'archived_on')
# -> HES APC dataset latest update : 2024-10-01

#tmp = spark.table(f'{dbc}.hes_op_all_years_archive')
#tmpt = tab(tmp, 'archived_on')
# -> HES Outpatient dataset latest update : 2024-10-01

#tmp = spark.table(f'{dbc}.hes_ae_all_years_archive')
#tmpt = tab(tmp, 'archived_on')
# -> HES A&E dataset latest update : 2024-10-01

#tmp = spark.table(f'{dbc}.sgss_{db}_archive')
#tmpt = tab(tmp, 'archived_on')
# -> SGSS dataset latest update : 2024-09-02

#tmp = spark.table(f'{dbc}.chess_{db}_archive')
#tmpt = tab(tmp, 'archived_on')
# -> CHESS dataset latest update :  2024-09-02

#tmp = spark.table(f'{dbc}.primary_care_meds_{db}_archive')
#tmpt = tab(tmp, 'archived_on')
# -> Medication latest update: 2024-10-01

# COMMAND ----------

# MAGIC %md ## 2. Function to get data on archived_on date

# COMMAND ----------

# function to extract the batch corresponding to the pre-defined archived_on date from the archive for the specified dataset
from pyspark.sql import DataFrame
def extract_batch_from_archive(_df_datasets: DataFrame, _dataset: str):
  
  # get row from df_archive_tables corresponding to the specified dataset
  _row = _df_datasets[_df_datasets['dataset'] == _dataset]
  
  # check one row only
  assert _row.shape[0] != 0, f"dataset = {_dataset} not found in _df_datasets (datasets = {_df_datasets['dataset'].tolist()})"
  assert _row.shape[0] == 1, f"dataset = {_dataset} has >1 row in _df_datasets"
  
  # create path and extract archived on
  _row = _row.iloc[0]
  _path = _row['database'] + '.' + _row['table']  
  _archived_on = _row['archived_on']  
  print(_path + ' (archived_on = ' + _archived_on + ')')
  
  # check path exists # commented out for runtime
#   _tmp_exists = spark.sql(f"SHOW TABLES FROM {_row['database']}")\
#     .where(f.col('tableName') == _row['table'])\
#     .count()
#   assert _tmp_exists == 1, f"path = {_path} not found"

  # extract batch
  _tmp = spark.table(_path)\
    .where(f.col('archived_on') == _archived_on)  
  
  # check number of records returned
  _tmp_records = _tmp.count()
  print(f'  {_tmp_records:,} records')
  assert _tmp_records > 0, f"number of records == 0"

  # return dataframe
  return _tmp

# COMMAND ----------

# MAGIC %md ## 3. Print defined parameters

# COMMAND ----------

print(f'Project:')
print("  {0:<22}".format('proj') + " = " + f'{proj}') 
print(f'')
print(f'Databases:')
print("  {0:<22}".format('db') + " = " + f'{db}') 
print("  {0:<22}".format('dbc') + " = " + f'{dbc}') 
print("  {0:<22}".format('dsa') + " = " + f'{dsa}') 
print(f'')
print(f'Paths:')
print(f'')
tmp = vars().copy()
for var in list(tmp.keys()):
  if(re.match('^path_.*$', var)):
    print("  {0:<22}".format(var) + " = " + tmp[var])    
print(f'')

