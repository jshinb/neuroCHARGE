############## PLINK INPUT FILES ###############################################
#
# Run this script locally
# This script will create PLINK files for GWASs of WMH and [metabolite].
#
# This script is available on the GitHub repository:
# https://github.com/jshinb/neuroCHARGE/tree/main/UKB_Nightingale/R
#
################################################################################

#------------------------------------------------------------------------------
# load packages and set working directory
options(stringsAsFactors = F)
x = c('stringr','tidyverse','dplyr','readxl',
      'here','data.table','psych',#table
      'GenABEL','mgcv',#rntransform
      'tableone','arsenal',#table
      'ggplot2','patchwork','pheatmap',#visualization
      'ggforestplot')#˚NightingaleData

lapply(x,require,character.only=T);rm(x)

wd="/Users/jshin/OneDrive - SickKids/neuroCHARGE_BrainMS_CircMetabolism/scripts/neuroCHARGE/UKB_Nightingale"
setwd(wd)

#------------------------------------------------------------------------------
# read in table files
covdata = fread("data_local/gwas_ukb_covdata.tsv")
brain = fread("data_local/gwas_brain_data.tsv")
metabo = fread("data_local/gwas_metabo_data.tsv")