###############################################################################
#
# match Nightingale metabolite column names with the full names
# created: Nov.8, 2021 (see the email messages from ZP and ES)
# Run locally
# 
# Ran once: Nightingale metabolite metadata
# install.packages("remotes")
# remotes::install_github("NightingaleHealth/ggforestplot")
# NG_biomarker_metadata = (ggforestplot::df_NG_biomarker_metadata)
# NG_biomarker_metadata = NG_biomarker_metadata %>% dplyr::select(-units)
# write_tsv(NG_biomarker_metadata ,"NG_biomarker_metadata.tsv")
# 
# output: 'NG_UKB_mathced_metaboIDs.tsv'
###############################################################################

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

wd="/Users/jshin/OneDrive - SickKids/ukbb_LDL.D_gwas/scripts/ukbb_gwas_LDL.D/R"
setwd(wd)

#------------------------------------------------------------------------------
# read in table files
ukb_metabo_vars = read_xlsx('../../../DataVariableList.xlsx',sheet=2)
ukb_metabo_vars = ukb_metabo_vars %>% 
  mutate(tmpID = str_remove(tolower(Description)," percentage")) %>% 
  mutate(tmpID = str_remove(tmpID," ratio ")) %>% 
  mutate(tmpID = str_remove(tmpID," ratio"))

NG_metabo_vars = fread('../NG_biomarker_metadata.tsv')
NG_metabo_vars = NG_metabo_vars %>% 
  mutate(tmpID = tolower(description)) %>%
  mutate(tmpID = str_remove(tmpID,"ratio of ")) %>% 
  mutate(tmpID = str_remove(tmpID,"ratio "))
  
setdiff(ukb_metabo_vars$tmpID,NG_metabo_vars$tmpID)
intersect(names(ukb_metabo_vars),names(NG_metabo_vars))

ukb_metabo_vars = merge(ukb_metabo_vars,NG_metabo_vars)
write_tsv(ukb_metabo_vars,"../NG_UKB_mathced_metaboIDs.tsv")          
