###############################################################################
#
# match Nightingale metabolite column names with the full names
# created: Nov.8, 2021 (see the email messages from ZP and ES)
# Run locally
# 
# Ran once: Nightingale metabolite metadata
# install.packages("remotes")
# remotes::install_github("NightingaleHealth/ggforestplot")
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

wd="/Users/jshin/OneDrive - SickKids/neuroCHARGE_BrainMS_CircMetabolism/scripts/neuroCHARGE/UKB_Nightingale"
setwd(wd)

#------------------------------------------------------------------------------
# read in table files
ukb_metabo_vars = fread("../UKB_Nightingale/metadata/Nightingale_UKB_FieldID.txt")#this does not work
head(ukb_metabo_vars)
# ukb_metabo_vars = read_xlsx('../UKB_Nightingale/metadata/DataVariableList.xlsx',sheet=2)
ukb_metabo_vars = ukb_metabo_vars %>% 
  # dplyr::rename(Description=title) %>%
  mutate(tmpID = str_remove(tolower(Description)," percentage")) %>%
  mutate(tmpID = str_remove(tmpID," ratio ")) %>%
  mutate(tmpID = str_remove(tmpID," ratio"))
head(ukb_metabo_vars)

NG_biomarker_metadata = (ggforestplot::df_NG_biomarker_metadata)
NG_biomarker_metadata = NG_biomarker_metadata %>% dplyr::select(-unit)

NG_biomarker_metadata = NG_biomarker_metadata %>% 
  mutate(tmpID = tolower(description)) %>%
  mutate(tmpID = str_remove(tmpID,"ratio of ")) %>% 
  mutate(tmpID = str_remove(tmpID,"ratio "))
  
setdiff(ukb_metabo_vars$tmpID,NG_biomarker_metadata$tmpID)
intersect(names(ukb_metabo_vars),names(NG_biomarker_metadata))

ukb_metabo_vars2 = merge(ukb_metabo_vars,NG_biomarker_metadata)
# unlist(lapply(lapply(ukb_metabo_vars$alternative_names,str_length),which.min))
# sapply(ukb_metabo_vars$alternative_names, function(x) x[which.min(str_length(x))])

write_tsv(ukb_metabo_vars2 %>% 
            dplyr::select(-tmpID,-unit,-alternative_names,-unit) %>%
            arrange(`Field ID`),
          "metadata/NG_UKB_mathced_metaboIDs2.tsv")          
