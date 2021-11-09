############## Covariate Data Extraction #################################################
#
# Run this script in SciNet cluster
# 
# This script will create metabo data with independent European participants:
# [or should I run family data script?].
#
# Before running the script, execute the following code:
# module load gcc/9.2.0 r/4.0.3;R
#
# This script is available on the GitHub repository:
# https://github.com/jshinb/ukbb_gwas_LDL.D/tree/main/R
#
##########################################################################################

# load packages and define functions -----------------------------------------------------
library(tidyverse)
library(data.table)

# create subsets with the variables of interest
extract_variables = function(fname,fieldID,fieldName){
  if (!file.exists(fname)){stop ("file does not exist")
  }else{
    d = fread(fname)    
    d.sub = subset(d, select=c('eid',fieldID))
    names(d.sub) = c('eid',fieldName)
    d.sub
  }
}

# set working directory ------------------------------------------------------------------
wd = '/gpfs/fs1/home/t/tpaus/jshinb/ukbb_LDL.D_gwas'
setwd(wd)

#-------------------------------------------------------------------------------
# IDs to be excluded
#==============================================================================#
# e0. dropouts
#==============================================================================#
excl = fread('/gpfs/fs1/home/t/tpaus/jshinb/ukbb/exclusion_sample_lists/w43688_20210201.csv')$V1
print(p0 <- length(excl))#50

#==============================================================================#
# disease status for exclusion
#==============================================================================#
fe1=file.path(ukbb_data_dir,'ukb42388_18062020/ukb42388.csv')
var.names.e1 = c(
  "AD" = '42020-0.0',
  "dementia" = '42018-0.0',
  "stroke" = '42006-0.0')
e1 = extract_variables(fe1,fieldID=var.names.e1,fieldName=names(var.names.e1))
head(e2)
excl = unique(c(excl,e1$eid[!is.na(e1$dementia)]));print(p1 <- length(excl));print((p1-p0))#2782:2732
excl = unique(c(excl,e1$eid[!is.na(e1$stroke)]));print(p2 <- length(excl));print((p2-p1))#12707
n.excl = c(p0,p1,p2)
print(n.excl)
#[1]    50  2782 15489

#==============================================================================#
# IDs of participants that are British whites
#==============================================================================#
fe2 = file.path(ukbb_data_dir,'ukb37194/ukb37194.csv')
var.names.e2 = c(white.british = '22006-0.0')
e2 = extract_variables(fe2,var.names.e2,names(var.names.e2))
incl = subset(e2, white.british==1)$eid

#========================================================================================#
#NMR
#========================================================================================#
metabo1="23400-0.0"
metabon='23648-0.0'
metabo_cols = paste(23400:23748,"-0.0",sep='')
metabo_colnames = ""
f5 = '/gpfs/fs1/home/t/tpaus/jshinb/ukbb/ukb48959/ukb48959.csv'
var.names=c("avg.LDL.D" = '23432-0.0',
            "QC1" = '23704-0.0') #Clinical LDL Cholesterol, QC Flag
d5.sub = extract_variables(f5,var.names,names(var.names))
