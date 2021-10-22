#*****************************************************************************#
#
#  Written by Jean Shin (jean.shin@sickkids.ca) and Eeva Sliz (eeva.sliz@oulu.fi) 
#  for the CHARGE Consortium (2019).

#  This R script that will generate a text file with 
#  numerical summaries that will need to be uploaded 
#  to the NeuroCHARGE working group.
# 
#*****************************************************************************#

#-----------------------------------------------------------------------------#
#  Run this script in the same directory where your data files are: 
#  All of the output will go into output folder 'COHORT_R2calc_ANCESTRY'
#-----------------------------------------------------------------------------#

op <- options(nwarnings = 10000)
# --------------------------------------
# Specify working directory where the script and data files are - EDIT!
# --------------------------------------
WorkingDirectory = "./"

# --------------------------------------
# Set working directory
# --------------------------------------
setwd(WorkingDirectory)
input_specification_file = 'Example_cohort_specific_inputs_2021-10-21.txt'
source(input_specification_file)

# --------------------------------------
# Install and/or load libraries
# --------------------------------------
options(stringsAsFactors = F)
options(scipen=999)
install.packages.if.not.avail <- function(pkg){
  if (!require(pkg,character.only = TRUE))
  {install.packages(pkg,dep=TRUE)}
}

options(repos="https://cran.rstudio.com" )
required.packages = c("corrplot","stringr","ggplot2","MASS","dplyr","e1071","tidyr","psych",
                      "FactoMineR","readr","car","missMDA","reshape","cowplot",'ggpubr','gridExtra',
                      'data.table','relaimpo')
lapply(required.packages,install.packages.if.not.avail)
lapply(required.packages,require,character.only=T)
source('calc.relimp.lm.js.r')

## function to writing output --------------
createOut = function(fit,ofile){
  capture.output(summary(fit),file=ofile,append = TRUE)
  capture.output(calc.relimp(fit),file=ofile,append = TRUE)
  r2_res = calc.relimp(fit)$lmg
  r2_res
}

# ------------------------------------------
# Create output directory and copy the 'cohort_specific_inputs.txt' file into the output foder
# ------------------------------------------
outfolder = paste(cohort_name,ancestry,sep="_")
outfolder = paste(outfolder,"/",sep="")
dir.create(outfolder)
file.copy(input_specification_file,outfolder)

##
brain_columns <- c(WMH,ICV_or_BrainSize)
names(brain_columns) <- c("WMH","ICV_or_BrainSize")

cov_columns <- c(age,sex,years,fasting_duration,statin_use,other_lipid_lowering_med_use,BMI,current_smoking,eGFR,hypertension,diabetes)
cov_columns_names <- unlist(str_split("age,sex,years,fasting_duration,statin_use,other_lipid_lowering_med_use,BMI,current_smoking,eGFR,hypertension,diabetes",","))
cov_columns_names <- c(cov_columns_names)
cov_columns_names <- cov_columns_names[!is.na(cov_columns)]
cov_columns = cov_columns[!is.na(cov_columns)]
names(cov_columns) = cov_columns_names

# --------------------------------------
# Load data and create a merged file (for brain and covariate variables)
# --------------------------------------
#Check if the data and metabolite list files have been created and stored in the working directory
braindata_file=paste(WorkingDirectory,"/",BrainDataFile,sep="")
covdata_file=paste(WorkingDirectory,"/",CovDataFile,sep="")
if(!file.exists(braindata_file)) stop("brain data file does not exist: Please read the instructions.")
if(!file.exists(covdata_file)) stop("covariate data file does not exist: Please read the instructions.")

## brain data
braindata = read_tsv(braindata_file)
braindata_colnames = names(braindata)
braindata = data.frame(braindata)
names(braindata) = braindata_colnames

## covariate data
covdata = read_tsv(covdata_file)
covdata_colnames = names(covdata)
covdata = data.frame(covdata)
names(covdata) = covdata_colnames
covdata = covdata %>% mutate(Batch=rbinom(nrow(covdata),size=1,prob=0.5),mri_scanner=rbinom(nrow(covdata),size=2,prob=0.5))

# recoding categorical variables
tmp = covdata[[sex]]
tmp[tmp==code_male] <- "M"
tmp[tmp==code_female] <- "F"
covdata[[sex]] <- tmp; rm(tmp)

tmp = covdata[[statin_use]]
tmp[tmp==code_on_statin] <- "yes"
tmp[tmp==code_not_on_statin] <- "no"
covdata[[statin_use]] <- tmp; rm(tmp)

if(!is.na(other_lipid_lowering_med_use)){
  tmp = covdata[[other_lipid_lowering_med_use]]
  tmp[tmp==code_on_other_lipid_lowering_med] <- "yes"
  tmp[tmp==code_not_on_other_lipid_lowering_med] <- "no"
  covdata[[other_lipid_lowering_med_use]]  <- tmp; rm(tmp)
}

if(!is.na(current_smoking)){
  tmp = covdata[[current_smoking]]
  tmp[!is.na(tmp) & tmp==code_current_smoking_yes] <- "yes"
  tmp[!is.na(tmp) & tmp==code_current_smoking_no] <- "no"
  covdata[[current_smoking]]  <- tmp; rm(tmp)
}

if(!is.na(hypertension)){
  tmp = covdata[[hypertension]]
  tmp[!is.na(tmp) & tmp==code_hypertensive_yes] <- "yes"
  tmp[!is.na(tmp) & tmp==code_hypertensive_no] <- "no"
  covdata[[hypertension]]  <- tmp; rm(tmp)
}

if( !is.na(diabetes) ){
  tmp = covdata[[diabetes]]
  tmp[!is.na(tmp) & tmp==code_diabetes_yes] <- "yes"
  tmp[!is.na(tmp) & tmp==code_diabetes_no] <- "no"
  covdata[[diabetes]]  <- tmp; rm(tmp)
}

## merge data
data <- merge(braindata,covdata,by=IID,sort=F,all.y=T)

# --------------------------------------
#Select variables available 
# --------------------------------------
outcome_vec = c("logWMH"="logWMH")
brain_vec = brain_columns[!is.na(brain_columns)]
cov_vec = cov_columns[!is.na(cov_columns)]

# transformation WMH
data[["logWMH"]] <- log(1+data[[WMH]])## following 2011 GWAS of WMH paper

# remove individuals with no outcome variables
data <- data[!is.na(data[[outcome_vec]]),]

select_vars = c(IID,outcome_vec,brain_vec,cov_vec)
if(any(!select_vars %in% names(data))){
  stop("Check the variable names and come back\n")
}

pheno_OUT<-subset(data, select = select_vars)
colnames(pheno_OUT)<-c("IID",
                       outcome_vec,
                       names(brain_vec),
                       names(cov_vec))

#Indicate the column names of metabolites
metabo_list_file=paste(WorkingDirectory,MetaboListFiles,sep="/")

## 
# calculate R2 ----------------------------------------------------------------
for(metaboi in 1:length(platforms)){
  metabo_outfolder = paste(outfolder,platforms[metaboi],"_",biosamples[metaboi],"/",sep="")
  dir.create(metabo_outfolder)
  metaboID_ordered <- read.table(metabo_list_file[metaboi],stringsAsFactors = F)[,1]
  metaboID_ind = metaboID_ordered == metaboID
  outfile=file.path(metabo_outfolder,"R2_results.txt")
  
  if(!any(metaboID_ind)){
    cat(paste(metaboID,"is not available.\n"),file=outfile)
  }else if(sum(metaboID_ind)>1){
    stop(paste(sum(metaboID_ordered == metaboID),metaboID,"columns exist: Check the data for duplicated columns.\n"))
  }else{#
    ## Define regression models to be fitted to data: 
    additional.covs = c("BMI","current_smoking","eGFR","statin_use")[!is.na(c(BMI,current_smoking,eGFR,statin_use))]
    add.hypertension = "hypertension"[!is.na(hypertension)]
    add.diabetes = "diabetes"[!is.na(diabetes)]
    if(typeof(cohort_specific_cov_list)=="list" & length(cohort_specific_cov_list)==length(platforms)){
      cohort_specific_covs = cohort_specific_cov_list[[metaboi]]
    }else if(typeof(cohort_specific_cov_list)=="character" | typeof(cohort_specific_cov_list)=="NULL"){
      cohort_specific_covs = cohort_specific_cov_list
    }
    covs = names(cov_vec)[!names(cov_vec)%in%c('age','sex','other_lipid_lowering_med_use',additional.covs,add.hypertension,add.diabetes)]
    
    # define model ----------------------------------------------------------------
    #1. combined
    M1 <- paste0(c("y ~ x + age*sex + ICV_or_BrainSize",covs),collapse=" + ")
    M2 <- paste0(c(M1,additional.covs),collapse=" + ")
    
    #2. sex-stratified
    M3 <- paste0(c("y ~ x + age + ICV_or_BrainSize",covs),collapse=" + ")
    M4 <- paste0(c(M3,additional.covs),collapse=" + ")

    #4. hypertension-adjusted
    M2.HTN <- paste0(c(M2,add.hypertension),collapse =" + ")
    M4.HTN <- paste0(c(M4,add.hypertension),collapse=" + ")

    #5. hypertension- and DM-adjusted
    M2.HTN.DM <- paste0(c(M2,add.hypertension,add.diabetes),collapse =" + ")
    M4.HTN.DM <- paste0(c(M4,add.hypertension,add.diabetes),collapse=" + ")

    if(!is.null(cohort_specific_covs)){
      M1 <- paste0(c(M1,cohort_specific_covs),collapse = " + ")
      M2 <- paste0(c(M2,cohort_specific_covs),collapse = " + ")
      M3 <- paste0(c(M3,cohort_specific_covs),collapse = " + ")
      M4 <- paste0(c(M4,cohort_specific_covs),collapse = " + ")
      # HTN
      M2.HTN <- paste0(c(M2.HTN,cohort_specific_covs),collapse = " + ")
      M4.HTN <- paste0(c(M4.HTN,cohort_specific_covs),collapse = " + ")

      # HTN+DM
      M2.HTN.DM <- paste0(c(M2.HTN.DM,cohort_specific_covs),collapse = " + ")
      M4.HTN.DM <- paste0(c(M4.HTN.DM,cohort_specific_covs),collapse = " + ")
    }
    
    mod.list <- list(M1=M1,
                     M2=M2,
                     M3=M3,
                     M4=M4,
                     #HTN
                     M2.HTN=M2.HTN,
                     M4.HTN=M4.HTN,
                     #HTN and DM
                     M2.HTN.DM=M2.HTN.DM,
                     M4.HTN.DM=M4.HTN.DM)#B
    print(mod.list)
    
    cat("\n-------------------------------------------------------\n",
        file=paste(outfolder,input_specification_file,sep="/"),append=T)
    capture.output(mod.list,
                   file=paste(outfolder,input_specification_file,sep="/"),append=T)
    cat("\n-------------------------------------------------------\n",
        file=paste(outfolder,input_specification_file,sep="/"),append=T)
  
  # --------------------------------------
  # Load data and create a merged file (for brain and covariate variables)
  # --------------------------------------
  metabodata_file=paste(WorkingDirectory,MetaboDataFiles[metaboi],sep="/")
  if(!file.exists(metabodata_file)) stop("metabolite data file does not exist: Please read the instructions.")
  
  ## read in metabolomic data
  metabodata = read_tsv(metabodata_file)
  metabodata_colnames = names(metabodata)
  metabodata = data.frame(metabodata)
  names(metabodata) = metabodata_colnames
  
  ## check the column names are OK!
  check.metabo_colnames = all(metaboID_ordered %in% metabodata_colnames)
  if(!check.metabo_colnames) {
    stop("\n The column names of the \'",MetaboDataFiles[metaboi], "\' do not match those in \'",
         MetaboDataFiles[metaboi],"\': Please check!",sep="")
  }
  ## merge data
  #--------------------------------  participants with brain and metabolite data --------------------------------#
  rm(list=ls(pattern = "pheno"))
  data3 <- merge(braindata[!is.na(braindata[[WMH]]),],covdata,by=IID,sort=F)
  data3 <- merge(data3,metabodata,by=IID,sort=F)
  data3[["logWMH"]] <- log(1+data3[[WMH]])
  select_vars = c(IID,outcome_vec,brain_vec,cov_vec,cohort_specific_covs,metaboID)
  if(any(!select_vars %in% names(data3))){
    stop("Check the variable names and come back\n")
  }
  
  pheno_OUT<-subset(data3, select = select_vars)
  colnames(pheno_OUT)<-c("IID",
                         names(outcome_vec),
                         names(brain_vec),
                         names(cov_vec),
                         cohort_specific_covs,
                         metaboID)
  
  # --------------------------------------
  # Association tests via linear regression
  # --------------------------------------
  ## standardize both outcomes and metabolites in **all individuals** (cf.,'correlation')
  dta <- pheno_OUT
  dta[,c(names(outcome_vec),metaboID)] <- apply(dta[,c(names(outcome_vec),metaboID)],2,scale)
  dta[['y']] <- dta[[names(outcome_vec)]]
  dta[['x']] <- dta[[metaboID]]
  
  ## combined
  ##
  fit_M1_combined = lm(as.formula(mod.list[["M1"]]),data=dta)
  outlierID = names(outlierTest(fit_M1_combined)$rstudent[outlierTest(fit_M1_combined)$bonf.p <0.05])
  dta2 = dta
  while(length(outlierID)>0){
    dta2 = subset(dta2,!IID%in% dta2[outlierID,"IID"])
    fit_M1_combined = lm(as.formula(mod.list[["M1"]]),
                         data = dta2)
    outlierID = names(outlierTest(fit_M1_combined)$rstudent[outlierTest(fit_M1_combined)$bonf.p <0.05])
  }
  
  fit_M2_combined = lm(as.formula(mod.list[["M2.HTN.DM"]]),data=dta2)
  
  ##sex_stratified
  fit_M3_males = lm(as.formula(mod.list[["M3"]]),data=dta2,subset=sex=="M")
  fit_M3_females = lm(as.formula(mod.list[["M3"]]),data=dta2,subset=sex=="F")
  fit_M4_males = lm(as.formula(mod.list[["M4.HTN.DM"]]),data=dta2,subset=sex=="M")
  fit_M4_females = lm(as.formula(mod.list[["M4.HTN.DM"]]),data=dta2,subset=sex=="F")
  
  cat("*--------combined_base--------*\n",file=outfile)
  r2_combined_base = createOut(fit_M1_combined,outfile)
  #
  cat("\n*--------combined_full--------*\n",file=outfile,append = TRUE)
  r2_combined_full = createOut(fit_M2_combined,outfile)
  #
  cat("\n*--------males_base--------*\n",file=outfile,append = TRUE)
  r2_males_base = createOut(fit_M3_males,outfile)
  #
  cat("\n*--------females_base--------*\n",file=outfile,append = TRUE)
  r2_females_base = createOut(fit_M3_females,outfile)
  #
  cat("\n*--------males_full--------*\n",file=outfile,append = TRUE)
  r2_males_full = createOut(fit_M4_males,outfile)
  #
  cat("\n*--------females_full--------*\n",file=outfile,append = TRUE)
  r2_females_full = createOut(fit_M4_females,outfile)
  
  # save lmg metrics for all subsets and models
  save(r2_combined_base,r2_combined_full,r2_males_base,r2_males_full,r2_females_base,r2_females_full,
       file=file.path(metabo_outfolder,'r2_lmg_res.Rdata'))
  }
  ## any warnings
  cat("\n\n",file=outfile,append = TRUE)
  capture.output(summary(warnings()),file=outfile,append=T)
}

capture.output(summary(warnings()),file=paste(outfolder,input_specification_file,sep=""),append=T)
options(op)
