#################################################
#  This R script that will generate a histogram plot (in .png format) 
#
#  In addition it will output a text file with numerical summaries that will need to be uploaded to the NeuroCHARGE working group 
#################################################

#---------------------------------------------------------------------------------------#
#  Run this script in the same directory where your data files are: 
#  All of the output will go into output folder 'COHORT_ANCESTRY'
#
#  Written by Jean Shin (jean.shin@sickkids.ca) and Eeva Sliz (XXX) for the CHARGE Consortium (2019)
#---------------------------------------------------------------------------------------#
op <- options(nwarnings = 10000)
# --------------------------------------
# Specify working directory where the script and data files are - EDIT!
# --------------------------------------
WorkingDirectory = './' #change if necessary

# --------------------------------------
# Set working directory
# --------------------------------------
setwd(WorkingDirectory)
input_specification_file = 'Example_cohort_specific_inputs_wi_HTN.DM_R2calc.txt'
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
                      "FactoMineR","readr","car","missMDA","reshape","cowplot","data.table",
                      'ggpubr','gridExtra',"coxme","MuMIn",
                      'data.table','relaimpo')
lapply(required.packages,install.packages.if.not.avail)
lapply(required.packages,require,character.only=T)
source('calc.relimp.lm.js.r')

# --------------------------------------
# Create output directory and copy the 'cohort_specific_inputs.txt' file into the output foder
# --------------------------------------
outfolder = paste(cohort_name,ancestry,sep="_")
outfolder = paste(outfolder,"/",sep="")
dir.create(outfolder)
file.copy(input_specification_file,outfolder,recursive = TRUE)
warning.f = paste(outfolder,'warning_model_',input_specification_file,sep="")

##
brain_columns <- c(WMH,ICV_or_BrainSize)
brain_columns_names <- unlist(str_split("WMH,ICV_or_BrainSize",","))
brain_columns_names <- brain_columns_names[!is.na(brain_columns)]
brain_columns = brain_columns[!is.na(brain_columns)]
names(brain_columns) = brain_columns_names

cov_columns <- c(FID,age,sex,years,fasting_duration,statin_use,other_lipid_lowering_med_use,BMI,current_smoking,eGFR,hypertension,diabetes)
cov_columns_names <- unlist(str_split("FID,age,sex,years,fasting_duration,statin_use,other_lipid_lowering_med_use,BMI,current_smoking,eGFR,hypertension,diabetes",","))
cov_columns_names <- c(cov_columns_names)
cov_columns_names <- cov_columns_names[!is.na(cov_columns)]
cov_columns = cov_columns[!is.na(cov_columns)]
names(cov_columns) = cov_columns_names

# Load data and create a merged file (for brain and covariate variables)
# --------------------------------------
#Check if the data and metabolite list files have been created and stored in the working directory
braindata_file=paste(WorkingDirectory,"/",BrainDataFile,sep="")
covdata_file=paste(WorkingDirectory,"/",CovDataFile,sep="")
gkmat_file=paste(WorkingDirectory,"/",kinMatrixFile,sep="")

if(!file.exists(braindata_file)) stop("brain data file does not exist: Please read the instructions.")
if(!file.exists(covdata_file)) stop("covariate data file does not exist: Please read the instructions.")

## brain data
braindata = fread(braindata_file)
braindata_colnames = names(braindata)
braindata = data.frame(braindata)
names(braindata) = braindata_colnames

## covariate data
covdata = fread(covdata_file)
covdata_colnames = names(covdata)
covdata = data.frame(covdata)
names(covdata) = covdata_colnames

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

if(!is.na(diabetes)){
  tmp = covdata[[diabetes]]
  tmp[!is.na(tmp) & tmp==code_diabetes_yes] <- "yes"
  tmp[!is.na(tmp) & tmp==code_diabetes_no] <- "no"
  covdata[[diabetes]]  <- tmp; rm(tmp)
}

## kinship matrix# object name=kmat
load(kinMatrixFile)

## merge data
braindata2 = braindata[,c(IID,brain_columns)]
covdata2 = covdata[,c(IID,cov_columns)]
data <- merge(braindata2,covdata2,by=IID,sort=F)

# --------------------------------------
#Select variables available
# --------------------------------------
outcome_vec = c("logWMH"="logWMH")
brain_vec = brain_columns[!is.na(brain_columns)]
cov_vec = cov_columns

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

####
## Define regression models to be fitted to data: 
additional.covs = c("BMI","current_smoking","eGFR","statin_use")[!is.na(c(BMI,current_smoking,eGFR,statin_use))]
add.hypertension = "hypertension"[!is.na(hypertension)]
add.diabetes = "diabetes"[!is.na(diabetes)]

#Indicate the column names of metabolites
metabo_list_file=paste(WorkingDirectory,MetaboListFiles,sep="/")
## 
for(metaboi in 1:length(platforms)){
  metabo_outfolder = paste(outfolder,platforms[metaboi],"_",biosamples[metaboi],"/",sep="")
  dir.create(metabo_outfolder)
  metaboID_ordered <- read.table(metabo_list_file[metaboi],stringsAsFactors = F)[,1]
  metaboID_ind = metaboID_ordered == metaboID
  if(!any(metaboID_ind)){
    cat(paste(metaboID,"is not available.\n"),file=file.path(metabo_outfolder,"R2_results.txt"))
  }else if(sum(metaboID_ind)>1){
    stop(paste(sum(metaboID_ordered == metaboID),metaboID,"columns exist: Check the data for duplicated columns.\n"))
  }else{#
    # --------------------------------------
    # Load data and create a merged file (for brain and covariate variables)
    # --------------------------------------
    if(typeof(cohort_specific_cov_list)=="list" & length(cohort_specific_cov_list)==length(platforms)){
      cohort_specific_covs = cohort_specific_cov_list[[metaboi]]
    }else if(typeof(cohort_specific_cov_list)=="character" | typeof(cohort_specific_cov_list)=="NULL"){
      cohort_specific_covs = cohort_specific_cov_list
    }
    
    ## 
    covs = names(cov_vec)[!names(cov_vec)%in%c('FID','age','sex','other_lipid_lowering_med_use',additional.covs,add.hypertension,add.diabetes)]
    
    #1. combined
    M1 <- paste0(c("y ~ x + age*sex + ICV_or_BrainSize",covs),collapse=" + ")
    M2 <- paste0(c(M1,additional.covs),collapse=" + ")
    
    #2. sex-stratified
    M3 <- paste0(c("y ~ x + age + ICV_or_BrainSize",covs),collapse=" + ")
    M4 <- paste0(c(M3,additional.covs),collapse=" + ")
    
    #6. hypertension- and DM-adjusted
    M2.HTN.DM <- paste0(c(M2,add.hypertension,add.diabetes),collapse =" + ")
    M4.HTN.DM <- paste0(c(M4,add.hypertension,add.diabetes),collapse=" + ")

    if(!is.null(cohort_specific_covs)){
      M1 <- paste0(c(M1,cohort_specific_covs),collapse = " + ")
      M3 <- paste0(c(M3,cohort_specific_covs),collapse = " + ")
      # HTN+DM
      M2.HTN.DM <- paste0(c(M2.HTN.DM,cohort_specific_covs),collapse = " + ")
      M4.HTN.DM <- paste0(c(M4.HTN.DM,cohort_specific_covs),collapse = " + ")
    }
    
    mod.list <- list(M1=M1,
                     M3=M3,
                     # HTN+DM
                     M2.HTN.DM=M2.HTN.DM,
                     M4.HTN.DM=M4.HTN.DM)#B
    
    if(!is.null(kmat)){#family-data
      M1_fam <- paste(M1," + (1|IID)",sep="")
      M3_fam <- paste(M3," + (1|IID)",sep="")
      M2.HTN.DM_fam <- paste(M2.HTN.DM," + (1|IID)",sep="")
      M4.HTN.DM_fam<- paste(M4.HTN.DM," + (1|IID)",sep="")

      mod.list <- c(mod.list,
                    list(M1_fam=M1_fam,
                         M3_fam=M3_fam,
                         #HTN+DM
                         M2.HTN.DM_fam=M2.HTN.DM_fam,
                         M4.HTN.DM_fam=M4.HTN.DM_fam))
    }#if(!is.null(kmat))
    print(mod.list)
    
    cat("\n-------------------------------------------------------\n",
        file=paste(outfolder,input_specification_file,sep="/"),append=T)
    capture.output(mod.list,file=warning.f,append=T)
    cat("\n-------------------------------------------------------\n",
        file=paste(outfolder,input_specification_file,sep="/"),append=T)
    
    metabodata_file=paste(WorkingDirectory,MetaboDataFiles[metaboi],sep="/")
    if(!file.exists(metabodata_file)) stop("metabolite data file does not exist: Please read the instructions.")
    
    ## read in metabolomic data
    metabodata = fread(metabodata_file)
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
    data3 <- merge(subset(braindata,!is.na(WMH),select=c(IID,brain_vec)),covdata,by=IID,sort=F)
    data3 <- merge(data3,subset(metabodata,select=c(IID,metaboID)),by=IID,sort=F)
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
    
    if(!is.null(FID)) {
      fam_tab = table(as.character(pheno_OUT[["FID"]]))
      write_tsv(as.data.frame(fam_tab),paste(metabo_outfolder,'family_freq.tsv',sep=""))
    }
    # --------------------------------------
    # Association tests via linear regression
    # --------------------------------------
    ## standardize both outcomes and metabolites in **all individuals** (cf.,'correlation')
    dta <- pheno_OUT
    dta[,c(names(outcome_vec),metaboID)] <- apply(dta[,c(names(outcome_vec),metaboID)],2,scale)
    dta[['y']] <- dta[[names(outcome_vec)]]
    dta[['x']] <- dta[[metaboID]]
    
    ## combined
    fit_M1_combined = lm(as.formula(mod.list[["M1"]]),data=dta)
    outlierID = names(outlierTest(fit_M1_combined)$rstudent[outlierTest(fit_M1_combined)$bonf.p <0.05])
    dta2 = dta
    while(length(outlierID)>0){
      dta2 = subset(dta2,!IID%in% dta2[outlierID,"IID"])
      fit_M1_combined = lm(as.formula(mod.list[["M1"]]),
                           data = dta2)
      outlierID = names(outlierTest(fit_M1_combined)$rstudent[outlierTest(fit_M1_combined)$bonf.p <0.05])
    }
    cat("nrow(dta)-nrow(dta2)=",(nrow(dta)-nrow(dta2)),"\n",file=warning.f,append=T)
    print(nrow(dta)-nrow(dta2))
    
    dta2 = na.omit(dta2)
    rownames(dta2) = dta2[["IID"]]
    fit_M1_combined = lm(as.formula(mod.list[["M1"]]),data=dta2)
    fit_M1_combined_fam = lmekin(as.formula(mod.list[["M1_fam"]]),data=dta2,
           varlist=kmat, na.action=na.omit)
    fit_M1_combined_fam0 = lmekin(as.formula(str_remove(mod.list[["M1_fam"]],"x [+]")), data=dta2,
                                  varlist=kmat, na.action=na.omit)
    #
    fit_M2_combined = lm(as.formula(mod.list[["M2.HTN.DM"]]),data=dta2)
    fit_M2_combined_fam = lmekin(as.formula(mod.list[["M2.HTN.DM_fam"]]),data=dta2,
                                 varlist=kmat, na.action=na.omit)
    fit_M2_combined_fam0 = lmekin(as.formula(str_remove(mod.list[["M2.HTN.DM_fam"]],"x [+]")), data=dta2,
                                  varlist=kmat, na.action=na.omit)
    # sex-stratified
    fit_M3_males = lm(as.formula(mod.list[["M3"]]),data=dta2,subset=sex=="M")
    fit_M3_males_fam = lmekin(as.formula(mod.list[["M3_fam"]]),
                              data=dta2,subset=sex=="M",
                                 varlist=kmat, na.action=na.omit)
    fit_M3_males_fam0 = lmekin(as.formula(str_remove(mod.list[["M3_fam"]],"x [+]")), 
                               data=dta2,subset=sex=="M",
                                  varlist=kmat, na.action=na.omit)
    
    fit_M3_females = lm(as.formula(mod.list[["M3"]]),data=dta2,subset=sex=="F")
    fit_M3_females_fam = lmekin(as.formula(mod.list[["M3_fam"]]),
                              data=dta2,subset=sex=="F",
                              varlist=kmat, na.action=na.omit)
    fit_M3_females_fam0 = lmekin(as.formula(str_remove(mod.list[["M3_fam"]],"x [+]")), 
                               data=dta2,subset=sex=="F",
                               varlist=kmat, na.action=na.omit)
    
    fit_M4_males = lm(as.formula(mod.list[["M4.HTN.DM"]]),data=dta2,subset=sex=="M")
    fit_M4_males_fam = lmekin(as.formula(mod.list[["M4.HTN.DM_fam"]]),
                              data=dta2,subset=sex=="M",
                              varlist=kmat, na.action=na.omit)
    fit_M4_males_fam0 = lmekin(as.formula(str_remove(mod.list[["M4.HTN.DM_fam"]],"x [+]")), 
                                  data=dta2,subset=sex=="M",
                                  varlist=kmat, na.action=na.omit)
    
    fit_M4_females = lm(as.formula(mod.list[["M4.HTN.DM"]]),
                        data=dta2,subset=sex=="F")
    fit_M4_females_fam = lmekin(as.formula(mod.list[["M4.HTN.DM_fam"]]),
                              data=dta2,subset=sex=="F",
                              varlist=kmat, na.action=na.omit)
    fit_M4_females_fam0 = lmekin(as.formula(str_remove(mod.list[["M4.HTN.DM_fam"]],"x [+]")), 
                                  data=dta2,subset=sex=="F",
                                  varlist=kmat, na.action=na.omit)
    
    ## writing output
    #
    cat("*--------combined_base--------*\n",file=file.path(metabo_outfolder,"R2_results.txt"))
    capture.output(summary(fit_M1_combined),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(calc.relimp(fit_M1_combined),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    r2_combined_base = calc.relimp(fit_M1_combined)$lmg
    R_AB_M1_combined_fam = r.squaredLR(fit_M1_combined_fam)
    R_A_M1_combined_fam = r.squaredLR(fit_M1_combined_fam0)
    #R_AB_M1_combined_fam = attr(R_AB,"adj.r.squared") 
    #R_A_M1_combined_fam = attr(R_A,"adj.r.squared")
    #
    cat("\n*--------combined_full--------*\n",file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(summary(fit_M2_combined),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(calc.relimp(fit_M2_combined),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    r2_combined_full = calc.relimp(fit_M2_combined)$lmg
    R_AB_M2_combined_fam = r.squaredLR(fit_M2_combined_fam)
    R_A_M2_combined_fam = r.squaredLR(fit_M2_combined_fam0)
    #
    cat("\n*--------males_base--------*\n",file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(summary(fit_M3_males),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(calc.relimp(fit_M3_males),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    r2_males_base = calc.relimp(fit_M3_males)$lmg
    R_AB_M3_males_fam = r.squaredLR(fit_M3_males_fam)
    R_A_M3_males_fam = r.squaredLR(fit_M3_males_fam0)
    
    #
    cat("\n*--------females_base--------*\n",file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(summary(fit_M3_females),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(calc.relimp(fit_M3_females),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    r2_females_base = calc.relimp(fit_M3_females)$lmg
    R_AB_M3_females_fam = r.squaredLR(fit_M3_females_fam)
    R_A_M3_females_fam = r.squaredLR(fit_M3_females_fam0)
    
    #
    cat("\n*--------males_full--------*\n",file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(summary(fit_M4_males),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(calc.relimp(fit_M4_males),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    r2_males_full = calc.relimp(fit_M4_males)$lmg
    R_AB_M4_males_fam = r.squaredLR(fit_M4_males_fam)
    R_A_M4_males_fam = r.squaredLR(fit_M4_males_fam0)
    
    #
    cat("\n*--------females_full--------*\n",file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(summary(fit_M4_females),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    capture.output(calc.relimp(fit_M4_females),file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
    r2_females_full = calc.relimp(fit_M4_females)$lmg
    R_AB_M4_females_fam = r.squaredLR(fit_M4_females_fam)
    R_A_M4_females_fam = r.squaredLR(fit_M4_females_fam0)
    
    
    # save lmg metrics for all subsets and models
    save(r2_combined_base,r2_combined_full,r2_males_base,r2_males_full,r2_females_base,r2_females_full,
         file=file.path(metabo_outfolder,'r2_lmg_res.Rdata'))
    save(R_A_M1_combined_fam,
         R_A_M2_combined_fam,
         R_A_M3_females_fam,
         R_A_M3_males_fam,
         R_A_M4_females_fam,
         R_A_M4_males_fam,
         R_AB_M1_combined_fam,
         R_AB_M2_combined_fam,
         R_AB_M3_females_fam,
         R_AB_M3_males_fam,
         R_AB_M4_females_fam,
         R_AB_M4_males_fam,
         file=file.path(metabo_outfolder,'psuedoR2_mods_wi_and_wo_x.Rdata'))
  }
  ## any warnings
  cat("\n\n",file=file.path(metabo_outfolder,"R2_results.txt"),append = TRUE)
  capture.output(summary(warnings()),file=file.path(metabo_outfolder,"R2_results.txt"),append=T)
}

capture.output(summary(warnings()),file = warning.f,append=T)
options(op)
