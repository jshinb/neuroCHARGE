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
WorkingDirectory = 

# --------------------------------------
# Set working directory
# --------------------------------------
setwd(WorkingDirectory)
input_specification_file = 'cohort_specific_inputs.txt'
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
                      "FactoMineR","readr","car","missMDA","reshape","cowplot")
lapply(required.packages,install.packages.if.not.avail)

library(corrplot)
library(stringr)
library(ggplot2)
library(MASS)
library(dplyr)
library(e1071)
library(tidyr)
library(psych)
library(FactoMineR)
library(readr)
library(car)
library(reshape)
library(cowplot)

# --------------------------------------
# Define functions
# --------------------------------------
source('define_functions.r')

# --------------------------------------
# Create output directory and copy the 'cohort_specific_inputs.txt' file into the output foder
# --------------------------------------
outfolder = paste(cohort_name,ancestry,sep="_")
outfolder = paste(outfolder,"/",sep="")
dir.create(outfolder)
file.copy(input_specification_file,outfolder)

##
brain_columns <- c(WMH,ICV_or_BrainSize)
names(brain_columns) <- c("WMH","ICV_or_BrainSize")
cov_columns <- c(age,sex,years,fasting_duration,statin_use,other_lipid_lowering_med_use,BMI,current_smoking,eGFR)
cov_columns_names <- unlist(str_split("age,sex,years,fasting_duration,statin_use,other_lipid_lowering_med_use,BMI,current_smoking,eGFR",","))
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

# --------------------------------------
# Summary statistics: Brain variables
# --------------------------------------
source('obtain_desc_brain_variables.r')

# create histograms for brainb variables: WMH, ICV or brain size
png(paste(outfolder,"brain_cov_variable_histograms%03d.png",sep=""),
    width=8.5,height=11, units = "in", res=300)
par(mfrow=c(3,2))
pnames = names(pheno_OUT)[-1]
tmp = data.frame(pheno_OUT[,-1],stringsAsFactors = T)
for(i in which(unlist(lapply(tmp,mode))=="numeric")) {
  distr.simple(tmp,pnames[i],pnames[i])
}
dev.off()

# Pairwise correlation among brain varialbles
## Pearson correlation
correlations <- cor(tmp[,which(unlist(lapply(tmp,mode))=="numeric")], use = "pairwise.complete.obs", method="pearson")

filename1 = "pearsoncorr_rawpheno_combined.tsv"
write.table(correlations, paste(outfolder,filename1,sep=""), col.names = T, row.names=T, sep="\t")

png(paste(outfolder,"brain_cov_variable_correlation_pearson.png",sep=""),width=5,height=5,units="in",res=300)
corrplot(correlations,order="original",type="lower",tl.pos="tp", tl.col="black", tl.cex=0.7)
corrplot(correlations,add=TRUE, type="upper", method="number",order="original", col="black", diag=FALSE,tl.pos="n", cl.pos="n", number.cex = 0.7)
dev.off()

## Spearman correlation
rnk_correlations<-cor(tmp[,which(unlist(lapply(tmp,mode))=="numeric")], use = "pairwise.complete.obs", method="spearman")

filename2 = "rnkcorr_rawpheno_combined.tsv"
write.table(rnk_correlations, paste(outfolder,filename2,sep=""), sep="\t", col.names=TRUE, row.names=TRUE)

png(paste(outfolder,"brain_cov_variable_correlation_spearman.png",sep=""),
    width=5,height=5, units="in", res=300)
corrplot(rnk_correlations,order="original",type="lower",tl.pos="tp", tl.col="black", tl.cex=0.7)
corrplot(rnk_correlations,add=TRUE, type="upper", method="number",order="original", col="black", diag=FALSE,tl.pos="n", cl.pos="n", number.cex = 0.7)
dev.off()

## Pairwise scatter plots and correlation coefficients
png(paste(outfolder,"brain_cov_variable_correlation_pairwise_scatter.png",sep=""),
    width=11,height=11,units="in",res=300)
pairs.panels(data.frame(pheno_OUT[,-1],stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()
capture.output(summary(warnings()),file=paste(outfolder,input_specification_file,sep=""),append=T)

rm(list=ls(pattern = "pheno"))
####
#Indicate the column names of metabolites
metabo_list_file=paste(WorkingDirectory,MetaboListFiles,sep="/")

## 

for(metaboi in 1:length(platforms)){
  
  ## Define regression models to be fitted to data: 
  additional.covs = c("BMI","current_smoking","eGFR","statin_use")[!is.na(c(BMI,current_smoking,eGFR,statin_use))]
  if(typeof(cohort_specific_cov_list)=="list" & length(cohort_specific_cov_list)==length(platforms)){
    cohort_specific_covs = cohort_specific_cov_list[[metaboi]]
  }else if(typeof(cohort_specific_cov_list)=="character" | typeof(cohort_specific_cov_list)=="NULL"){
    cohort_specific_covs = cohort_specific_cov_list
  }
  covs = names(cov_vec)[!names(cov_vec)%in%c('age','sex','other_lipid_lowering_med_use',additional.covs,cohort_specific_covs)]
  
  #1. combined
  M1 <- paste0(c("y ~ x + age*sex + ICV_or_BrainSize",covs),collapse=" + ")
  M2 <- paste0(c(M1,additional.covs),collapse=" + ")
  
  #2. sex-stratified
  M3 <- paste0(c("y ~ x + age + ICV_or_BrainSize",covs),collapse=" + ")
  M4 <- paste0(c(M3,additional.covs),collapse=" + ")
  
  #3. statin-straitified
  M5 <- paste0(c(M1,additional.covs[additional.covs != "statin_use"]),collapse=" + ")
  
  if(!is.null(cohort_specific_covs)){
    M1 <- paste0(c(M1,cohort_specific_covs),collapse = " + ")
    M2 <- paste0(c(M2,cohort_specific_covs),collapse = " + ")
    M3 <- paste0(c(M3,cohort_specific_covs),collapse = " + ")
    M4 <- paste0(c(M4,cohort_specific_covs),collapse = " + ")
    M5 <- paste0(c(M5,cohort_specific_covs),collapse = " + ")
  }
  
  mod.list <- list(M1=M1,
                   M2=M2,
                   M3=M3,
                   M4=M4,
                   M5=M5)#B
  
  print(mod.list)
  
  cat("\n-------------------------------------------------------\n",
      file=paste(outfolder,input_specification_file,sep="/"),append=T)
  capture.output(mod.list,
                 file=paste(outfolder,input_specification_file,sep="/"),append=T)
  cat("\n-------------------------------------------------------\n",
      file=paste(outfolder,input_specification_file,sep="/"),append=T)
  
  metabo_outfolder = paste(outfolder,platforms[metaboi],"_",biosamples[metaboi],"/",sep="")
  
  dir.create(metabo_outfolder)
  metaboID_ordered <- read.table(metabo_list_file[metaboi],stringsAsFactors = F)[,1]
  
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
  data2 <- merge(metabodata,covdata,by=IID,sort=F)
  # --------------------------------------
  #Select variables available 
  # --------------------------------------
  select_vars = c(IID,cov_vec,metaboID_ordered)
  if(any(!select_vars %in% names(data2))){
    stop("Check the variable names and come back\n")
  }
  
  pheno_OUT<-subset(data2, select = select_vars)
  colnames(pheno_OUT)<-c("IID",
                         names(cov_vec),
                         metaboID_ordered)
  # --------------------------------------
  # Summary statistics: with metabolite varialbes
  # --------------------------------------
  source('obtain_desc_metabo_variables.r')
  
  # --------------------------------------
  # Metabolite variables only
  # --------------------------------------
  metabo_OUT<-subset(pheno_OUT, select = metaboID_ordered) #Metabolites only
  metabo_OUT_males<-subset(pheno_OUT, sex=="M", select = metaboID_ordered) #Male metabolites only
  metabo_OUT_females<-subset(pheno_OUT, sex=="F", select = metaboID_ordered) #Female metabolites only
  
  # --------------------------------------
  # Metabolite varialbe - 'per-metabolite' rates of individuals with not-analysed or below-detection-level
  # --------------------------------------
  metabo_missing_count <- apply(metabo_OUT,2,calculate_missing_count)
  metabo_zero_count <- apply(metabo_OUT,2,calculate_zero_count)
  metabo_missing_count <- data.frame(metabo_ID=names(metabo_missing_count),n_missing=metabo_missing_count)
  metabo_zero_count <- data.frame(metabo_ID=names(metabo_zero_count),n_zero=metabo_zero_count)
  metabo_missing_zero_info <- merge(metabo_missing_count,metabo_zero_count,sort=F)
  metabo_missing_zero_info$n_sample <- nrow(metabo_OUT)
  metabo_missing_zero_info$N <-  metabo_missing_zero_info$n_sample-metabo_missing_zero_info$n_missing
  metabo_missing_zero_info$missingness_pct <-  metabo_missing_zero_info$n_missing/nrow(metabo_OUT)
  metabo_missing_zero_info$rate_zero <-  metabo_missing_zero_info$n_zero/nrow(metabo_OUT)
  
  filename3 = "metabo_missing_zero_freqs.tsv"
  write_tsv(metabo_missing_zero_info,paste(metabo_outfolder,filename3,sep=""))
  
  # --------------------------------------
  # Metabolite variable PCA
  # --------------------------------------
  incl.ind = metabo_missing_zero_info$missingness_pct < 0.05 #remove metabolites with missing rate >= 5%
  ncp = ncol(metabo_OUT[,incl.ind]) #top 10% PC
  
  ### combined
  PCA_metabo <- PCA(na.omit(metabo_OUT[,incl.ind]),scale.unit = TRUE, graph=F, ncp=ncp)
  PCA_metabo_var <- data.frame(PCA_metabo$eig)[,"cumulative.percentage.of.variance",drop=F]
  PCA_metabo_var$group <- "combined"
  PCA_metabo_var$x <- 1:nrow(PCA_metabo_var)
  PCA_metabo_loadings = data.frame(PCA_metabo$var$coord)
  PCA_metabo_loadings$group <- "combined"
  
  ### females
  PCA_metabo_F <- PCA(na.omit(metabo_OUT_females[,incl.ind]),scale.unit = TRUE, graph=F, ncp=ncp)
  PCA_metabo_F_var <- data.frame(PCA_metabo_F$eig)[,"cumulative.percentage.of.variance",drop=F]
  PCA_metabo_F_var$group <- "females"
  PCA_metabo_F_var$x <- 1:nrow( PCA_metabo_F_var)
  PCA_metabo_loadings_F = data.frame(PCA_metabo_F$var$coord)
  PCA_metabo_loadings_F$group = "females"
  
  ### males
  PCA_metabo_M <- PCA(na.omit(metabo_OUT_males[,incl.ind]),scale.unit = TRUE, graph=F, ncp=ncp)
  PCA_metabo_M_var <- data.frame(PCA_metabo_M$eig)[,"cumulative.percentage.of.variance",drop=F]
  PCA_metabo_M_var$group <- "males"
  PCA_metabo_M_var$x <- 1:nrow( PCA_metabo_M_var)
  PCA_metabo_loadings_M = data.frame(PCA_metabo_M$var$coord)
  PCA_metabo_loadings_M$group = "males"
  
  ### plot of cumulative proportion of variance explained by PC
  df = rbind(PCA_metabo_var,PCA_metabo_F_var,PCA_metabo_M_var)
  df$group <- factor(df$group,levels = c("combined","females","males"))
  
  p <- ggplot(data=df,aes(x=x,y=cumulative.percentage.of.variance,group=group,linetype=group,color=group))
  p <- p + geom_line() + xlab("PC dimension")
  p <- p + geom_hline(yintercept = 95, colour=grey(0.2))
  
  png(paste(metabo_outfolder,"metabolite_PC_curves.png",sep=""),width=8.5,height=4,units="in",res=300)
  print(p)
  dev.off()
  
  metabo_PCA_var_table = rbind(PCA_metabo_var,PCA_metabo_F_var,PCA_metabo_M_var)
  write_tsv(metabo_PCA_var_table, paste(metabo_outfolder,'metabo_PCA_var_table.tsv',sep=""))
  write_tsv(PCA_metabo_loadings, paste(metabo_outfolder,'PCA_metabo_loading_combined.tsv',sep=""))
  write_tsv(PCA_metabo_loadings_F, paste(metabo_outfolder,'PCA_metabo_loading_females.tsv',sep=""))
  write_tsv(PCA_metabo_loadings_M, paste(metabo_outfolder,'PCA_metabo_loading_males.tsv',sep=""))
  
  # --------------------------------------
  # Metabolite variable - histograms
  # --------------------------------------
  plot_metabo_histograms = TRUE
  if(plot_metabo_histograms){
    ## plotting histograms and save the plots (60 metabolites on each page)
    generate_histograms(metabo_OUT = metabo_OUT,
                        summarydata_Rdata = paste(metabo_outfolder,"desc_stats_combined_wi_metabo.Rdata",sep=""),
                        metabo_missing_zero_info=metabo_missing_zero_info,
                        metabo_columns=metaboID_ordered)
  }
  
  #--------------------------------  participants with brain and metabolite data --------------------------------#
  rm(list=ls(pattern = "pheno"))
  data3 <- merge(braindata[!is.na(braindata[[WMH]]),],covdata,by=IID,sort=F)
  data3 <- merge(data3,metabodata,by=IID,sort=F)
  data3[["logWMH"]] <- log(1+data3[[WMH]])
  select_vars = c(IID,outcome_vec,brain_vec,cov_vec,cohort_specific_covs,metaboID_ordered)
  if(any(!select_vars %in% names(data3))){
    stop("Check the variable names and come back\n")
  }
  
  pheno_OUT<-subset(data3, select = select_vars)
  colnames(pheno_OUT)<-c("IID",
                         names(outcome_vec),
                         names(brain_vec),
                         names(cov_vec),
                         cohort_specific_covs,
                         metaboID_ordered)
  # --------------------------------------
  # Summary statistics in individuals with brain + metabolite data
  # --------------------------------------
  source('obtain_desc_brain_metabo_variables.r')
  
  # --------------------------------------
  # PCA of metabolite variables
  # --------------------------------------
  metabo_OUT<-subset(pheno_OUT, select = metaboID_ordered) #Metabolites only
  
  # --------------------------------------
  # Per-metabolite' rates of individuals with not-analysed or below-detection-level
  # --------------------------------------
  metabo_missing_count <- apply(metabo_OUT,2,calculate_missing_count)
  metabo_zero_count <- apply(metabo_OUT,2,calculate_zero_count)
  metabo_missing_count <- data.frame(metabo_ID=names(metabo_missing_count),n_missing=metabo_missing_count)
  metabo_zero_count <- data.frame(metabo_ID=names(metabo_zero_count),n_zero=metabo_zero_count)
  metabo_missing_zero_info <- merge(metabo_missing_count,metabo_zero_count,sort=F)
  metabo_missing_zero_info$n_sample <- nrow(metabo_OUT)
  metabo_missing_zero_info$N <-  metabo_missing_zero_info$n_sample-metabo_missing_zero_info$n_missing
  metabo_missing_zero_info$missingness_pct <-  metabo_missing_zero_info$n_missing/nrow(metabo_OUT)
  metabo_missing_zero_info$rate_zero <-  metabo_missing_zero_info$n_zero/nrow(metabo_OUT)
  
  filename3 = "metabo_missing_zero_freqs_wi_brain.tsv"
  write_tsv(metabo_missing_zero_info,paste(metabo_outfolder,filename3,sep=""))
  
  # --------------------------------------
  # Association tests via linear regression
  # --------------------------------------
  ## standardize both outcomes and metabolites in **all individuals** (cf.,'correlation')
  dta <- pheno_OUT
  dta[,c(names(outcome_vec),metaboID_ordered)] <- apply(dta[,c(names(outcome_vec),metaboID_ordered)],2,scale)
  rownames(dta) <- dta[[IID]]
  resdir = paste(metabo_outfolder,"fit_results/",sep='')
  dir.create(resdir)
  
  ## combined
  M1_combined = fit_get_results("M1",group = "combined")
  M2_combined = fit_get_results("M2",group = "combined")
  
  ##sex_stratified
  M3_males = fit_get_results("M3",group = "males")
  M3_females = fit_get_results("M3",group = "females")
  M4_males = fit_get_results("M4",group = "males")
  M4_females = fit_get_results("M4",group = "females")
  
  ## statin_use_stratified
  M1_on_statin = fit_get_results("M1",group = "on_statin")
  M1_not_on_statin = fit_get_results("M1",group = "not_on_statin")
  M5_on_statin = fit_get_results("M5",group = "on_statin")
  M5_not_on_statin = fit_get_results("M5",group = "not_on_statin")
  capture.output(summary(warnings()),file=paste(outfolder,input_specification_file,sep=""),append=T)
}
capture.output(summary(warnings()),file=paste(outfolder,input_specification_file,sep=""),append=T)

options(op)
