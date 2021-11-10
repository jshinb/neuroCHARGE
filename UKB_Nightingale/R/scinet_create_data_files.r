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
# https://github.com/jshinb/neuroCHARGE/tree/main/UKB_Nightingale/R
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

#========================================================================================#
# data directory
#========================================================================================#
ukbb_data_dir = '/project/t/tpaus/tpaus/UKBB/datasets'
#ukb26444
#ukb30069
#ukb37194
#ukb38158
#ukb40646_20-02-2020
#ukb40646_23-03-2020
#ukb41763
#ukb42388_18062020
#ukb_01-04-2020
#gene_data

cat(dir(file.path(ukbb_data_dir,'ukb_01-04-2020')),sep='\n')
#41448
#41449
#41450

# set working directory ------------------------------------------------------------------
wd = '/gpfs/fs1/home/t/tpaus/jshinb/ukbb_LDL.D_gwas'
setwd(wd)

#-------------------------------------------------------------------------------
# IDs to be excluded
#==============================================================================#
# e0. dropouts
#==============================================================================#
excl = fread('/gpfs/fs1/home/t/tpaus/jshinb/ukbb/exclusion_sample_lists/w43688_20210201.csv')$V1

#==============================================================================#
# disease status for exclusion
#==============================================================================#
fe1=file.path(ukbb_data_dir,'ukb42388_18062020/ukb42388.csv')
var.names.e1 = c(
  "AD" = '42020-0.0',
  "dementia" = '42018-0.0',
  "stroke" = '42006-0.0')
e1 = extract_variables(fe1,fieldID=var.names.e1,fieldName=names(var.names.e1))
excl.dementia <- e1$eid[!is.na(e1$dementia)]
excl.stroke <- e1$eid[!is.na(e1$stroke)]

#==============================================================================#
# IDs of participants that are British whites
#==============================================================================#
fe2 = file.path(ukbb_data_dir,'ukb37194/ukb37194.csv')
var.names.e2 = c(white.british = '22006-0.0')
e2 = extract_variables(fe2,var.names.e2,names(var.names.e2))
incl = subset(e2, white.british==1)$eid

#-------------------------------------------------------------------------------
# brain, metabolite and covariate data
#==============================================================================#
# get brain data
#==============================================================================#
f1 = file.path(ukbb_data_dir,'ukb41763/ukb41763.csv')
varnames1 = c("WMH"="25781-2.0","ICV"="26521-2.0")
d1 = extract_variables(f1,varnames1,names(varnames1))
sum.na = apply(apply(subset(d1,select=-eid),2,is.na),1,sum)
table(sum.na)
d1 = subset(d1, sum.na==0);print(dim(d1))#37873
rm(sum.na)

#==============================================================================#
# get metabolite data
#==============================================================================#
f2 = '/gpfs/fs1/home/t/tpaus/jshinb/ukbb/ukb48959/ukb48959.csv'
metabo1=23400; metabon=23648
varnames2=paste(metabo1:metabon,"-0.0",sep='')
names(varnames2) = str_remove(varnames2,"-0.0")
d2 = extract_variables(f2,varnames2,names(varnames2))
sum.na = apply(apply(subset(d2,select=-eid),2,is.na),1,sum)
table(sum.na)
d2 = subset(d2,sum.na<249)
rm(sum.na)

#==============================================================================#
# select individuals with both brain and metabolite data; n=9129
#==============================================================================#
d1 = subset(d1, eid %in% d2$eid);print(dim(d1))#9129
d2 = subset(d2, eid %in% d1$eid);print(dim(d2))#9129

#==============================================================================#
# exclude withdrawals: n_remove=0
#==============================================================================#
d1 = subset(d1, !eid %in% excl);print(dim(d1))
d2 = subset(d2, !eid %in% excl);print(dim(d2))

#==============================================================================#
# exclude individuals with dementia or stroke (n_remove=106), n=9023
#==============================================================================#
length(excl.dementia)#2732
length(excl.stroke)#13134
length(intersect(excl.stroke,excl.dementia))#427
excl.dementia.stroke = d1$eid[d1$eid%in%unique(c(excl.dementia,excl.stroke))]#15439(all), 106

d1 = subset(d1, !eid %in% excl.dementia.stroke);print(dim(d1))
d2 = subset(d2, !eid %in% excl.dementia.stroke);print(dim(d2))

#==============================================================================#
# select white British individuals: n=7866
#==============================================================================#
d1 = subset(d1, eid %in% incl);print(dim(d1))
d2 = subset(d2, eid %in% incl);print(dim(d2))

# ------------------------------------------------------------------------------
#========================================================================================#
# exclude related individuals based on genetic relatedness: n=7805 (-61)
#========================================================================================#
f='/project/t/tpaus/tpaus/UKBB/datasets/ukb40646_23-03-2020/ukb43688_rel_s488264.dat'
d = fread(f)
d = d[order(d$ID1,d$ID2),]
d = subset(d, ID1>0)

kinship.info = d
kinship.info = subset(kinship.info,Kinship>0)
ind2 = kinship.info$ID1 %in%kinship.info$ID2
ind3 = kinship.info$ID1 > kinship.info$ID2
ID1 = kinship.info$ID1
ID2 = kinship.info$ID2
kinship.info$ID1[ind3] <- ID2[ind3]
kinship.info$ID2[ind3] <- ID1[ind3]
kinship.info = kinship.info[order(kinship.info$ID1,kinship.info$ID2),]

kinship.anal = subset(kinship.info,ID1 %in% d1$eid & kinship.info$ID2 %in% d1$eid)
print(dim(kinship.anal))#62

# 
ids = unique(c(kinship.anal$ID1,kinship.anal$ID2));print(length(ids))#120
fam.ids = list()
i=1
while(length(ids)>0){
  id=ids[1]
  ind.id = kinship.anal$ID1 == id | kinship.anal$ID2 == id
  fam.id = unique(c(id,kinship.anal$ID1[ind.id],kinship.anal$ID2[ind.id]))
  fam.id.tmp = NULL
  for(id in fam.id[-1]){
    ind.id = kinship.anal$ID1 == id | kinship.anal$ID2 == id
    fam.id.tmp = c(fam.id.tmp,unique(c(id,kinship.anal$ID1[ind.id],kinship.anal$ID2[ind.id])))
  }
  fam.ids[[i]] = unique(c(fam.id,fam.id.tmp))
  
  ids = ids[!ids %in% fam.ids[[i]]];print(length(ids))
  i = i+1
}
tmp.fam.ids = sort(unique(c(fam.id,fam.id.tmp)))
subset(kinship.anal,ID1%in%tmp.fam.ids&ID2%in% tmp.fam.ids)
mat = matrix(NA,nrow=length(tmp.fam.ids),ncol=length(tmp.fam.ids))
rownames(mat) <- colnames(mat) <- tmp.fam.ids

# construct kinship mat
for(i in 1:nrow(mat)){
  id1 = as.numeric(rownames(mat)[i]);print(id1)
  for(j in c(1:nrow(mat))[-i]){
    id2 = as.numeric(rownames(mat)[j]);print(id2)
    ind = (kinship.anal$ID1==id1 & kinship.anal$ID2==id2);print(sum(ind))
    ind = ind | (kinship.anal$ID1==id2 & kinship.anal$ID2==id1)
    print(sum(ind))#1
    if(sum(ind)==1){
      mat[i,j] <- mat[j,i] <- kinship.anal$Kinship[ind]
    }
  }
}
#
tmp.fam.ids = sort(unique(c(kinship.anal$ID1,kinship.anal$ID2)))
print(length(tmp.fam.ids))
mat = matrix(NA,nrow=length(tmp.fam.ids),ncol=length(tmp.fam.ids))
rownames(mat) <- colnames(mat) <- tmp.fam.ids

#
for(i in 1:nrow(mat)){
  id1 = as.numeric(rownames(mat)[i]);#print(id1)
  for(j in c(1:nrow(mat))[-i]){
    id2 = as.numeric(rownames(mat)[j]);#print(id2)
    ind = (kinship.anal$ID1==id1 & kinship.anal$ID2==id2);#print(sum(ind))
    ind = ind | (kinship.anal$ID1==id2 & kinship.anal$ID2==id1)
    #print(sum(ind))#1
    if(sum(ind)==1){
      mat[i,j] <- mat[j,i] <- kinship.anal$Kinship[ind]
    }
  }
  cat("*",sep='')
}
#
nrel = c()
for(i in 1:nrow(mat)){
  nrel = c(nrel,sum(!is.na(mat[,i])))
}
diag(mat) <- 1

o = order(nrel,decreasing = T)
ids = rownames(mat)[o]
fam.ids = list()
i = 1
while(length(ids)>0){
  print(i)
  id = ids[1]
  fam.idsj = rownames(mat)[!is.na(mat[,id])]
  for(idj in fam.idsj[fam.idsj!=id]){
    print(idj)
    fam.idsj = c(fam.idsj,rownames(mat)[!is.na(mat[,idj])])
  }
  fam.idsj = unique(fam.idsj)
  mat[fam.idsj,fam.idsj]
  fam.ids[[i]] = fam.idsj
  ids = ids[!ids %in% fam.idsj]
  i = i+1
}
#
dup.inds = unlist(fam.ids)[duplicated(unlist(fam.ids))]#none
detect.dup.id = function(x,id){
  ret = any(x == id)
  ret
}
for(i in 1:length(dup.inds)){
  print(which(sapply(fam.ids,detect.dup.id,id=dup.inds[i])))  
}

# the following part is manual...can I automate it?
do.manual <- function(dup.inds,fam.ids){
  if(any(dup.inds)){
    fam.ids_wo_dup = fam.ids
    # do the following manually
    #fam.ids_wo_dup[[1]] = unique(c(fam.ids_wo_dup[[1]],fam.ids_wo_dup[[466]],fam.ids_wo_dup[[533]]))
    print( length(fam.ids_wo_dup[[1]]) )#15
    # do the following manually
    #fam.ids_wo_dup[[466]] <- fam.ids_wo_dup[[533]] <- NA
    print( length(unlist(fam.ids_wo_dup)[!is.na(unlist(fam.ids_wo_dup))]))
    print( length(unlist(fam.ids_wo_dup[!is.na(fam.ids_wo_dup)])) )#2064
    fam.ids_wo_dup = fam.ids_wo_dup[!is.na(fam.ids_wo_dup)]#998
    print(length(unique(unlist(fam.ids_wo_dup))))
    fam.ids = fam.ids_wo_dup
  }
  rm(fam.ids_wo_dup)
  fam.ids
}
#
set.seed(20211109)
include.ids = c()
for(i in 1:length(fam.ids)){#59
  sample.ids = fam.ids[[i]]
  include.ids = c(include.ids,sample(sample.ids,size = 1))
  print(length(include.ids))
}
remove.ids = unlist(fam.ids)[!unlist(fam.ids)%in% include.ids]#n=61

brain_data = subset(d1,!eid %in% remove.ids);print(dim(brain_data))#7805
metabo_data = subset(d2,!eid %in% remove.ids);print(dim(metabo_data))#7805

print(all(brain_data$eid %in% metabo_data$eid))
print(all(metabo_data$eid %in% brain_data$eid))

rm(d1,d2)

# ------------------------------------------------------------------------------
#==============================================================================#
# get covariate data: NEED TO GET DIABETES + HTN data
#==============================================================================#
# 0. genetic sex 
f.c0 = file.path(ukbb_data_dir,'ukb_01-04-2020/41449/ukb41449.csv')
varnames.c0=c("genetic.sex"="22001-0.0")
c0 = extract_variables(f.c0,fieldID=varnames.c0,fieldName=names(varnames.c0))
head(c0,3)

# 1. fasting time, medication
###############################################################################
# coding
# 1	Cholesterol lowering medication #need the baseline value 
# 2	Blood pressure medication # for HTN (at imaging)
# 3	Insulin # for T1D and severe T2-DM (at imaging)
# 4	Hormone replacement therapy
# 5	Oral contraceptive pill or minipill
# -7	None of the above
# -1	Do not know
# -3	Prefer not to answer
###############################################################################
f.c1 = file.path(ukbb_data_dir,'ukb40646_20-02-2020/ukb40646.csv')
varnames.c1 = c(
  fasting_time = '74-0.0',
  medication01 = '6153-0.0',
  medication02 = '6153-0.1',
  medication03 = '6153-0.2'
#  medication1 = '6153-2.0',
#  medication2 = '6153-2.1',
#  medication3 = '6153-2.2'
)
c1 = extract_variables(f.c1,fieldID=varnames.c1,fieldName=names(varnames.c1))
c1 <- c1 %>% 
  mutate(medication01 = ifelse(is.na(medication01),0,medication01)) %>%
  mutate(medication01 = ifelse(medication01==-7,0,medication01)) %>%
  mutate(medication01 = ifelse(medication01 %in% c(-1,-3),NA,medication01)) %>%
  mutate(medication02 = ifelse(is.na(medication02),0,medication02)) %>%
  mutate(medication03 = ifelse(is.na(medication03),0,medication03)) %>%
  mutate(chol_lowering_med = ifelse(medication01==1,1,0)) %>%
  mutate(bp_med = ifelse((medication01==2|medication02==2),1,0)) %>%
  mutate(insulin = ifelse(((medication01==3|medication02==3)|medication03==3),1,0))
for(i in c('chol_lowering_med','bp_med','insulin')){
  print(tab <- table(c1[[i]]))
  print(table(c1[[i]],useNA = 'a'))
  print(prop.table(tab))
}

#2. age, sex, bmi, smoking, ethnicity
################################################################################
# smoking code
# 
# Coding	Meaning
# -3	Prefer not to answer
# 0	Never
# 1	Previous
# 2	Current
################################################################################
f.c2 = file.path(ukbb_data_dir,'ukb_01-04-2020/41448/ukb41448.csv')
varnames.c2 = c("age01" = '21022-0.0', "age0" = '21003-0.0',#are these the same?
                "age1" = '21003-2.0',
                "sex" = '31-0.0',
                "BMI0"="21001-0.0",
                "BMI1"="21001-2.0",
                "smoking0" = "20116-0.0",
                "smoking" = "20116-2.0",
                "SBP_manual0" = "93-0.0","SBP_auto0" = "4080-0.0",# at the time of blood sample
                "DBP_manual0" = "94-0.0","DBP_auto0" = "4079-0.0",# at the time of blood sample
                'ethnicity' = '21000-0.0')#may not need this
c2 = extract_variables(f.c2,fieldID=varnames.c2,fieldName=names(varnames.c2))
c2 <- c2 %>% 
  mutate(SBP0=ifelse(is.na(SBP_auto0)&!is.na(SBP_manual0),SBP_manual0,SBP_auto0)) %>% 
  mutate(DBP0=ifelse(is.na(DBP_auto0)&!is.na(DBP_manual0),DBP_manual0,DBP_auto0)) %>%
  mutate(highSBP = ifelse(SBP0>=140,1,0),highDBP=ifelse(DBP0>=90,1,0)) %>%
  mutate(highBP = ifelse(highSBP+highDBP>0,1,0)) %>%
  mutate(smoking0 = ifelse(smoking0==-3,NA,smoking0)) %>%
  mutate(current_smoking0 = ifelse(smoking0==1,1,0)) %>%
  mutate(smoking = ifelse(smoking==-3,NA,smoking)) %>%
  mutate(current_smoking = ifelse(smoking==1,1,0))
table(c2$smoking,c2$current_smoking,useNA='a')
c2 <- c2 %>% dplyr::select(-highSBP,-highDBP,-SBP_auto0,-DBP_auto0,-SBP_manual0,-DBP_manual0)
summary(c2)

#3. derive T2D
################################################################################
Diabetes <- c("E110", "E111", "E112", "E113", "E114", "E115", "E116",
              "E117", "E118", "E119", 
              "E121", "E128", "E129",
              "E131", "E132", "E133", "E134", "E135", "E136", "E137",
              "E138", "E139",
              "E140", "E141", "E142", "E143", "E144", "E145", "E146",
              "E147", "E148", "E149",
              "O241", "O243", "O244", "O249")
################################################################################

f.c3 = file.path(ukbb_data_dir,'ukb42388_18062020/ukb42388.csv')
c3= fread(f.c3)
varnames.c3 = names(c3)[str_detect(names(c3),"41270-0")]
c3 = subset(c3, select=c('eid',varnames.c3))
c3 <- c3 %>% mutate(T2D = apply(c3, 1, function(x)as.integer(
  any(grep(paste(Diabetes,collapse="|"),x)))))
c3 <- dplyr::select(c3, eid, T2D)#1-yes; 0-no

#4. merge cov data sets 
covdata = merge(c0,c1);print(dim(covdata))
covdata = merge(covdata,c2);print(dim(covdata))
covdata = merge(covdata,c3);print(dim(covdata))
covdata = subset(covdata ,eid %in% brain_data$eid);print(dim(covdata))
# define hypertension
covdata = covdata %>% mutate(HTN = ifelse(bp_med==1|highBP==1,1,0))
table(covdata$HTN,useNA='a')
#https://academic.oup.com/eurheartj/article/39/suppl_1/ehy563.3028/5080338
prop.table(table(covdata$HTN))#46% with HTN? (in the individuals with brain MRI)

brain_data = brain_data %>% dplyr::arrange(eid)
metabo_data = metabo_data %>% dplyr::arrange(eid)
covdata = covdata %>% dplyr::arrange(eid)

identical(brain_data$eid,metabo_data$eid)
identical(brain_data$eid,covdata$eid)

# write files ------------------------------------------------------------------
write_tsv(brain_data,"../data/ukb_brain_data.tsv")
write_tsv(metabo_data,"../data/ukb_metabo_data.tsv")
write_tsv(covdata,"../data/ukb_covdata.tsv")
