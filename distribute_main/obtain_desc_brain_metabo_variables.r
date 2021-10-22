brain_metabo_varnames = c(names(outcome_vec),names(brain_vec),names(cov_vec),metaboID_ordered)

# Combined
group='combined'
outfile = paste(metabo_outfolder,"desc_stats_",group,"_wi_brain_metabo.Rdata",sep="")
desc_stats_wi_brain_metabo <- save_desc_output(phenoOUT=pheno_OUT,
                                         varnames = brain_metabo_varnames,
                                         group = group)
save(desc_stats_wi_brain_metabo,file=outfile)

# Sex-stratified
## females
group='females'
outfile = paste(metabo_outfolder,"desc_stats_",group,"_wi_brain_metabo.Rdata",sep="")
desc_stats_wi_brain_metabo <- save_desc_output(phenoOUT=subset(pheno_OUT,sex=="F"),
                                         varnames = brain_metabo_varnames,
                                         group = group)
save(desc_stats_wi_brain_metabo,file=outfile)

## males
group='males'
outfile = paste(metabo_outfolder,"desc_stats_",group,"_wi_brain_metabo.Rdata",sep="")
desc_stats_wi_brain_metabo <- save_desc_output(phenoOUT=subset(pheno_OUT,sex=="M"),
                                         varnames = brain_metabo_varnames,
                                         group = group)
save(desc_stats_wi_brain_metabo,file=outfile)

# Statin-use stratified
## on statin 
group="on_statin"
outfile = paste(metabo_outfolder,"desc_stats_",group,"_wi_brain_metabo.Rdata",sep="")
tmp = subset(pheno_OUT,statin_use=="yes")
if(!is.na(other_lipid_lowering_med_use)){
  tmp=subset(tmp,other_lipid_lowering_med_use=="no")
}
desc_stats_wi_brain_metabo <- save_desc_output(phenoOUT=tmp,
                                         varnames = brain_metabo_varnames,
                                         group = group)
save(desc_stats_wi_brain_metabo,file=outfile)
rm(tmp)

## not on statin
group="not_on_statin"
outfile = paste(metabo_outfolder,"desc_stats_",group,"_wi_brain_metabo.Rdata",sep="")
tmp = subset(pheno_OUT,statin_use=="no")
if(!is.na(other_lipid_lowering_med_use)){
  tmp=subset(tmp,other_lipid_lowering_med_use=="no")
}
desc_stats_wi_brain_metabo <- save_desc_output(phenoOUT=tmp,
                                         varnames = brain_metabo_varnames,
                                         group = group)

save(desc_stats_wi_brain_metabo,file=outfile)
rm(tmp)