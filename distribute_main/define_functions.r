replace.missString.wi.NA <- function(x,na.string="NA"){
  ind.na = x %in% na.string
  ind.na = is.na(x) | (ind.na & !is.na(ind.na))
  x[ind.na] <- NA
  x
}
descriptive_stats_continuous <- function(x){
  x.mean <- mean(x,na.rm=T)
  x.median <- median(x,na.rm = T)
  x.min <- min(x,na.rm=T)
  x.max <- max(x,na.rm=T)
  x.q.0.25 <- quantile(x,probs=0.25,na.rm = T)
  x.q.0.75 <- quantile(x,probs=0.75,na.rm = T)
  x.sd <- sd(x,na.rm=T)
  x.var <- var(x,na.rm=T)
  x.skew <- skewness(x,na.rm = T,type = 1) #typical
  x.kurt <- kurtosis(x,na.rm=T,type = 1)#typical 
  n_complete <- sum(!is.na(x))
  ret = c(x.mean,x.median,x.min,x.max,x.q.0.25,x.q.0.75,x.sd,x.var,x.skew,x.kurt,n_complete)
  names(ret) = c("mean","median","min","max","qu1","qu3","sd","var","skew","kurt","n_complete")
  ret = format(ret,scientific = T, digits = 4)
  ret
}#descriptive_stats_continuous

descriptive_stats_discrete <- function(x){
  ret = table(x)
  ret
}#descriptive_stats_discrete

distr.simple<-function(dset,vec,tle){
  # plotting
  dset<-subset(dset, select=vec)
  dset<-na.exclude(dset)
  obj<-eval(parse(text=paste("dset",vec,sep="$")))
  truehist(obj, col="grey80",xlab = "",main="",cex.axis=1)
  title(main=tle)
  
  pdist<-density(obj, na.rm=TRUE)
  #lines(density(obj, na.rm=TRUE))
  m<-mean(obj, na.rm=TRUE)
  md<-median(obj, na.rm=TRUE)
  d<-sd(obj, na.rm=TRUE)
  ys<-c(0,max(pdist$y))
  xs<-c(m,m)
  #lines(xs, ys, lty = 1, lwd = 1.5)  
  xs<-c(md,md)
  #lines(xs, ys, lty = 2, lwd = 1.5)    
  plot(function(x) dnorm(x, mean=m, sd=d, log = FALSE), min(obj, na.rm=T),max(obj, na.rm=T), add=TRUE, lty=2)
  #legend("topleft",c("Normal","Density"), lty=c(2,1),bty="n")
}#calculate_missing_coun

calculate_missing_count <- function(x){
  ret = sum(is.na(x))
  ret
}#calculate_missing_coun

calculate_zero_count <- function(x){
  ret = x==0
  ret = sum(ret,na.rm=T)
  ret
}#calculate_zero_count

save_desc_output <- function(phenoOUT,varnames,group){
  tmp = subset(phenoOUT,select=varnames)
  coltypes = unlist(lapply(tmp, mode))
  ret_cont = apply(tmp[,coltypes=="numeric",drop=F],2,descriptive_stats_continuous)
  ret_discrete = lapply(tmp[,coltypes!="numeric",drop=F],descriptive_stats_discrete)
  ret = list(continuous=ret_cont,discrete=ret_discrete)
  ret
}#save_desc_output

generate_histograms <- function(metabo_OUT,summarydata_Rdata,metabo_missing_zero_info,metabo_columns){
  #create a sub-directory to store the pdf files
  dir.create(paste(metabo_outfolder,"plots/",sep=''))
  
  d <- metabo_OUT
  names(metabo_missing_zero_info)[1] <- "Abbreviation2"
  load(summarydata_Rdata)#load 'desc_stats_wi_metabo'
  metabo_names = colnames(desc_stats_wi_metabo[[1]][,metabo_columns])
  summarydata <- t(desc_stats_wi_metabo[[1]][,metabo_columns])
  summarydata <- apply(summarydata,2,as.numeric)
  summarydata <- data.frame(Abbreviation2=metabo_names,summarydata)
  summarydata <- merge(metabo_missing_zero_info[,c("Abbreviation2","N",'missingness_pct')],
                       summarydata,by="Abbreviation2",sort=F)
  summarydata$Abbreviation2 <- factor(summarydata$Abbreviation2,levels=metabo_columns)
  
  myvars <- paste(summarydata$Abbreviation2[which(summarydata$N>=2)])
  i=1
  for(l in seq(1,length(myvars), by=60)){
    temp_vars <- myvars[l:(l+59)]
    temp <- d[,which(colnames(d) %in% temp_vars)]
    tmp.colnames = colnames(temp)
    temp = data.frame(temp)
    names(temp) = tmp.colnames
    temp_long <- melt(temp)
    temp_info <- summarydata[which(summarydata$Abbreviation2 %in% temp_vars),]
    names(temp_info)[1] <- "variable"
    
    myplot <- ggplot()+
      geom_histogram(data=temp_long, aes(x=value), fill="lightsteelblue2", bins=30)+
      geom_text(data=temp_info,
                aes(x=Inf, y=Inf, label=ifelse(mean<100, paste("mean", formatC(mean, format="f", digits=3)), paste("mean", round(mean, digits=0)))),
                hjust=1, vjust=2, size=2.3)+
      geom_text(data=temp_info,
                aes(x=Inf, y=Inf, label=ifelse(median<100, paste("median", formatC(median, format="f", digits=3)), paste("median", round(median, digits=0)))),
                hjust=1, vjust=3.5, size=2.3)+
      geom_text(data=temp_info,
                aes(x=Inf, y=Inf, label=paste("N", N)),
                hjust=1, vjust=5, size=2.3)+
      scale_y_continuous(expand=c(0,0))+
      facet_wrap(~variable, scales="free", ncol=10)+
      theme(strip.background=element_blank(),
            strip.text=element_text(size=8, face="bold"),
            axis.text.y=element_text(size=7),
            axis.text.x=element_text(size=7, angle=30, hjust=1, vjust=1),
            axis.title=element_text(size=8),
            plot.title=element_text(size=9, face="bold"))
    nvars <- length(temp_vars[which(!is.na(temp_vars))])
    h <- ifelse(nvars<11, 2.1, ifelse(nvars>=11 & nvars<21, 3.2, ifelse(nvars>=21 & nvars<31, 4.3, ifelse(nvars>=31 & nvars<41, 5.4, ifelse(nvars>=41 & nvars<51, 6.5, ifelse(nvars>=51, 7.6, NA))))))
    ggsave(paste0(paste(metabo_outfolder,"plots/",sep=''),"histo",i,".pdf"), myplot, width=15, height=h)
    print(paste("variables",l, "-",(l+59),"done")) ; i=i+1  
  }
}#generate_histograms

fit_mod <- function(data, yname, xname, mod){
  # extract terms from the provided model
  require(stringr)
  if(yname!="logWMH"){# if brain outcome is not WMH, do not adjust for ICV or brainsize
    mod = str_remove(mod," [+] ICV_or_BrainSize")
  }
  regmod = as.formula(mod)
  
  mod.terms = row.names(attr(terms(regmod,attr="var"),"factors"))
  vnames=c(yname,xname,mod.terms[-c(1,2)])
  if(any(vnames=="statin_use")){
    if(!is.na(other_lipid_lowering_med_use)){
      vnames=c(vnames,'other_lipid_lowering_med_use')
    }
  }
  
  # subsetting data just to include the variables of interest
  anal_data = subset(data,select=vnames)
  anal_data = na.omit(anal_data)
  anal_data$y <- anal_data[,yname]
  anal_data$x <- anal_data[,xname]
  n_sample_x = sum(!is.na(anal_data$x))
  
  if(any(vnames=="sex")){
    #cat("Recoding \'Sex\': 0-male and 1-female\n")
    Sex_tmp <- anal_data$sex
    Sex_tmp[anal_data$sex=="F"] <- 1
    Sex_tmp[anal_data$sex=="M"] <- 0
    anal_data$sex <- Sex_tmp 
  }
  if(any(vnames=="statin_use")){
    #cat("Recoding \'statin_use\': 0-off statin and 1-on statin\n")
    # cat("Recoding \'OnOtherMed\': 0-off other lipid-lowering medication and 1-on other lipid-lowering medication\n")
    statin_use_tmp <- anal_data$statin_use
    statin_use_tmp[anal_data$statin_use=="yes"] <- 1
    statin_use_tmp[anal_data$statin_use=="no"] <- 0
    anal_data$statin_use <- statin_use_tmp
    if(!is.na(other_lipid_lowering_med_use)){
      OnOtherMed_tmp <- anal_data[,"other_lipid_lowering_med_use"]
      OnOtherMed_tmp[anal_data[,"other_lipid_lowering_med_use"]=="yes"] <- 1
      OnOtherMed_tmp[anal_data[,"other_lipid_lowering_med_use"]=="no"] <- 0
      anal_data[,"other_lipid_lowering_med_use"] <- OnOtherMed_tmp
      anal_data <- subset(anal_data,other_lipid_lowering_med_use==0)
    }
  }
  var.ind = NULL
  vnames = names(anal_data)[names(anal_data)!='other_lipid_lowering_med_use']
  for(varname in vnames){
    x <- anal_data[[varname]]
    if(is.numeric(x)) {
      ind = var(x) > 0
    }else{
      ind = length(table(x))>1
      if(!ind) x<- as.numeric(as.factor(x))
      anal_data[[varname]] <- x
    }
    var.ind = c(var.ind,ind)
  }
  rm(x)
  
  n_sample_ok <- (n_sample_x > (length(mod.terms)+10))#
  if(!n_sample_ok) warning("Only ",n_sample_x," has information on ",xname,": less than length(mod.terms)+10.\n")
  x_y_var <- all(c('y','x') %in% vnames[var.ind])# y and x are variable
  if(!x_y_var){
    if(!'y'%in% vnames[var.ind]) warning(yname," is not variable.\n",sep='') 
    if(!'x'%in% vnames[var.ind]) warning(xname," is not variable.\n",sep="") 
  } 
  ## print(xname)##
  
  if(n_sample_ok & x_y_var){
    # fit the linear regression model to the data
    fit <- lm(regmod, data=anal_data)
    
    # check if there is any influential points (re-fit without them)
    outtest = outlierTest(fit)
    ret = data.frame(sub_id = names(outtest$rstudent),
                     rstudent = outtest$rstudent,
                     p = outtest$p,
                     bonf.p = outtest$bonf.p, 
                     signif = outtest$signif,
                     cutoff = outtest$cutoff,
                     stringsAsFactors = F)
    bonf.p.0.05 <- outtest$bonf.p < 0.05
    if(any(bonf.p.0.05)){
      # cat(yname,"-",xname,": Removing ",sum(bonf.p.0.05)," outliers\n")
      outliers = as.numeric(ret$sub_id)
      anal_data_screen = anal_data[-outliers,]
      fit <- lm(regmod,data=anal_data_screen)
    }#if(any(bonf.p.0.05))
  }else{
    fit <- NA
  }
  fit
}#fit_mod

fit_mod_tmp <- function(x,data,yname,mod){
  ret = fit_mod(data,yname,x,mod)
  ret
}#fit_mod_tmp

fit_get_results <- function(modName,group){
  mod = mod.list[[modName]]
  
  if(group == "combined"){
    anal_subset = dta
  }else if(group == "males"){
    anal_subset = subset(dta,sex=="M")
  }else if(group == "females"){
    anal_subset = subset(dta,sex=="F")
  }else if(group == "on_statin"){
    if(!is.na(other_lipid_lowering_med_use)){
      anal_subset = subset(dta,statin_use=="yes" & other_lipid_lowering_med_use=="no")
    }else{
      anal_subset = subset(dta,statin_use=="yes")
    }
  }else if(group == "not_on_statin"){
    if(!is.na(other_lipid_lowering_med_use)){
    anal_subset = subset(dta,statin_use=="no" & other_lipid_lowering_med_use=="no")
    }else{
      anal_subset = subset(dta,statin_use=="no")
    }
  }else{
    stop('Error [fit_get_results]: not a recognized group.\n')
  }
  
  fit_res <- vector(mode="list",length=length(outcome_vec))
  names(fit_res) <- names(outcome_vec)
  
  #fitting
  for(yi in 1:length(outcome_vec)){
    yname=names(outcome_vec)[yi]
    fit_res[[yi]] <- lapply(metaboID_ordered,fit_mod_tmp,data=anal_subset,yname=yname,mod=mod)
    names(fit_res[[yi]]) <- metaboID_ordered
  }
  
  #organize results
  res_modi <- NULL
  for(yi in 1:length(outcome_vec)){
    res = NULL
    nsample <- xname <- c()
    for(xi in 1:length(metaboID_ordered)){
      fit_avail1 = any(!is.na(fit_res[[yi]][[xi]])) 
      if(fit_avail1){
        fit_avail2 = !is.na(coef(fit_res[[yi]][[xi]])[["x"]])
        if(fit_avail2){ 
          res = rbind(res,summary(fit_res[[yi]][[xi]])$coef["x",])
          nsample = c(nsample,nobs(fit_res[[yi]][[xi]]))
          xname = c(xname,names(fit_res[[yi]])[xi])
        }#if(fit_avail2)
      }#if(fit_avail1)
    }#for(xi in 1:length(metaboID_ordered))
    if(!is.null(res)){
      res = data.frame(res)
      res = res[,names(res)!="t.value"]
      names(res) = c("BETA","SE","P")
      res = data.frame(MetaboID=xname,res)
    }else{#if(!is.null(res))
      res = data.frame(MetaboID=metaboID_ordered)
    }#else
    res$N <- nsample
    res$pheno <- names(fit_res)[yi]
    res_modi <- rbind(res_modi,res)
  }#for(yi in 1:length(outcome_vec))
  res_modi$model<- mod
  res_modi$modelName <- modName
  res_modi$group <- group
  
  fname=paste(resdir,modName,"_",group,sep="")
  fname=paste(fname,".tsv",sep="")
  write_tsv(res_modi,fname)
  
  res_modi
  
}#fit_get_results
