library(tidyverse)
library(gamm4)
library(MuMIn)

# function for gamm and r-squared output
# data contains columns for the dependent variable (dv), independent variables (ivs [usually k and m]) ,
#          covariates (covs--demographic info), ABCD site id (site_id_l), and family id (rel_family_id)
#          name gives the filepath and name to name the output csv (not including csv)
gamm.abcd<-function(dv,ivs,covs,data,name){
  all_complete=complete.cases(data %>% select(dv,ivs,covs))
  data=data[all_complete,]
  
  gamm<-gamm4(as.formula(paste(paste0("scale(",dv,")")," ~ ",paste0("scale(",ivs,")",collapse='+')," + ",
                               paste(covs,collapse=' + '))),
              random = ~(1|site_id_l/rel_family_id), 
              data = data, na.action='na.omit')
  sum<-as.data.frame(summary(gamm$gam)$p.table)
  sum$rsqrd<-NA
  sum$n<-gamm$gam$df.null
  
  for (i in ivs){
    if(length(ivs>1)){
      null<-gamm4(as.formula(paste(paste0("scale(",dv,")")," ~ ",paste0("scale(",ivs[!ivs%in%i],")",collapse=' + ')," + ",
                                   paste(covs,collapse=' + '))),
                  random = ~(1|site_id_l/rel_family_id), 
                  data = data, na.action='na.omit')} else {
                    null<-gamm4(as.formula(paste(dv," ~ ",
                                                 paste(covs,collapse=' + '))),
                                random = ~(1|site_id_l/rel_family_id), 
                                data = data, na.action='na.omit')}
    
    sum[i,]$rsqrd<-round(as.numeric(r.squaredLR(gamm$mer,null$mer)),5)
  }
  
  write.csv(sum,file=paste0(name,".csv"))
  tmp=summary(gamm$gam)
  ps=tmp$p.pv
  r=tmp$r.sq
  names(r)="r.sq"
  stats=c(ps,r)
  return(stats)
  
}