options(stringsAsFactors=FALSE,survey.lonely.psu='adjust')
suppressPackageStartupMessages(library(survey))
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_dir <- if(length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork=FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
source_dir <- file.path(output_root, "dynamic_analysis")
out <- file.path(output_root, "submission_extensions")
for(e in parse(file.path(root,'code/meta/07_methodological_sensitivity.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e)
for(e in parse(file.path(root,'code/meta/09_interval_onset_and_prediction.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('metrics')))eval(e)
set.seed(20261004)
pm<-read.csv(file.path(out,'prediction_metrics.csv'));di<-read.csv(file.path(out,'prediction_increments.csv'))
write.csv(pm,file.path(out,'prediction_metrics_refit_bootstrap.csv'),row.names=FALSE)
write.csv(di,file.path(out,'prediction_increments_refit_bootstrap.csv'),row.names=FALSE)
for(cohort in c('HRS','CHNS','ELSA')) {
 d<-prepare_data(cohort);h<-read.csv(file.path(out,paste0(tolower(cohort),'_history.csv')));d<-merge(d,h,by='person_id',all.x=TRUE,sort=FALSE)
 d<-complete_model_data(d,cohort,'hypertension_onset','excess_bmi25_z',c(minimal_covars(cohort),'start_sbp10','start_dbp10','simple_bmi','bmi_mean_z'))
 obj<-readRDS(file.path(out,paste0(tolower(cohort),'_validation_local.rds')))
 stopifnot(identical(d$person_id,obj$person_id))
 strata<-if(cohort=='HRS')d$strata_value else rep(1,nrow(d));cluster<-interaction(strata,d$psu_id,drop=TRUE)
 cltab<-unique(data.frame(cl=as.character(cluster),st=strata));boot<-array(NA_real_,c(4,6,1000))
 for(b in 1:1000) {
  mult<-setNames(rep(0,nrow(cltab)),cltab$cl)
  for(st in unique(cltab$st)){cc<-cltab$cl[cltab$st==st];tb<-table(sample(cc,length(cc),replace=TRUE));mult[names(tb)]<-as.numeric(tb)}
  ww<-d$model_weight*mult[as.character(cluster)]
  boot[,,b]<-sapply(1:6,function(j)metrics(d$hypertension_onset,obj$point$p[,j],ww))
 }
 for(j in 1:6)for(k in 1:4){ix<-pm$cohort==cohort&pm$model==colnames(obj$point$metrics)[j]&pm$metric==rownames(obj$point$metrics)[k];ci<-quantile(boot[k,j,],c(.025,.975));pm$low[ix]<-ci[1];pm$high[ix]<-ci[2];pm$successful[ix]<-1000}
 for(pair in list(c(1,2),c(3,5),c(4,6)))for(k in 1:2){ix<-di$cohort==cohort&di$comparison==paste(colnames(obj$point$metrics)[pair],collapse=' -> ')&di$metric==rownames(obj$point$metrics)[k];ci<-quantile(boot[k,pair[2],]-boot[k,pair[1],],c(.025,.975));di$low[ix]<-ci[1];di$high[ix]<-ci[2];di$successful[ix]<-1000}
 cat(cohort,'conditional validation intervals complete\n')
}
write.csv(pm,file.path(out,'prediction_metrics.csv'),row.names=FALSE)
write.csv(di,file.path(out,'prediction_increments.csv'),row.names=FALSE)


