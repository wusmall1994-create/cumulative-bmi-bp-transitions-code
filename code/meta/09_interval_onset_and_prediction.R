options(stringsAsFactors=FALSE, survey.lonely.psu='adjust')
suppressPackageStartupMessages(library(survey))
set.seed(20261003)
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_dir <- if(length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork=FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
source_dir <- file.path(output_root, "dynamic_analysis")
out <- file.path(output_root, "submission_extensions")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
# Load only function assignments: never execute the historical output-writing workflow.
for(e in parse(file.path(root,'code/meta/07_methodological_sensitivity.R'))) {
  if(is.call(e) && identical(e[[1]],as.name('<-')) && is.call(e[[3]]) && identical(e[[3]][[1]],as.name('function'))) eval(e)
}
endpoint_labels <- c(hypertension_onset='Interval-level hypertension onset',untreated_bp_improvement='Lower untreated BP state')
effects<-list(); descriptives<-list(); flows<-list(); dests<-list(); checks<-list(); predictions<-list(); increments<-list(); calbins<-list()
fit_custom <- function(z,cohort,label,bp=FALSE,exposure='excess_bmi25_z',link='cloglog',offset=TRUE,endpoint='hypertension_onset',more=character()) {
  covars<-c(minimal_covars(cohort),if(bp)c('start_sbp10','start_dbp10'),more)
  z<-complete_model_data(z,cohort,endpoint,exposure,covars); z<-droplevels(z)
  covars<-covars[vapply(covars,function(v) length(unique(z[[v]]))>1,logical(1))]
  f<-as.formula(paste(endpoint,'~',paste(c(exposure,covars,if(offset)'offset(log(interval_years))'),collapse='+')))
  m<-svyglm(f,design=make_design(z,cohort),family=quasibinomial(link))
  b<-coef(m)[exposure]; s<-sqrt(vcov(m)[exposure,exposure])
  data.frame(cohort=cohort,analysis=label,scale=if(link=='logit')'OR' else 'cloglog ratio',log_tir=b,se=s,tir=exp(b),ci_low=exp(b-1.96*s),ci_high=exp(b+1.96*s),n_intervals=nrow(z),participants=length(unique(z$person_id)),events=sum(z[[endpoint]]))
}
metrics <- function(y,p,w) {
  keep<-w>0 & is.finite(p);y<-y[keep];p<-pmin(pmax(p[keep],1e-7),1-1e-7);w<-w[keep];w<-w/mean(w)
  o<-order(p);g<-match(p[o],unique(p[o]));w0<-rowsum(w[o]*(1-y[o]),g,reorder=FALSE)[,1];w1<-rowsum(w[o]*y[o],g,reorder=FALSE)[,1]
  auc<-sum(w1*(cumsum(w0)-w0/2))/(sum(w1)*sum(w0))
  lp<-qlogis(p)
  ci<-suppressWarnings(glm(y~1+offset(lp),weights=w,family=quasibinomial()))
  cs<-suppressWarnings(glm(y~lp,weights=w,family=quasibinomial()))
  c(AUC=auc,Brier=weighted.mean((y-p)^2,w),calibration_intercept=unname(coef(ci)[1]),calibration_slope=unname(coef(cs)[2]))
}
for(cohort in c('HRS','CHNS','ELSA')) {
  cat('Cohort',cohort,'\n');flush.console()
  d<-prepare_data(cohort)
  h<-read.csv(file.path(out,paste0(tolower(cohort),'_history.csv')))
  h$prior_positive<-tolower(as.character(h$prior_positive))=='true'
  h$history_observed<-tolower(as.character(h$history_observed))=='true'
  d<-merge(d,h,by='person_id',all.x=TRUE,sort=FALSE)
  d$excess_per10<-d$excess_bmi25/10
  cov<-minimal_covars(cohort)
  z<-complete_model_data(d,cohort,'hypertension_onset','excess_bmi25_z',cov)
  zb<-complete_model_data(z,cohort,'hypertension_onset','excess_bmi25_z',c(cov,'start_sbp10','start_dbp10'))
  # Reproduce separate transition coefficient, then retain the published joint-fit SE for main results.
  effects[[length(effects)+1]]<-fit_custom(z,cohort,'primary_reproduction')
  effects[[length(effects)+1]]<-fit_custom(zb,cohort,'BP_adjusted')
  effects[[length(effects)+1]]<-fit_custom(zb,cohort,'primary_BP_complete')
  effects[[length(effects)+1]]<-fit_custom(z,cohort,'common_unit',exposure='excess_per10')
  effects[[length(effects)+1]]<-fit_custom(zb,cohort,'common_unit_BP',TRUE,'excess_per10')
  # Correct explicit BP flag for this named model.
  effects[[length(effects)-3]]<-fit_custom(zb,cohort,'BP_adjusted',TRUE)
  effects[[length(effects)+1]]<-fit_custom(z,cohort,'no_offset',offset=FALSE)
  effects[[length(effects)+1]]<-fit_custom(z,cohort,'logistic',link='logit',offset=FALSE)
  initial<-d[order(d$person_id,d$start_year),];initial<-initial[!duplicated(initial$person_id),]
  initial_ids<-initial$person_id[initial$origin %in% c(1,2)]
  first_ids<-initial$person_id[initial$origin %in% c(1,2) & !is.na(initial$prior_positive) & !initial$prior_positive & initial$history_observed]
  # Stop at first observed hypertension (including intervals unusable in the fitted model).
  first_event<-aggregate(end_year~person_id,data=d[d$destination %in% c(3,4,5),],FUN=min)
  first_time<-first_event$end_year[match(d$person_id,first_event$person_id)]
  first_keep<-is.na(first_time)|d$start_year<first_time
  first<-d[d$person_id %in% first_ids & first_keep,]
  effects[[length(effects)+1]]<-fit_custom(first,cohort,'first_observed')
  effects[[length(effects)+1]]<-fit_custom(first,cohort,'first_observed_BP',TRUE)
  effects[[length(effects)+1]]<-fit_custom(d[d$person_id %in% initial_ids & first_keep,],cohort,'first_followup_only')
  effects[[length(effects)+1]]<-fit_custom(initial,cohort,'initial_interval')
  flows[[cohort]]<-data.frame(cohort=cohort,exposure_complete=length(unique(d$person_id)),observed_transition=length(unique(d$person_id[!is.na(d$origin)&!is.na(d$destination)])),onset_eligible=length(unique(d$person_id[!is.na(d$hypertension_onset)])),primary_people=length(unique(z$person_id)),primary_intervals=nrow(z),events=sum(z$hypertension_onset),history_eligible_people=length(unique(first$person_id)),history_positive_at_initial=sum(initial$origin %in% c(1,2)&initial$prior_positive,na.rm=TRUE))
  zz<-z[order(z$person_id,z$start_year),];zz<-zz[!duplicated(zz$person_id),]
  excluded<-initial[!initial$person_id %in% zz$person_id,]
  for(group in c('primary','not_in_primary')) {
    v<-if(group=='primary')zz else excluded
    for(nm in c('age_start10','female','start_sbp','start_dbp','simple_bmi','bmi_mean','excess_bmi25')) {
      descriptives[[length(descriptives)+1]]<-data.frame(cohort=cohort,group=group,variable=nm,n=sum(!is.na(v[[nm]])),mean=mean(v[[nm]],na.rm=TRUE),sd=sd(v[[nm]],na.rm=TRUE))
    }
  }
  for(state in c(1,2,3,4,5)) dests[[length(dests)+1]]<-data.frame(cohort=cohort,destination=state,n=sum(z$destination==state))
  for(type in c('untreated_htn','treated_htn')) {
    z$dest_event<-as.integer(if(type=='untreated_htn')z$destination==3 else z$destination %in% c(4,5))
    effects[[length(effects)+1]]<-fit_custom(z,cohort,type,endpoint='dest_event')
  }
  for(iv in unique(z$interval)) checks[[length(checks)+1]]<-data.frame(cohort=cohort,interval=as.character(iv),durations=paste(sort(unique(z$interval_years[z$interval==iv])),collapse=','))
  # Conditional association beyond mean BMI, not a duration-effect estimate.
  effects[[length(effects)+1]]<-fit_custom(zb,cohort,'mean_BMI_plus_burden',TRUE,more='bmi_mean_z')
  # Common sample, identical participant-grouped five-fold split for every prediction model.
  pdat<-complete_model_data(d,cohort,'hypertension_onset','excess_bmi25_z',c(cov,'start_sbp10','start_dbp10','simple_bmi','bmi_mean_z'))
  pdat<-droplevels(pdat);covp<-c(cov,'start_sbp10','start_dbp10');covp<-covp[vapply(covp,function(v)length(unique(pdat[[v]]))>1,logical(1))]
  specs<-list(clinical=character(),clinical_burden='excess_bmi25_z',clinical_current='simple_bmi',clinical_mean='bmi_mean_z',clinical_current_burden=c('simple_bmi','excess_bmi25_z'),clinical_mean_burden=c('bmi_mean_z','excess_bmi25_z'))
  X<-lapply(specs,function(s)model.matrix(as.formula(paste('~',paste(c(covp,s),collapse='+'))),pdat))
  ids<-unique(pdat$person_id);fold<-sample(rep(1:5,length.out=length(ids)));fv<-fold[match(pdat$person_id,ids)]
  y<-pdat$hypertension_onset;w<-pdat$model_weight;off<-log(pdat$interval_years)
  cv<-function(ww) {
    pp<-matrix(NA_real_,nrow(pdat),length(X),dimnames=list(NULL,names(X)))
    for(j in seq_along(X))for(k in 1:5) {
      tr<-fv!=k & ww>0;te<-fv==k
      ff<-suppressWarnings(glm.fit(X[[j]][tr,,drop=FALSE],y[tr],weights=ww[tr]/mean(ww[tr]),offset=off[tr],family=quasibinomial('cloglog'),control=glm.control(maxit=100)))
      if(!ff$converged)stop('CV convergence failed')
      bb<-ff$coefficients;bb[is.na(bb)]<-0
      pp[te,j]<-1-exp(-exp(pmin(30,as.vector(X[[j]][te,,drop=FALSE]%*%bb)+off[te])))
    }
    mm<-sapply(seq_along(X),function(j)metrics(y,pp[,j],ww));colnames(mm)<-names(X)
    list(metrics=mm,p=pp)
  }
  point<-cv(w)
  B<-200L;boot<-array(NA_real_,c(4,length(X),B))
  strata<-if(cohort=='HRS')pdat$strata_value else rep(1,nrow(pdat))
  cluster<-interaction(strata,pdat$psu_id,drop=TRUE)
  cltab<-unique(data.frame(cl=as.character(cluster),st=strata))
  for(b in seq_len(B)) {
    mult<-setNames(rep(0,nrow(cltab)),cltab$cl)
    for(st in unique(cltab$st)) {cc<-cltab$cl[cltab$st==st];tb<-table(sample(cc,length(cc),replace=TRUE));mult[names(tb)]<-as.numeric(tb)}
    wb<-w*mult[as.character(cluster)]
    rr<-try(cv(wb),silent=TRUE);if(!inherits(rr,'try-error'))boot[,,b]<-rr$metrics
    if(b%%50==0){cat(cohort,'validation bootstrap',b,'/',B,'\n');flush.console()}
  }
  for(j in seq_along(X))for(m in 1:4) {
    vv<-boot[m,j,];ci<-quantile(vv,c(.025,.975),na.rm=TRUE)
    predictions[[length(predictions)+1]]<-data.frame(cohort=cohort,model=names(X)[j],metric=rownames(point$metrics)[m],estimate=point$metrics[m,j],low=ci[1],high=ci[2],successful=sum(is.finite(vv)),n_intervals=nrow(pdat),participants=length(ids),events=sum(y))
  }
  for(pair in list(c(1,2),c(3,5),c(4,6)))for(m in 1:2) {
    vv<-boot[m,pair[2],]-boot[m,pair[1],];ci<-quantile(vv,c(.025,.975),na.rm=TRUE)
    increments[[length(increments)+1]]<-data.frame(cohort=cohort,comparison=paste(names(X)[pair],collapse=' -> '),metric=rownames(point$metrics)[m],estimate=point$metrics[m,pair[2]]-point$metrics[m,pair[1]],low=ci[1],high=ci[2],successful=sum(is.finite(vv)))
  }
  for(j in seq_along(X)) {
    pp<-point$p[,j];grp<-cut(rank(pp,ties.method='first'),breaks=seq(0,length(pp),length.out=6),include.lowest=TRUE,labels=FALSE)
    for(g in 1:5)calbins[[length(calbins)+1]]<-data.frame(cohort=cohort,model=names(X)[j],group=g,n=sum(grp==g),predicted=weighted.mean(pp[grp==g],w[grp==g]),observed=weighted.mean(y[grp==g],w[grp==g]))
  }
  saveRDS(list(point=point,bootstrap=boot,fold=fv,person_id=pdat$person_id),file.path(out,paste0(tolower(cohort),'_validation_local.rds')))
  write.csv(do.call(rbind,effects),file.path(out,'effects.csv'),row.names=FALSE)
  write.csv(do.call(rbind,predictions),file.path(out,'prediction_metrics.csv'),row.names=FALSE)
  write.csv(do.call(rbind,increments),file.path(out,'prediction_increments.csv'),row.names=FALSE)
}
ef<-do.call(rbind,effects)
meta<-lapply(split(ef,ef$analysis),function(q){if(nrow(q)!=3)return(NULL);cbind(analysis=q$analysis[1],meta_estimates(q))})
write.csv(do.call(rbind,meta),file.path(out,'meta.csv'),row.names=FALSE)
write.csv(do.call(rbind,descriptives),file.path(out,'descriptives.csv'),row.names=FALSE)
write.csv(do.call(rbind,flows),file.path(out,'sample_flow.csv'),row.names=FALSE)
write.csv(do.call(rbind,dests),file.path(out,'destinations.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(out,'duration_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,calbins),file.path(out,'calibration_bins.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
cat('ANALYSIS COMPLETE\n')


