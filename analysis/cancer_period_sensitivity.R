# Exploratory diagnosis-period check for every named histology-site group.
# Preserve primary estimates; use shared patient resamples for paired contrasts.
suppressPackageStartupMessages({library(data.table);library(survival)})
p<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study");out<-file.path(p,'revision/cancer');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-fread(file.path(p,'primary_analysis_cohort.csv'));x[,id:=.I]
src<-readLines('analysis/population_primary.R',warn=FALSE)
i<-grep('^prep<-function',src);j<-grep('^groups<-list',src)
eval(parse(text=paste(src[i:(j-1)],collapse='\n')))
groups<-unique(fread(file.path(p,'primary_conditional_risks.csv'))[kind=='Subtype x site',group])
cache<-list()
for(g in groups){a<-strsplit(g,' | ',fixed=TRUE)[[1]];z<-x[subtype==a[1]&site==a[2]]
 for(period in c('Primary 2004-2023','Restricted 2004-2015')){
  zz<-if(period=='Restricted 2004-2015')z[year<=2015] else z
  cache[[length(cache)+1]]<-prep(zz,60,g,'Subtype x site',period)
 }}
point<-t(vapply(cache,calc,numeric(8),w=rep(1,nrow(x))))
err<-0
for(k in seq_along(cache)){
 cc<-cache[[k]];zz<-x[cc$id]
 fit<-survfit(Surv(time-60,factor(event,levels=0:3))~1,data=zz)
 ss<-summary(fit,times=36,extend=TRUE)$pstate[1,2:4]
 err<-max(err,max(abs(ss-point[k,1:3])))
}
stopifnot(err<1e-10)
B<-1000L;set.seed(20261017);boot<-array(NA_real_,c(length(cache),8,B))
for(b in seq_len(B)){
 w<-tabulate(sample.int(nrow(x),nrow(x),replace=TRUE),nbins=nrow(x))
 boot[,,b]<-t(vapply(cache,calc,numeric(8),w=w))
 if(b%%250==0)cat('Joint period bootstrap',b,'/',B,'\n')
}
res<-rbindlist(lapply(seq_along(cache),function(k){cc<-cache[[k]]
 rbindlist(lapply(c(1,2,4),function(v)data.table(group=cc$group,period=cc$scenario,landmark=5,horizon=3,n=cc$n,
 cancer_events=sum(cc$event==1&cc$relative<=36),other_events=sum(cc$event==2&cc$relative<=36),
 observed_beyond_horizon=sum(cc$relative>36),outcome=c('Cancer','Other','Unknown','Difference')[v],
 estimate=point[k,v],lower=quantile(boot[k,v,],.025),upper=quantile(boot[k,v,],.975),bootstrap=B)))
}))
delta<-rbindlist(lapply(seq_along(groups),function(k){a<-2*k-1;b<-2*k
 rbindlist(lapply(c(1,2,4),function(v){vv<-boot[b,v,]-boot[a,v,]
 data.table(group=groups[k],outcome=c('Cancer','Other','Unknown','Difference')[v],contrast='Restricted minus primary, percentage points',
 estimate=100*(point[b,v]-point[a,v]),lower=100*quantile(vv,.025),upper=100*quantile(vv,.975))}))
}))
fwrite(res,file.path(out,'joint_period_risks.csv'));fwrite(delta,file.path(out,'joint_period_paired_differences.csv'))
writeLines(c('Exploratory analysis added during Cancer submission preparation, 2026-10-17.',
 'All 12 named joint groups; five-year survivors; future three-year risks; diagnoses restricted to 2004-2015.',
 '1000 shared resamples of the full primary cohort; percentile pointwise intervals; no multiplicity correction.',
 'Paired differences preserve overlap. Period restriction changes case mix and is not independent validation.',
 paste('Independent survival::survfit maximum absolute error',format(err,scientific=TRUE)),capture.output(sessionInfo())),file.path(out,'period_analysis_log.txt'))
print(res[period=='Restricted 2004-2015'&outcome=='Cancer',.(group,n,cancer_events,estimate,lower,upper)])
