suppressPackageStartupMessages({library(data.table);library(survival)})
p<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study");out<-file.path(p,'survivorship')
x<-fread(file.path(p,'primary_analysis_cohort.csv'),colClasses=c(patient_id='character'))
x[,id:=.I]
r<-fread(Sys.getenv("SEER_EXPORT", unset="inputs/seer_export.tsv"),colClasses='character');setnames(r,names(r),paste0('v',1:23))
r[,rank:=fifelse(v10=='One primary only',1L,suppressWarnings(as.integer(sub('^(\\d+).*','\\1',v10))))];r[,yr:=as.integer(v2)];r[,ord:=.I];setorder(r,v1,rank,yr,ord)
r<-r[!duplicated(v1),.(patient_id=v1,cod_label=v23,raw_sequence=v10,raw_year=yr,vital=v20)]
x<-merge(x,r,by='patient_id',sort=FALSE)
stopifnot(nrow(x)==12219,all(x$sequence==x$raw_sequence),all(x$year==x$raw_year),all((x$vital=='Dead')==(x$event>0)))
map<-fread(file.path(out,'cod_label_mapping.csv'));x<-merge(x,map,by='cod_label',sort=FALSE);setorder(x,id)
stopifnot(!anyNA(x$cod_code),all(x[event==0,cod_code]==0),all(x[event==3,cod_code] %in% c(41000,99999)))
# Preserve the SEER attribution; do not label malignant-coded other death as a second primary.
x[,detailed:=fcase(event==0,0L,event==1,1L,event==3,6L,certificate_group=='Malignant cancer',2L,certificate_group=='Cardiovascular',3L,certificate_group=='Other specified nonmalignant',4L,certificate_group=='Residual or ill-defined',5L,default=6L)]
x[,certificate:=fcase(event==0,0L,certificate_group=='Malignant cancer',1L,certificate_group=='Cardiovascular',2L,certificate_group=='Other specified nonmalignant',3L,certificate_group=='Residual or ill-defined',4L,default=5L)]
names_d<-c('SEER-attributed cancer','Malignant-coded other death','Cardiovascular','Other specified nonmalignant','Residual or ill-defined','Unknown COD')
names_c<-c('All malignant cancer','Cardiovascular','Other specified nonmalignant','Residual or ill-defined','Unknown COD')
fwrite(x[event>0,.(n=.N),by=.(sequence,site,subtype,event,cod_label,cod_code,certificate_group)],file.path(out,'cod_attribution_crosswalk_counts.csv'))
fwrite(x[event>0,.(n=.N),by=.(event,certificate_group)],file.path(out,'cod_audit_summary.csv'))
fwrite(x[time>60 & time<=96 & event>0,.(n=.N),by=.(event,cod_label,certificate_group)],file.path(out,'five_year_window_cod_counts.csv'))
fwrite(x[,.(patient_id,detailed,certificate,cod_code)],file.path(out,'local_patient_event_mapping.csv'))
N<-nrow(x);cache<-list()
add<-function(z,L,group,kind,scheme='Detailed'){
 if(!nrow(z))return(NULL)
 ev<-if(scheme=='Detailed')z$detailed else z$certificate
 list(id=z$id,bin=as.integer(pmin(round((z$time-L)*2),73)),event=ev,relative=z$time-L,n=nrow(z),L=L,group=group,kind=kind,scheme=scheme,K=if(scheme=='Detailed')6L else 5L)
}
for(L in c(0,12,36,60)){
 z<-if(L==0)copy(x) else copy(x[time>L]);z[,attained:=age+L/12];z[,agegroup:=fifelse(attained<65,'<65',fifelse(attained<75,'65-74','75+'))]
 cache[[length(cache)+1]]<-add(z,L,'Overall','Overall')
 cache[[length(cache)+1]]<-add(z,L,'Overall','Overall','Certificate')
 for(a in c('<65','65-74','75+')){
  cache[[length(cache)+1]]<-add(z[agegroup==a],L,a,'Attained age')
  cache[[length(cache)+1]]<-add(z[agegroup==a],L,a,'Attained age','Certificate')
 }
 if(L==60){
  for(v in c('subtype','site'))for(g in unique(z[[v]])){
   cache[[length(cache)+1]]<-add(z[get(v)==g],L,g,v)
   for(a in c('<65','65-74','75+'))cache[[length(cache)+1]]<-add(z[get(v)==g & agegroup==a],L,paste(g,a,sep=' | '),paste(v,'x age'))
  }
 }
}
calc<-function(cc,w){
 ww<-w[cc$id];count<-function(sel){ans<-numeric(73);if(any(sel)){v<-rowsum(ww[sel],cc$bin[sel],reorder=FALSE);ans[as.integer(rownames(v))]<-v[,1]};ans}
 removals<-count(rep(TRUE,length(ww)));risk<-rev(cumsum(rev(removals)))
 deaths<-sapply(seq_len(cc$K),function(k)count(cc$event==k))
 haz<-rowSums(deaths)/pmax(risk,1);before<-c(1,cumprod(1-haz)[1:72])
 cif<-colSums((deaths/pmax(risk,1)*before)[1:72,,drop=FALSE])
 if(cc$scheme=='Detailed')return(c(cif,sum(cif[2:5]),sum(cif[2:5])-cif[1],sum(cif),sum(cif[c(3,4)])))
 c(cif,rep(NA_real_,5))
}
point<-t(vapply(cache,calc,numeric(10),w=rep(1,N)))
# Independent AJ check for EVERY analyzed stratum and scheme, event-before-censor ties.
maxerr<-0
for(i in seq_along(cache)){
 cc<-cache[[i]];fit<-survfit(Surv(cc$relative,factor(cc$event,levels=0:cc$K))~1)
 v<-summary(fit,times=36,extend=TRUE)$pstate[1,-1]
 maxerr<-max(maxerr,max(abs(v-point[i,seq_len(cc$K)])))
}
stopifnot(maxerr<1e-10)
original<-fread(file.path(p,'primary_conditional_risks.csv'))
for(i in which(vapply(cache,function(cc)cc$scheme=='Detailed'&&cc$kind=='Overall',logical(1)))){
 L<-cache[[i]]$L/12;old<-original[kind=='Overall' & landmark==L & horizon==3]
 stopifnot(abs(point[i,1]-old[outcome=='Cancer',estimate])<1e-10,abs(point[i,7]-old[outcome=='Other',estimate])<1e-10)
}
B<-1000L;set.seed(20261008);boot<-array(NA_real_,c(length(cache),10,B))
for(b in 1:B){w<-tabulate(sample.int(N,N,replace=TRUE),nbins=N);boot[,,b]<-t(vapply(cache,calc,numeric(10),w=w));if(b%%100==0)cat('Survivorship bootstrap',b,'/',B,'\n')}
result<-rbindlist(lapply(seq_along(cache),function(i){
 cc<-cache[[i]];nm<-if(cc$scheme=='Detailed')c(names_d,'All other attributed','Other minus attributed cancer','All deaths','Specified nonmalignant total') else names_c
 rbindlist(lapply(seq_along(nm),function(j){
  sels<-if(j<=cc$K)j else switch(as.character(j),'7'=2:5,'8'=integer(),'9'=1:6,'10'=c(3,4))
  d<-sum(cc$event %in% sels & cc$relative<=36);support<-sum(cc$relative>36)
  data.table(index=i,scheme=cc$scheme,kind=cc$kind,group=cc$group,landmark=cc$L/12,n=cc$n,outcome=nm[j],events=d,observed_beyond_horizon=support,sparse=cc$n<100 || (j!=8 && d<20) || support<20,estimate=point[i,j],lower=quantile(boot[i,j,],.025),upper=quantile(boot[i,j,],.975),bootstrap=B)
 }))
}))
fwrite(result,file.path(out,'conditional_cause_age_risks.csv'));saveRDS(boot,file.path(out,'bootstrap.rds'))
# Age-standardize subtype/site contrasts among five-year survivors using a common 3-band distribution.
ref<-x[time>60];ref[,ag:=fifelse(age+5<65,'<65',fifelse(age+5<75,'65-74','75+'))]
weights<-ref[,.N,by=ag];weights[,weight:=N/sum(N)];fwrite(weights,file.path(out,'five_year_age_reference.csv'))
standard<-list()
for(v in c('subtype','site'))for(g in unique(ref[[v]])){
 ix<-vapply(weights$ag,function(a)which(vapply(cache,function(cc)cc$kind==paste(v,'x age')&&cc$group==paste(g,a,sep=' | '),logical(1))),integer(1))
 for(j in c(1,3,7,8)){
  draws<-colSums(boot[ix,j,]*weights$weight)
  standard[[length(standard)+1]]<-data.table(kind=v,group=g,outcome=c(names_d,'All other attributed','Other minus attributed cancer')[j],estimate=sum(point[ix,j]*weights$weight),lower=quantile(draws,.025),upper=quantile(draws,.975),reference='Observed five-year survivors; fixed empirical age-band weights',bootstrap=B)
 }
}
fwrite(rbindlist(standard),file.path(out,'age_standardized_five_year.csv'))
writeLines(c(paste('Independent AJ maximum discrepancy:',maxerr),'Primary cancer and other-cause overall estimates reproduced at all four landmarks.','1000 shared patient bootstrap samples. Fixed weights in age standardization; intervals conditional on chosen reference weights.',capture.output(sessionInfo())),file.path(out,'verification.txt'))
print(result[scheme=='Detailed' & landmark==5 & kind %chin% c('Overall','Attained age') & outcome %chin% c('SEER-attributed cancer','Cardiovascular','All other attributed','Malignant-coded other death')]);cat('COMPLETE\n')
