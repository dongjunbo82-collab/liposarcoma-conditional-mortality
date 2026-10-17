suppressPackageStartupMessages({library(data.table);library(survival)})
p <- Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study"); out <- file.path(p,'revision')
x <- fread(file.path(p,'primary_analysis_cohort.csv'),colClasses=c(patient_id='character'))
x[,id:=.I]
r <- fread(Sys.getenv("SEER_EXPORT", unset="inputs/seer_export.tsv"),colClasses='character');setnames(r,paste0('v',1:23))
r[,`:=`(rank=fifelse(v10=='One primary only',1L,suppressWarnings(as.integer(sub('^(\\d+).*','\\1',v10)))),yr=as.integer(v2),ord=.I)]
setorder(r,v1,rank,yr,ord);r<-r[!duplicated(v1)]
x<-merge(x,r[,.(patient_id=v1,raw_year=yr,raw_sequence=v10,grade_old=v11,grade_new=v12,radiation=v17,chemotherapy=v18)],by='patient_id',sort=FALSE)
stopifnot(nrow(x)==12219,all(x$year==x$raw_year),all(x$sequence==x$raw_sequence))
m<-fread(file.path(p,'survivorship/local_patient_event_mapping.csv'),colClasses=c(patient_id='character'))
x<-merge(x,m,by='patient_id',sort=FALSE);setorder(x,id)
stopifnot(identical(x$id,1:12219),sum(x$time>60)==5980)
x[,era:=fifelse(year<=2017,'2004-2017','2018-2023')]
x[,grade_raw:=fifelse(year<=2017,grade_old,grade_new)]
x[,grade_missing:=grepl('Unknown|Blank|cannot be assessed',grade_raw)]
x[,grade_group:=fcase(grade_old=='Well differentiated; Grade I','I',grade_old=='Moderately differentiated; Grade II','II',grade_old=='Poorly differentiated; Grade III','III',grade_old=='Undifferentiated; anaplastic; Grade IV','IV',default='Unknown')]
missing<-rbindlist(lapply(c('Overall','subtype','site'),function(v){z<-copy(x);z[,group:=if(v=='Overall')'Overall' else get(v)];z[,dimension:=v];z[,.(n=.N,missing=sum(grade_missing),missing_pct=100*mean(grade_missing),five_year_n=sum(time>60),five_year_missing=sum(grade_missing & time>60)),by=.(era,dimension,group)]}))
fwrite(missing,file.path(out,'grade_missingness.csv'))
fwrite(x[,.(n=.N,five_year_n=sum(time>60)),by=.(era,grade_raw,grade_missing)],file.path(out,'grade_distributions.csv'))
fwrite(x[,.(n=.N,five_year_n=sum(time>60)),by=.(era,subtype,site,grade_raw,grade_missing)],file.path(out,'grade_distribution_by_histology_site.csv'))
rt_yes<-c('Beam radiation','Radiation, NOS  method or source not specified','Combination of beam with implants or isotopes','Radioactive implants (includes brachytherapy) (1988+)','Radioisotopes (1988+)')
x[,rt:=fifelse(radiation%chin%rt_yes,'Recorded radiation',radiation)]
fwrite(rbindlist(list(x[,.(n=.N,five_year_n=sum(time>60)),by=.(category=radiation)][,variable:='Radiation'],x[,.(n=.N,five_year_n=sum(time>60)),by=.(category=chemotherapy)][,variable:='Chemotherapy'])),file.path(out,'treatment_field_distributions.csv'))
cache<-list();add<-function(z,L,kind,g){z<-z[time>L];if(nrow(z))cache[[length(cache)+1]]<<-list(id=z$id,t=z$time-L,event=z$detailed,n=nrow(z),L=L,kind=kind,group=g,bin=as.integer(pmin(round((z$time-L)*2),73)))}
for(g in c('I','II','III','IV','Unknown'))add(x[year<=2017 & grade_group==g],60,'Grade 2004-2017',g)
for(g in c('I','II','III','IV','Unknown'))add(x[year<=2017 & subtype=='Dedifferentiated' & grade_group==g],60,'Dedifferentiated grade 2004-2017',g)
for(g in unique(x$rt))add(x[rt==g],60,'Radiation',g)
for(g in unique(x$chemotherapy))add(x[chemotherapy==g],60,'Chemotherapy',g)
for(L in c(0,12,36,60)){
 add(x[subtype=='Myxoid/round cell'],L,'Myxoid definition','8852 and 8853')
 add(x[histology=='8852'],L,'Myxoid definition','8852 only')
}
calc<-function(cc,w){ww<-w[cc$id];cnt<-function(s){a<-numeric(73);if(any(s)){q<-rowsum(ww[s],cc$bin[s],reorder=FALSE);a[as.integer(rownames(q))]<-q[,1]};a};risk<-rev(cumsum(rev(cnt(rep(TRUE,cc$n)))));deaths<-sapply(1:6,function(k)cnt(cc$event==k));h<-rowSums(deaths)/pmax(risk,1);before<-c(1,cumprod(1-h)[1:72]);v<-colSums((deaths/pmax(risk,1)*before)[1:72,,drop=FALSE]);c(v[1],sum(v[2:5]),v[3])}
N<-nrow(x);point<-t(vapply(cache,calc,numeric(3),w=rep(1,N)));err<-0
for(i in seq_along(cache)){cc<-cache[[i]];sf<-survfit(Surv(cc$t,factor(cc$event,levels=0:6))~1);v<-summary(sf,times=36,extend=TRUE)$pstate[1,-1];err<-max(err,abs(point[i,]-c(v[1],sum(v[2:5]),v[3])))}
stopifnot(err<1e-10)
set.seed(20261015);B<-1000L;boot<-array(NA_real_,c(length(cache),3,B))
for(b in 1:B){w<-tabulate(sample.int(N,N,replace=TRUE),nbins=N);boot[,,b]<-t(vapply(cache,function(cc)if(sum(w[cc$id])>0)calc(cc,w)else rep(NA_real_,3),numeric(3)));if(b%%200==0)cat('Clinical context bootstrap',b,'\n')}
res<-rbindlist(lapply(seq_along(cache),function(i){cc<-cache[[i]];ev<-cc$event[cc$t<=36];events<-c(sum(ev==1),sum(ev%in%2:5),sum(ev==3));rbindlist(lapply(1:3,function(j)data.table(kind=cc$kind,group=cc$group,landmark=cc$L/12,n=cc$n,outcome=c('Cancer','Other','Cardiovascular')[j],events=events[j],all_deaths=sum(ev>0),beyond=sum(cc$t>36),sparse=cc$n<100||events[j]<20||sum(cc$t>36)<20,estimate=point[i,j],lower=quantile(boot[i,j,],.025,na.rm=TRUE),upper=quantile(boot[i,j,],.975,na.rm=TRUE),valid_bootstrap=sum(is.finite(boot[i,j,])))))}))
fwrite(res,file.path(out,'clinical_context_risks.csv'))
fwrite(res[kind%chin%c('Radiation','Chemotherapy') & outcome=='Cardiovascular',.(kind,group,n,cvd_deaths=events,all_deaths,cvd_share_of_deaths_pct=100*events/all_deaths,estimate,lower,upper,sparse)],file.path(out,'treatment_cardiovascular_composition.csv'))
delta<-rbindlist(lapply(c(0,1,3,5),function(L){ii<-which(vapply(cache,function(cc)cc$kind=='Myxoid definition'&&cc$L==L*12,logical(1)));v<-boot[ii[2],,]-boot[ii[1],,];data.table(landmark=L,outcome=c('Cancer','Other','Cardiovascular'),difference=point[ii[2],]-point[ii[1],],lower=apply(v,1,quantile,.025),upper=apply(v,1,quantile,.975))}))
fwrite(delta,file.path(out,'myxoid_definition_contrasts.csv'))
writeLines(c(paste('Maximum absolute AJ discrepancy against survival::survfit:',err),'1000 shared patient bootstrap samples, seed 20261015. Pointwise descriptive intervals. Zero-event intervals uninformative.', 'Grade fields analyzed separately by era; historic grade I-IV not equated to FNCLCC. Grade analyses restricted to 2004-2017; unknown grade retained.', 'Radiation received uses explicitly administered categories; none/unknown, refused and recommended/unknown remain distinct. No/unknown chemotherapy is not an untreated control.',capture.output(sessionInfo())),file.path(out,'clinical_context_verification.txt'))
print(missing[dimension=='Overall']);print(res[landmark==5 & (kind=='Myxoid definition'|kind=='Chemotherapy'|kind=='Radiation')]);cat('COMPLETE\n')

