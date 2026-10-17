# Independent half-month-grid integration. All primary age/calendar discontinuities
# lie on this grid; numerical quadrature is exact for these piecewise rates.
suppressPackageStartupMessages(library(data.table))
p<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study");o<-file.path(p,'revision/smr')
x<-fread(file.path(p,'primary_analysis_cohort.csv'),colClasses=c(patient_id='character'))[time>60]
m<-fread(file.path(p,'survivorship/local_patient_event_mapping.csv'),colClasses=c(patient_id='character'))
x<-merge(x,m,by='patient_id');r<-fread(file.path(p,'revision/wonder_20261016/cdc_cardiovascular_reference_2009_2023.csv'))
r[,age0:=as.integer(sub('[-+].*','',age_group))]
rates<-setNames(r$rate_per_person_year,paste(r$year,r$sex,r$age0,sep='|'))
ev<-numeric(nrow(x));py<-numeric(nrow(x))
for(i in seq_len(nrow(x))){
 a<-x[i];end<-min(a$time/12,8);py[i]<-end-5
 mids<-seq(5+1/48,end-1/48,by=1/24)
 dx<-(0:11+.5)/12;dx<-dx[a$year+dx+a$time/12<=2024+1e-9]
 ages<-pmin(85,5*floor((a$age+.5+mids)/5))
 ev[i]<-mean(vapply(dx,function(d){
  yy<-floor(a$year+d+mids)
  rr<-rates[paste(yy,a$sex,ages,sep='|')]
  stopifnot(!anyNA(rr));sum(rr)/24
 },numeric(1)))
}
z<-fread(file.path(o,'local_person_time_expected.csv'),colClasses=c(patient_id='character'))
v<-merge(data.table(patient_id=x$patient_id,r_expected=ev,r_py=py),z,by='patient_id')
stopifnot(nrow(v)==5980,max(abs(v$r_expected-v$expected))<1e-10,max(abs(v$r_py-v$person_years))<1e-10)
observed<-x[time<=96 & cod_code %in% c(50060,50070,50080,50090,50100,50110),.N]
pt<-poisson.test(observed,T=sum(ev))
main<-fread(file.path(o,'primary_and_subgroups.csv'))[group=='Overall']
stopifnot(observed==91,abs(main$smr-unname(pt$estimate))<1e-12,
 max(abs(c(main$smr_lower,main$smr_upper)-pt$conf.int))<1e-10)
writeLines(c('PASS: independent R half-month-grid expected-death integration agrees for all 5980 people.',
 paste('Maximum person-level expected-count discrepancy:',max(abs(v$r_expected-v$expected))),
 'PASS: person-years and 91 observed cardiovascular deaths independently reproduced.',
 'PASS: exact Poisson interval agrees with stats::poisson.test.',capture.output(sessionInfo())),file.path(o,'independent_validation.txt'))
cat('VALIDATED N=',nrow(v),'O=',observed,'E=',sum(ev),'PY=',sum(py),'\n')
