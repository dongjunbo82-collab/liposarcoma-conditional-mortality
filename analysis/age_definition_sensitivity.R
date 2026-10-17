suppressPackageStartupMessages({library(data.table);library(survival);library(splines)})
P<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study")
out<-file.path(P,'revision');x<-fread(file.path(P,'primary_analysis_cohort.csv'),colClasses=c(patient_id='character'))
m<-fread(file.path(P,'survivorship/local_patient_event_mapping.csv'),colClasses=c(patient_id='character'))
stopifnot(!anyDuplicated(x$patient_id),identical(x$patient_id,m$patient_id))
x[,detailed:=m$detailed]
# Describe calendar composition, not incidence or temporal treatment effects.
x[,era:=fcase(year<=2009,'2004-2009',year<=2015,'2010-2015',default='2016-2023')]
comp<-rbindlist(lapply(c(0,1,3,5),function(L){z<-x[time>L*12];z[,.(landmark=L,n=.N,percent=100*.N/nrow(z)),by=era]}))
fwrite(comp,file.path(out,'landmark_diagnosis_composition.csv'))
z<-x[time>60];z[,id:=.I];z[,age5:=age+5];z[,tm:=pmin(time-60,36)]
z[,ev:=fcase(time-60>36,0L,event==1,1L,detailed==3,2L,event==2,3L,event==3,4L,default=0L)]
stopifnot(nrow(z)==5980,all(z$tm>0))
knots<-as.numeric(quantile(z$age5,c(1/3,2/3)));bounds<-range(z$age5)
parts<-rbindlist(lapply(1:2,function(j){lo<-c(0,12)[j];hi<-c(12,36)[j];a<-copy(z[tm>lo]);a[,exposure:=pmin(tm,hi)-lo];a[,interval:=factor(j,levels=1:2)];a[,e:=ifelse(tm<=hi,ev,0L)];a}))
form<-~interval+ns(age5,knots=knots,Boundary.knots=bounds)
X<-model.matrix(form,parts);agegrid<-35:90
predX<-lapply(1:2,function(j)model.matrix(form,data.frame(age5=agegrid,interval=factor(j,levels=1:2))))
fitpred<-function(w,linear=FALSE){
 xx<-if(linear)model.matrix(~interval+age5,parts) else X
 pp<-if(linear)lapply(1:2,function(j)model.matrix(~interval+age5,data.frame(age5=agegrid,interval=factor(j,levels=1:2)))) else predX
 bet<-list()
 for(k in 1:3){f<-suppressWarnings(glm.fit(xx,as.numeric(parts$e==k),weights=w[parts$id],offset=log(parts$exposure),family=poisson(),control=glm.control(maxit=80)));if(!f$converged||any(!is.finite(f$coefficients)))stop('Fit failed');bet[[k]]<-f$coefficients}
 unk<-vapply(1:2,function(j){a<-parts$interval==j;sum(w[parts$id[a]]*(parts$e[a]==4))/sum(w[parts$id[a]]*parts$exposure[a])},numeric(1))
 S<-rep(1,length(agegrid));F<-matrix(0,length(agegrid),4)
 for(j in 1:2){H<-cbind(exp(pp[[j]]%*%bet[[1]]),exp(pp[[j]]%*%bet[[2]]),exp(pp[[j]]%*%bet[[3]]),unk[j]);h<-rowSums(H);F<-F+H/h*(S*(-expm1(-h*c(12,24)[j])));S<-S*exp(-h*c(12,24)[j])}
 stopifnot(all(is.finite(F)),all(F>=0),max(abs(rowSums(F)+S-1))<1e-10)
 cbind(Cancer=F[,1],Other=F[,2]+F[,3],Cardiovascular=F[,2],Unknown=F[,4])
}
N<-nrow(z);point<-fitpred(rep(1,N));lin<-fitpred(rep(1,N),TRUE)
set.seed(20261011);B<-1000L;boot<-array(NA_real_,c(length(agegrid),4,B));fail<-integer()
for(b in seq_len(B)){w<-tabulate(sample.int(N,N,TRUE),N);p<-try(fitpred(w),silent=TRUE);if(inherits(p,'try-error'))fail<-c(fail,b)else boot[,,b]<-p;if(b%%100==0)cat('Bootstrap',b,'\n')}
stopifnot(length(fail)<10)
lo<-apply(boot,c(1,2),quantile,.025,na.rm=TRUE);hi<-apply(boot,c(1,2),quantile,.975,na.rm=TRUE)
rr<-rbindlist(lapply(1:4,function(k)data.table(age=agegrid,outcome=colnames(point)[k],estimate=point[,k],lower=lo[,k],upper=hi[,k],linear_age_estimate=lin[,k],valid_bootstrap=B-length(fail))))
# Do not assign exact ages to top-coded records without checking the consequence.
w<-as.numeric(z$age<90);notop<-fitpred(w)
rr[,without_topcoded:=as.vector(notop)]
fwrite(rr,file.path(out,'continuous_age_risks.csv'))
fwrite(z[,.(n=.N,cancer=sum(ev==1),cardiovascular=sum(ev==2),other=sum(ev%in%c(2,3)),unknown=sum(ev==4),beyond=sum(time>96)),by=.(ageband=cut(age5,c(-Inf,64,74,Inf),labels=c('<65','65-74','75+')))],file.path(out,'continuous_age_support.csv'))
# Independent numerical integration verifies the closed-form competing hazard calculation.
# A saturated interval-only Poisson fit checks the likelihood exposure/event bookkeeping.
for(k in 1:3){f<-glm(as.numeric(e==k)~0+interval+offset(log(exposure)),data=parts,family=poisson());for(j in 1:2){a<-parts$interval==j;stopifnot(abs(exp(coef(f)[j])-sum(parts$e[a]==k)/sum(parts$exposure[a]))<1e-8)}}
writeLines(c('Exploratory age-definition sensitivity; five-year survivors; next 36 months.','Natural cubic spline with two interior knots (three age degrees of freedom); interval-specific baseline hazard (0-12,12-36 months); no age-by-time interaction.','Cause-specific Poisson models: SEER cancer, cardiovascular other death, remaining other death. Unknown cause has interval-specific aggregate rate. Competing hazards integrated jointly.','Unadjusted age-associated population risks, not independent age effects or individual predictions.',paste('Interior knots',paste(knots,collapse=','),'boundary knots',paste(bounds,collapse=',')),paste('Five-year survivors',N,'top-coded baseline age records',sum(z$age>=90)),paste('Bootstrap successes',B-length(fail),'of',B,'seed 20261011'),'Primary calculation assigns baseline 90+ to 90; exclusion sensitivity removes these observations. Report ages 50,65,75,85; retain full 35-90 grid for audit.','Compared with linear age; no search for an optimal cutpoint, no confirmatory nonlinearity test.','Probability conservation and intercept-only Poisson event/exposure identities passed.',capture.output(sessionInfo())),file.path(out,'continuous_age_methods.txt'))
print(rr[age%in%c(50,65,75,85)]);print(comp)
