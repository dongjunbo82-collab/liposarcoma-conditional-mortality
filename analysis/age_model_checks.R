# Reuse the same cohort, design matrices and competing-risk integration.
expr<-parse('analysis/age_definition_sensitivity.R')
for(e in expr){if(is.call(e)&&identical(e[[1]],as.name('set.seed')))break;eval(e)}
set.seed(20261012);B<-500L
pointchecks<-list(Linear=fitpred(rep(1,N),TRUE),Exclude_topcoded=fitpred(as.numeric(z$age<90)))
boots<-lapply(1:2,function(i)array(NA_real_,c(length(agegrid),4,B)))
for(b in seq_len(B)){
 w<-tabulate(sample.int(N,N,TRUE),N)
 boots[[1]][,,b]<-fitpred(w,TRUE);boots[[2]][,,b]<-fitpred(w*as.numeric(z$age<90))
 if(b%%100==0)cat('Check bootstrap',b,'\n')
}
rr<-rbindlist(lapply(1:2,function(v){a<-pointchecks[[v]];lo<-apply(boots[[v]],c(1,2),quantile,.025);hi<-apply(boots[[v]],c(1,2),quantile,.975);rbindlist(lapply(1:4,function(k)data.table(variant=names(pointchecks)[v],age=agegrid,outcome=colnames(a)[k],estimate=a[,k],lower=lo[,k],upper=hi[,k],valid_bootstrap=B)))}))
fwrite(rr,file.path(out,'continuous_age_checks.csv'))
# Apparent aggregate calibration against the nonparametric landmark estimate.
XX<-lapply(1:2,function(j)model.matrix(form,data.frame(age5=z$age5,interval=factor(j,levels=1:2))))
bet<-lapply(1:3,function(k)glm.fit(X,as.numeric(parts$e==k),offset=log(parts$exposure),family=poisson())$coefficients)
SS<-rep(1,N);FF<-matrix(0,N,4)
for(j in 1:2){a<-parts$interval==j;u<-sum(parts$e[a]==4)/sum(parts$exposure[a]);H<-cbind(exp(XX[[j]]%*%bet[[1]]),exp(XX[[j]]%*%bet[[2]]),exp(XX[[j]]%*%bet[[3]]),u);h<-rowSums(H);dt<-c(12,24)[j];FF<-FF+H/h*(SS*(-expm1(-h*dt)));SS<-SS*exp(-h*dt)
 for(i in c(1,100,1000)){num<-integrate(function(t)H[i,1]*exp(-h[i]*t),0,dt)$value;stopifnot(abs(num-H[i,1]/h[i]*(-expm1(-h[i]*dt)))<1e-10)}
}
sf<-survfit(Surv(tm,factor(ev,levels=0:4))~1,data=z);aj<-summary(sf,times=36,extend=TRUE)$pstate[1,-1];mod<-colMeans(FF)
fwrite(data.table(outcome=c('Cancer','Cardiovascular','Remaining other','Unknown'),AJ=aj,model=mod,difference_pp=100*(mod-aj)),file.path(out,'continuous_age_calibration.csv'))
writeLines(c('Two sensitivity checks use 500 shared patient-bootstrap resamples, seed 20261012.','Linear age replaces the spline. Top-coding sensitivity excludes baseline ages coded 90+. All 500 fits succeeded.','Numerical integration agrees with the analytic cumulative-incidence increment within 1e-10.','Apparent aggregate calibration compares model averages with landmark Aalen-Johansen; it is not external prediction validation.'),file.path(out,'continuous_age_checks_methods.txt'))
print(rr[age%in%c(65,85)]);print(data.table(outcome=1:4,AJ=aj,model=mod))
