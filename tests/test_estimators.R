# Independent small artificial example; contains no study records.
suppressPackageStartupMessages({library(data.table);library(survival)})
src<-readLines('analysis/population_primary.R',warn=FALSE)
i<-grep('^prep<-function',src);j<-grep('^groups<-list',src)
eval(parse(text=paste(src[i:(j-1)],collapse='\n')))
z<-data.table(id=1:5,time=c(6,12,18,24,48),event=c(1L,0L,2L,3L,0L))
c<-prep(z,0,'Artificial','Overall');v<-calc(c,rep(1,5))
stopifnot(max(abs(v[1:4]-c(.2,4/15,4/15,1/15)))<1e-12)
f<-survfit(Surv(time,factor(event,levels=0:3))~1,data=z)
a<-summary(f,times=36,extend=TRUE)$pstate[1,2:4]
stopifnot(max(abs(a-v[1:3]))<1e-12)
# Conditioning must exclude deaths/censors exactly at a positive landmark.
c2<-prep(z,12,'Artificial','Overall');stopifnot(c2$n==3)
v2<-calc(c2,rep(1,5));stopifnot(max(abs(v2[1:3]-c(0,1/3,1/3)))<1e-12)
# Preserve event-before-censor handling when times are tied.
z2<-data.table(id=1:4,time=c(12,12,24,48),event=c(1L,0L,2L,0L))
v3<-calc(prep(z2,0,'Artificial','Overall'),rep(1,4))
stopifnot(max(abs(v3[1:2]-c(.25,.375)))<1e-12)
cat('PASS: analytic probabilities, independent survfit check, landmark eligibility and tied times.\n')
