suppressPackageStartupMessages({library(data.table);library(survival)})
root <- '.'
p <- Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study")
x<-fread(file.path(p,'primary_analysis_cohort.csv'),colClasses=c(patient_id='character'))
x[,`:=`(follow_months=months,event_code=event,main_subtype=subtype,main_site=site)]
raw<-fread(Sys.getenv("SEER_EXPORT", unset="inputs/seer_export.tsv"),colClasses='character')
setnames(raw,names(raw),paste0('v',seq_along(raw)))
raw[,seqrank:=fifelse(v10=='One primary only',1L,suppressWarnings(as.integer(sub('^(\\d+).*','\\1',v10))))]
raw[,yr:=as.integer(v2)];raw[,ord:=.I];setorder(raw,v1,seqrank,yr,ord)
r<-raw[!duplicated(v1),.(patient_id=v1,stage_raw=v14,seqrank)]
x<-merge(x,r,by='patient_id',sort=FALSE)
main_ids<-fread(file.path(Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study"),"primary_analysis_cohort.csv"),colClasses=c(patient_id='character'))$patient_id
x<-x[patient_id %chin% main_ids]
x[,stage:=factor(fcase(stage_raw=='Localized only','Localized',grepl('^Regional',stage_raw),'Regional',grepl('^Distant',stage_raw),'Distant',default='Unknown'),levels=c('Localized','Regional','Distant','Unknown'))]
x[,agegroup:=cut(age,c(-Inf,49,64,74,Inf),labels=c('<50','50-64','65-74','75+'))]
x[,sex:=factor(sex)];x[,st:=factor(main_subtype,levels=c('Well-differentiated','Dedifferentiated','Myxoid/round cell','Pleomorphic','Other/NOS'))]
x[,site:=factor(main_site,levels=c('Extremity','Retroperitoneum','Trunk/pelvis soft tissue','Other'))]
x[,era:=factor(fcase(year<=2009,'2004-2009',year<=2015,'2010-2015',default='2016-2023'))]
x[,id:=.I]
# Fixed reference: five-year observed survivors diagnosed by 2015; supported at all landmarks.
ref<-x[follow_months>60 & year<=2015]
stopifnot(nrow(ref)>0,!anyNA(x$agegroup),!anyDuplicated(x$patient_id))
out<-file.path(Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study"),"standardized");dir.create(out,showWarnings=FALSE)
form<-~ interval*(agegroup+st+site+stage)+sex+era
form0<-~ interval+agegroup+st+site+stage+sex+era
make_parts<-function(L){
 z<-if(L==0)copy(x)else x[follow_months>L]
 z[,tm:=pmin(pmax(follow_months,0.5)-L,36)]
 z[,ev:=ifelse(pmax(follow_months,0.5)-L>36,0L,event_code)]
 parts<-rbindlist(lapply(1:2,function(j){v<-copy(z);lo<-c(0,12)[j];hi<-c(12,36)[j];v<-v[tm>lo];v[,exposure:=pmin(tm,hi)-lo];v[,interval:=factor(j,levels=1:2)];v[,e:=ifelse(tm<=hi,ev,0L)];v}))
 P<-lapply(1:2,function(j){v<-copy(ref);v[,interval:=factor(j,levels=1:2)];model.matrix(form,v)})
 list(z=z,parts=parts,X=model.matrix(form,parts),X0=model.matrix(form0,parts),P=P)
}
cache<-lapply(c(0,12,36,60),make_parts)
fit_one<-function(c,w,ph=FALSE){
 X<-if(ph)c$X0 else c$X
 bb<-list();warn<-character();dev<-numeric()
 for(k in 1:2){
  f<-withCallingHandlers(glm.fit(X,as.numeric(c$parts$e==k),family=poisson(),offset=log(c$parts$exposure),weights=w[c$parts$id],control=glm.control(maxit=60)),warning=function(e){warn<<-c(warn,conditionMessage(e));invokeRestart('muffleWarning')})
  if(!f$converged || anyNA(f$coefficients) || any(!is.finite(f$coefficients)))stop('nonconvergent or aliased model')
  bb[[k]]<-f$coefficients;dev[k]<-f$deviance
 }
 unknown<-vapply(1:2,function(j){ii<-c$parts$interval==j;sum(w[c$parts$id[ii]]*(c$parts$e[ii]==3))/sum(w[c$parts$id[ii]]*c$parts$exposure[ii])},numeric(1))
 list(beta=bb,unknown=unknown,dev=dev,warnings=warn)
}
predict_one<-function(c,f,w){
 S<-rep(1,nrow(ref));F<-matrix(0,nrow(ref),3)
 for(j in 1:2){
  H<-cbind(exp(c$P[[j]]%*%f$beta[[1]]),exp(c$P[[j]]%*%f$beta[[2]]),rep(f$unknown[j],nrow(ref)))
  h<-rowSums(H);inc<-S*(-expm1(-h*c(12,24)[j]));F<-F+H/h*inc;S<-S*exp(-h*c(12,24)[j])
 }
 stopifnot(all(is.finite(F)),max(abs(rowSums(F)+S-1))<1e-10)
 wref<-w[ref$id];p<-colSums(F*wref)/sum(wref)
 c(cancer=p[1],other=p[2],unknown=p[3],difference=p[2]-p[1])
}
N<-nrow(x);w<-rep(1,N);point<-list();diags<-list()
for(i in 1:4){
 f<-fit_one(cache[[i]],w);p0<-fit_one(cache[[i]],w,TRUE)
 point[[i]]<-data.table(landmark=c(0,1,3,5)[i],t(predict_one(cache[[i]],f,w)))
 for(k in 1:2){d<-p0$dev[k]-f$dev[k];df<-ncol(cache[[i]]$X)-ncol(cache[[i]]$X0);diags[[length(diags)+1]]<-data.table(landmark=c(0,1,3,5)[i],cause=k,n=nrow(cache[[i]]$z),events=sum(cache[[i]]$z$ev==k),parameters=ncol(cache[[i]]$X),LR=d,df=df,p=pchisq(d,df,lower.tail=FALSE),max_abs_beta=max(abs(f$beta[[k]])),warnings=paste(f$warnings,collapse=';'))}
 saveRDS(f,file.path(out,paste0('piecewise_model_',i,'.rds')))
}
fwrite(rbindlist(point),file.path(out,'standardized_point_estimates.csv'))
fwrite(rbindlist(diags),file.path(out,'time_varying_effect_diagnostics.csv'))
fwrite(ref[,.(n=.N),by=.(agegroup,st,site,stage,era)],file.path(out,'reference_distribution.csv'))
# Apparent calibration: aggregate individual model predictions over each observed landmark cohort.
calibration<-list()
for(i in 1:4){
 cc<-cache[[i]]; fit<-fit_one(cc,rep(1,N)); z<-cc$z
 S<-rep(1,nrow(z)); F<-matrix(0,nrow(z),3)
 for(j in 1:2){
  z[,interval:=factor(j,levels=1:2)]; XX<-model.matrix(form,z)
  H<-cbind(exp(XX%*%fit$beta[[1]]),exp(XX%*%fit$beta[[2]]),rep(fit$unknown[j],nrow(z)))
  h<-rowSums(H); F<-F+H/h*(S*(-expm1(-h*c(12,24)[j])));S<-S*exp(-h*c(12,24)[j])
 }
 sf<-survfit(Surv(tm,factor(ev,levels=0:3))~1,data=z)
 su<-summary(sf,times=36,extend=TRUE)
 for(k in 1:2)calibration[[length(calibration)+1]]<-data.table(landmark=c(0,1,3,5)[i],cause=k,AJ=su$pstate[1,k+1],model=mean(F[,k]))
}
fwrite(rbindlist(calibration),file.path(out,'apparent_aggregate_calibration.csv'))
args<-commandArgs(TRUE);B<-if(length(args))as.integer(args[1])else 200L
set.seed(20261007);boots<-list();fails<-list()
for(b in seq_len(B)){
 w<-tabulate(sample.int(N,N,replace=TRUE),nbins=N)
 result<-tryCatch(rbindlist(lapply(1:4,function(i){f<-fit_one(cache[[i]],w);data.table(replicate=b,landmark=c(0,1,3,5)[i],t(predict_one(cache[[i]],f,w)))})),error=function(e)e)
 if(inherits(result,'error'))fails[[length(fails)+1]]<-data.table(replicate=b,error=conditionMessage(result))else boots[[length(boots)+1]]<-result
 if(b%%10==0||b==1)cat('Bootstrap',b,'/',B,'complete\n')
}
br<-rbindlist(boots);fwrite(br,file.path(out,'patient_bootstrap_draws.csv'))
fwrite(if(length(fails))rbindlist(fails,fill=TRUE)else data.table(replicate=integer(),error=character()),file.path(out,'bootstrap_failures.csv'))
long<-melt(br,id.vars=c('replicate','landmark'),variable.name='outcome')
ci<-long[,.(lower=quantile(value,.025),upper=quantile(value,.975),successful_replicates=.N),by=.(landmark,outcome)]
po<-melt(rbindlist(point),id.vars='landmark',variable.name='outcome',value.name='estimate')
fwrite(merge(po,ci,by=c('landmark','outcome')),file.path(out,'standardized_risks_with_CI.csv'))
# Paired contrast across landmarks preserves same-patient bootstrap dependence.
delta<-br[landmark==5,.(replicate,end=difference)][br[landmark==0,.(replicate,start=difference)],on='replicate'][,end-start]
fwrite(data.table(contrast='Change in other-minus-cancer risk: landmark 5 vs 0',estimate=point[[4]]$difference-point[[1]]$difference,lower=quantile(delta,.025),upper=quantile(delta,.975)),file.path(out,'paired_landmark_contrast.csv'))
writeLines(c(paste('Model N',N),paste('Reference N',nrow(ref)),paste('Bootstrap attempts',B),paste('Successful',length(boots)),capture.output(sessionInfo())),file.path(out,'run_metadata.txt'))
print(rbindlist(point));print(rbindlist(diags))
