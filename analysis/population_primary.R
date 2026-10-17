suppressPackageStartupMessages({library(data.table);library(survival)})
out<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study");dir.create(out,showWarnings=FALSE,recursive=TRUE)
r<-fread(Sys.getenv("SEER_EXPORT", unset="inputs/seer_export.tsv"),colClasses='character')
setnames(r,names(r),paste0('v',seq_along(r)))
r[,age:=as.integer(sub('[^0-9].*','',v3))];r[,year:=as.integer(v2)]
r[,months:=suppressWarnings(as.integer(v19))]
r[,rank:=fifelse(v10=='One primary only',1L,suppressWarnings(as.integer(sub('^(\\d+).*','\\1',v10))))]
r[,ord:=.I];setorder(r,v1,rank,year,ord)
flow<-list(data.table(step='Exported records after SEER malignant/known-age filter',n=nrow(r)))
r<-r[!duplicated(v1)];flow[[2]]<-data.table(step='First-sequenced liposarcoma record per patient',n=nrow(r))
r<-r[year>=2004];flow[[3]]<-data.table(step='Diagnosis 2004-2023',n=nrow(r))
r<-r[age>=18];flow[[4]]<-data.table(step='Adults age 18+',n=nrow(r))
r<-r[v9 %chin% c('Positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies')]
flow[[5]]<-data.table(step='Histologically confirmed',n=nrow(r))
r[,first:=v10 %chin% c('One primary only','1st of 2 or more primaries')]
flow[[6]]<-data.table(step='First or only primary',n=r[first==TRUE,.N])
flow[[7]]<-data.table(step='Known survival months: main cohort',n=r[first==TRUE & !is.na(months),.N])
fwrite(rbindlist(flow),file.path(out,'cohort_flow.csv'))
r<-r[!is.na(months)]
r[,event:=fcase(grepl('attributable to this cancer',v21,fixed=TRUE),1L,grepl('causes other than this cancer',v22,fixed=TRUE),2L,grepl('missing/unknown COD',v21,fixed=TRUE),3L,default=0L)]
stopifnot(all((r$v20=='Dead')==(r$event>0)),all(r$v8=='Malignant'))
r[,subtype:=fcase(v7=='8851','Well-differentiated',v7=='8858','Dedifferentiated',v7 %chin% c('8852','8853'),'Myxoid/round cell',v7=='8854','Pleomorphic',default='Other/NOS')]
r[,sitecode:=sub('-.*','',v6)]
r[,site:=fcase(sitecode=='C48.0','Retroperitoneum',sitecode %chin% c('C49.1','C49.2'),'Extremity',sitecode %chin% c('C49.3','C49.4','C49.5','C49.6'),'Trunk/pelvis soft tissue',default='Other')]
r[,stage:=fcase(v14=='Localized only','Localized',grepl('^Regional',v14),'Regional',grepl('^Distant',v14),'Distant',default='Unknown')]
r[,time:=pmax(months,.5)];x<-copy(r[first==TRUE]);x[,id:=.I]
fwrite(x[,.(patient_id=v1,year,age,sex=v4,histology=v7,subtype,site,sitecode,stage,sequence=v10,months,time,event)],file.path(out,'primary_analysis_cohort.csv'))
# Table 1 and follow-up; no baseline significance tests.
base<-rbindlist(lapply(c('Overall',unique(x$subtype)),function(g){z<-if(g=='Overall')x else x[subtype==g];res<-list(data.table(group=g,variable='N',level='All',n=nrow(z),percent=100));for(v in c('site','stage','v4','event')){tab<-z[,.N,by=v];setnames(tab,c('level','n'));tab[,`:=`(group=g,variable=v,percent=100*n/nrow(z))];res[[length(res)+1]]<-tab};rbindlist(res,use.names=TRUE)}))
fwrite(base,file.path(out,'table1_counts.csv'))
fwrite(x[,.(n=.N,age_median=as.numeric(median(age)),age_q1=quantile(age,.25),age_q3=quantile(age,.75)),by=subtype],file.path(out,'table1_age.csv'))
rev<-survfit(Surv(time,event==0)~1,data=x)
writeLines(c(paste('Reverse Kaplan-Meier median follow-up (months):',summary(rev)$table['median']),paste('Main cohort N',nrow(x)),paste('Unknown-COD deaths',sum(x$event==3))),file.path(out,'cohort_metadata.txt'))
# Nonparametric cumulative incidence on a half-month grid; event-before-censor ties.
prep<-function(z,L,group,kind,scenario='Primary'){
 z<-if(L==0)z else z[time>L]
 if(!nrow(z))return(NULL)
 relative<-z$time-L
 bin<-as.integer(pmin(round(relative*2),121));stopifnot(all(bin>=1))
 list(id=z$id,bin=bin,event=z$event,n=nrow(z),relative=relative,L=L,group=group,kind=kind,scenario=scenario)
}
calc<-function(c,w){
 ww<-w[c$id]
 count<-function(sel){ans<-numeric(121);if(any(sel)){v<-rowsum(ww[sel],c$bin[sel],reorder=FALSE);ans[as.integer(rownames(v))]<-v[,1]};ans}
 removals<-count(rep(TRUE,length(ww)));risk<-rev(cumsum(rev(removals)))
 deaths<-sapply(1:3,function(k)count(c$event==k))
 hazard<-rowSums(deaths)/pmax(risk,1)
 before<-c(1,cumprod(1-hazard)[1:120])
 inc<-apply(deaths/pmax(risk,1)*before,2,cumsum)
 v<-inc[c(72,120),,drop=FALSE]
 c(v[1,1],v[1,2],v[1,3],v[1,2]-v[1,1],v[2,1],v[2,2],v[2,3],v[2,2]-v[2,1])
}
groups<-list(list(name='Overall',kind='Overall',z=x))
for(g in unique(x$subtype))groups[[length(groups)+1]]<-list(name=g,kind='Subtype',z=x[subtype==g])
for(g in unique(x$site))groups[[length(groups)+1]]<-list(name=g,kind='Site',z=x[site==g])
for(s in c('Well-differentiated','Dedifferentiated','Myxoid/round cell','Pleomorphic'))for(t in c('Retroperitoneum','Extremity','Trunk/pelvis soft tissue'))groups[[length(groups)+1]]<-list(name=paste(s,t,sep=' | '),kind='Subtype x site',z=x[subtype==s & site==t])
cache<-list()
for(g in groups)for(L in c(0,12,36,60)){cc<-prep(g$z,L,g$name,g$kind);if(!is.null(cc))cache[[length(cache)+1]]<-cc}
N<-nrow(x);w<-rep(1,N);point<-t(vapply(cache,calc,numeric(8),w=w))
# Independent validation against survival::survfit for all four overall landmarks.
for(i in which(vapply(cache,function(c)c$kind=='Overall',logical(1)))){
 cc<-cache[[i]];z<-x[cc$id];fit<-survfit(Surv(time-cc$L,factor(event,levels=0:3))~1,data=z)
 ss<-summary(fit,times=c(36,60),extend=TRUE)$pstate[,2:4]
 stopifnot(max(abs(as.numeric(t(ss))-point[i,c(1,2,3,5,6,7)]))<1e-10)
}
# Future-window support flags; retain estimates but flag sparse/unsupported inference.
metadata<-rbindlist(lapply(seq_along(cache),function(i){cc<-cache[[i]];rbindlist(lapply(c(36,60),function(H){ds<-cc$event[cc$relative<=H];data.table(index=i,kind=cc$kind,group=cc$group,landmark=cc$L/12,horizon=H/12,n=cc$n,cancer_events=sum(ds==1),other_events=sum(ds==2),unknown_events=sum(ds==3),observed_beyond_horizon=sum(cc$relative>H),sparse=cc$n<100||sum(ds==1)<20||sum(ds==2)<20||sum(cc$relative>H)<20)}))}))
fwrite(metadata,file.path(out,'landmark_support.csv'))
B<-1000L;set.seed(20261008)
boot<-array(NA_real_,c(length(cache),8,B))
for(b in seq_len(B)){
 w<-tabulate(sample.int(N,N,replace=TRUE),nbins=N)
 boot[,,b]<-t(vapply(cache,calc,numeric(8),w=w))
 if(b%%100==0)cat('Primary AJ bootstrap',b,'of',B,'\n')
}
result<-rbindlist(lapply(seq_along(cache),function(i){rbindlist(lapply(1:8,function(j){meta<-metadata[index==i & horizon==if(j<=4)3 else 5];cbind(meta,data.table(outcome=c('Cancer','Other','Unknown','Difference')[(j-1)%%4+1],estimate=point[i,j],lower=quantile(boot[i,j,],.025),upper=quantile(boot[i,j,],.975),bootstrap=B))}))}))
fwrite(result,file.path(out,'primary_conditional_risks.csv'));saveRDS(boot,file.path(out,'primary_bootstrap_array.rds'))
# Same-patient bootstrap contrasts across landmarks for each group.
contrasts<-rbindlist(lapply(seq_along(groups),function(g){ix<-which(vapply(cache,function(c)c$kind==groups[[g]]$kind&&c$group==groups[[g]]$name,logical(1)));if(length(ix)!=4)return(NULL);val<-boot[ix[4],4,]-boot[ix[1],4,];data.table(kind=groups[[g]]$kind,group=groups[[g]]$name,contrast='Change in risk difference: 5-year vs diagnosis',estimate=point[ix[4],4]-point[ix[1],4],lower=quantile(val,.025),upper=quantile(val,.975))}))
fwrite(contrasts,file.path(out,'paired_landmark_contrasts.csv'))
# Sensitivity point estimates with identical adult/pathology criteria; explicit scenarios.
sens<-list()
for(sc in c('All-primary-sequences','Diagnosis-2004-2015','Non-distant','Unknown-as-cancer','Unknown-as-other','Midpoint-all-months','Myxoid-without-round-cell','Known-stage')){
 z<-if(sc=='All-primary-sequences')copy(r)else copy(x)
 if(sc=='Diagnosis-2004-2015')z<-z[year<=2015]
 if(sc=='Non-distant')z<-z[stage!='Distant']
 if(sc=='Unknown-as-cancer')z[event==3,event:=1L]
 if(sc=='Unknown-as-other')z[event==3,event:=2L]
 if(sc=='Midpoint-all-months')z[,time:=months+.5]
 if(sc=='Myxoid-without-round-cell')z<-z[v7!='8853']
 if(sc=='Known-stage')z<-z[stage!='Unknown']
 z[,id:=.I]
 for(L in c(0,12,36,60)){cc<-prep(z,L,'Overall','Overall',sc);v<-calc(cc,rep(1,nrow(z)));sens[[length(sens)+1]]<-data.table(scenario=sc,landmark=L/12,n=cc$n,cancer=v[1],other=v[2],difference=v[4])}
}
fwrite(rbindlist(sens),file.path(out,'sensitivity_point_estimates.csv'))
writeLines(c('1000 bootstrap replicates; all strata and landmarks share each patient resample.',paste('Date:',Sys.Date()),capture.output(sessionInfo())),file.path(out,'session_info.txt'))
print(rbindlist(flow));print(result[kind=='Overall' & horizon==3 & outcome=='Difference']);cat('COMPLETE\n')

