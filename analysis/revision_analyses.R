suppressPackageStartupMessages({library(data.table);library(survival)})
p<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study");out<-file.path(p,'revision');dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Reuse only the tested estimators, without rerunning the primary analysis.
for(e in parse('analysis/population_primary.R'))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])[1]%in%c('prep','calc'))eval(e)
r<-fread(Sys.getenv("SEER_EXPORT", unset="inputs/seer_export.tsv"),colClasses='character');setnames(r,names(r),paste0('v',1:23))
r[,`:=`(age=as.integer(sub('[^0-9].*','',v3)),year=as.integer(v2),months=suppressWarnings(as.integer(v19)),rank=fifelse(v10=='One primary only',1L,suppressWarnings(as.integer(sub('^(\\d+).*','\\1',v10)))),ord=.I)]
setorder(r,v1,rank,year,ord);r<-r[!duplicated(v1)]
r<-r[year>=2004 & age>=18 & v9%chin%c('Positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies') & !is.na(months)]
r[,`:=`(first=v10%chin%c('One primary only','1st of 2 or more primaries'),event=fcase(grepl('attributable to this cancer',v21,fixed=TRUE),1L,grepl('causes other than this cancer',v22,fixed=TRUE),2L,grepl('missing/unknown COD',v21,fixed=TRUE),3L,default=0L),time=pmax(months,.5))]
r[,subtype:=fcase(v7=='8851','Well-differentiated',v7=='8858','Dedifferentiated',v7%chin%c('8852','8853'),'Myxoid/round cell',v7=='8854','Pleomorphic',default='Other/NOS')]
r[,sitecode:=sub('-.*','',v6)];r[,site:=fcase(sitecode=='C48.0','Retroperitoneum',sitecode%chin%c('C49.1','C49.2'),'Extremity',sitecode%chin%c('C49.3','C49.4','C49.5','C49.6'),'Trunk/pelvis soft tissue',default='Other')]
r[,stage:=fcase(v14=='Localized only','Localized',grepl('^Regional',v14),'Regional',grepl('^Distant',v14),'Distant',default='Unknown')]
r[,id:=.I];N<-nrow(r);x<-copy(r[first==TRUE]);stopifnot(nrow(x)==12219)
old<-fread(file.path(p,'primary_analysis_cohort.csv'),colClasses=c(patient_id='character'));stopifnot(identical(x$v1,old$patient_id),identical(x$event,old$event))
x[,surg_raw:=fifelse(year<=2022,v15,v16)]
x[,surg_code:=suppressWarnings(as.integer(fifelse(year<=2022,surg_raw,sub('^[AB]([0-9]{2})0$','\\1',surg_raw))))]
# Receipt of a procedure is broader than a documented resection.
x[,surgery:=fcase(surg_code==0,'No recorded surgery',surg_code>=10 & surg_code<=90,'Recorded surgery/procedure',default='Unknown/not classifiable')]
x[,strict_resection:=(site=='Retroperitoneum' & surg_code%in%c(30L,40L,60L)) | (site=='Extremity' & surg_code%in%c(25L,26L,30L,40L,41L,42L,50L,51L,52L,53L,54L))]
x[,broad_resection:=strict_resection | (site=='Retroperitoneum' & surg_code%in%c(20:27,50))]
x[,age5:=fifelse(age+5<65,'<65',fifelse(age+5<75,'65-74','75+'))]
fwrite(x[,.(n=.N,five_year_survivors=sum(time>60)),by=.(site,year,surg_raw,surg_code,surgery,strict_resection,broad_resection)],file.path(out,'surgery_coding_audit.csv'))
fwrite(x[,.(n=.N,five_year_survivors=sum(time>60)),by=surgery],file.path(out,'surgery_counts.csv'))
cache<-list();add<-function(z,L,g,k,sc='Revision'){cc<-prep(z,L,g,k,sc);if(!is.null(cc))cache[[length(cache)+1]]<<-cc}
scenarios<-c('All-primary-sequences','Diagnosis-2004-2015','Non-distant','Unknown-as-cancer','Unknown-as-other','Midpoint-all-months','Myxoid-without-round-cell','Known-stage')
for(sc in scenarios){
 z<-if(sc=='All-primary-sequences')copy(r)else copy(x)
 if(sc=='Diagnosis-2004-2015')z<-z[year<=2015]
 if(sc=='Non-distant')z<-z[stage!='Distant']
 if(sc=='Unknown-as-cancer')z[event==3,event:=1L]
 if(sc=='Unknown-as-other')z[event==3,event:=2L]
 if(sc=='Midpoint-all-months')z[,time:=months+.5]
 if(sc=='Myxoid-without-round-cell')z<-z[v7!='8853']
 if(sc=='Known-stage')z<-z[stage!='Unknown']
 for(L in c(0,12,36,60))add(z,L,sc,'Sensitivity',sc)
}
for(g in unique(x$surgery))add(x[surgery==g],60,g,'Surgery')
for(site_i in c('Retroperitoneum','Extremity')){
 for(L in c(0,12,36,60))add(x[site==site_i & stage%chin%c('Localized','Regional') & strict_resection],L,site_i,'Resected nonmetastatic')
 add(x[site==site_i & stage%chin%c('Localized','Regional') & broad_resection],60,site_i,'Broad resection nonmetastatic')
}
for(s in c('Well-differentiated','Dedifferentiated','Myxoid/round cell','Pleomorphic'))for(t in c('Retroperitoneum','Extremity'))for(a in c('<65','65-74','75+'))add(x[subtype==s & site==t & age5==a],60,paste(s,t,a,sep=' | '),'Joint age subtype site')
point<-t(vapply(cache,calc,numeric(8),w=rep(1,N)))
maxerr<-0
for(i in seq_along(cache)){
 cc<-cache[[i]];fit<-survfit(Surv(cc$relative,factor(cc$event,levels=0:3))~1)
 v<-summary(fit,times=36,extend=TRUE)$pstate[1,-1];maxerr<-max(maxerr,max(abs(v-point[i,1:3])))
}
stopifnot(maxerr<1e-10)
oldsen<-fread(file.path(p,'sensitivity_point_estimates.csv'))
for(i in which(vapply(cache,function(z)z$kind=='Sensitivity',logical(1)))){cc<-cache[[i]];v<-oldsen[scenario==cc$group & landmark==cc$L/12];stopifnot(max(abs(point[i,c(1,2,4)]-unlist(v[,.(cancer,other,difference)])))<1e-10)}
B<-1000L;set.seed(20261009);boot<-array(NA_real_,c(length(cache),4,B))
for(b in seq_len(B)){
 w<-tabulate(sample.int(N,N,replace=TRUE),nbins=N)
 boot[,,b]<-t(vapply(cache,function(cc)if(sum(w[cc$id])==0)rep(NA_real_,4)else calc(cc,w)[1:4],numeric(4)))
 if(b%%100==0)cat('Revision bootstrap',b,'/',B,'\n')
}
res<-rbindlist(lapply(seq_along(cache),function(i){cc<-cache[[i]];events<-tabulate(cc$event[cc$relative<=36],nbins=3);beyond<-sum(cc$relative>36);rbindlist(lapply(1:4,function(j){data.table(kind=cc$kind,group=cc$group,landmark=cc$L/12,n=cc$n,outcome=c('Cancer','Other','Unknown','Difference')[j],events=if(j<4)events[j]else NA_integer_,cancer_events=events[1],other_events=events[2],beyond=beyond,sparse20=cc$n<100 || min(events[1:2])<20 || beyond<20,sparse10=cc$n<100 || min(events[1:2])<10 || beyond<20,estimate=point[i,j],lower=quantile(boot[i,j,],.025,na.rm=TRUE),upper=quantile(boot[i,j,],.975,na.rm=TRUE),valid_bootstrap=sum(is.finite(boot[i,j,])))}))}))
fwrite(res,file.path(out,'revision_conditional_risks.csv'));saveRDS(boot,file.path(out,'bootstrap.rds'))
writeLines(c(paste('Independent survfit validation maximum absolute error',maxerr),'Existing sensitivity point estimates reproduced within 1e-10.','1000 shared patient bootstrap resamples of eligible all-sequence universe; seed 20261009.','Exploratory descriptive extension. Pointwise intervals, no formal multiple pairwise testing.','Surgery codes 30/40/60 for strict retroperitoneal resection; specified excision/resection/amputation codes for extremity.','Broad retroperitoneal sensitivity adds 20-27 and 50; no inference about margins or curative intent.',capture.output(sessionInfo())),file.path(out,'verification.txt'))
print(res[landmark==5 & kind%chin%c('Surgery','Resected nonmetastatic') & outcome!='Unknown']);cat('COMPLETE\n')
