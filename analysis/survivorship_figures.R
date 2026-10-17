# R quantitative grid: reveal age-related differences in future mortality, with uncertainty.
# Existing workflow; 183 mm, editable SVG and PNG; no PDF; no inferred crossover lines.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
p<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study");d<-fread(file.path(p,'survivorship/conditional_cause_age_risks.csv'))
z<-d[scheme=='Detailed' & landmark==5 & kind %chin% c('Overall','Attained age') & outcome %chin% c('SEER-attributed cancer','All other attributed','Cardiovascular','Malignant-coded other death')]
z[,group:=factor(group,levels=c('Overall','<65','65-74','75+'))]
z[,panel:=factor(outcome,levels=c('SEER-attributed cancer','All other attributed','Cardiovascular','Malignant-coded other death'),labels=c('SEER-attributed cancer death','All other attributed deaths','Cardiovascular death','Malignant-coded other death'))]
f<-ggplot(z,aes(x=group,y=100*estimate))+
 geom_errorbar(aes(ymin=100*lower,ymax=100*upper),width=.12,linewidth=.45,colour='#276789')+
 geom_point(aes(shape=sparse),size=2.3,colour='#276789')+
 facet_wrap(~panel,ncol=2)+scale_y_continuous(limits=c(0,24),breaks=seq(0,24,5))+
 scale_shape_manual(values=c('FALSE'=16,'TRUE'=1),labels=c('Adequate event support','Fewer than 20 events'))+
 labs(x='Estimated age at the five-year landmark (years)',y='Risk of death within the next 3 years (%)',shape=NULL,
 subtitle='Overall n = 5,980; <65: n = 3,368; 65-74: n = 1,473; 75+: n = 1,139',
 caption='Pointwise 95% intervals from 1,000 patient bootstrap samples. Hollow point: 9 cardiovascular deaths.\nBottom panels partition other-cause mortality; malignant-coded deaths do not establish second primaries.')+
 theme_classic(base_size=9,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='bottom',plot.subtitle=element_text(size=8),plot.caption=element_text(hjust=0,size=7))
ggsave(file.path(p,'figure3_survivorship_age.png'),f,width=183,height=145,units='mm',dpi=400)
ggsave(file.path(p,'figure3_survivorship_age.svg'),f,width=183,height=145,units='mm',device=svglite::svglite)
