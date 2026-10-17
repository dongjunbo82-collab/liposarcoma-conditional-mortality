# Complete histology-site display; no outcome-based selection.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
p <- Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study")
d <- fread(file.path(p,'primary_conditional_risks.csv'))
z <- d[horizon==3 & kind=='Subtype x site' & outcome %chin% c('Cancer','Other')]
# Fixed string splitting avoids regular-expression ambiguity in group labels.
z[,histology:=factor(vapply(strsplit(group,' | ',fixed=TRUE),`[`,character(1),1),levels=c('Well-differentiated','Dedifferentiated','Myxoid/round cell','Pleomorphic'))]
z[,site:=factor(vapply(strsplit(group,' | ',fixed=TRUE),`[`,character(1),2),levels=c('Retroperitoneum','Extremity','Trunk/pelvis soft tissue'),labels=c('Retroperitoneum','Extremity','Trunk/pelvis'))]
z[,outcome:=factor(outcome,levels=c('Cancer','Other'))]
flags<-unique(z[sparse==TRUE,.(histology,site,landmark)])
ytop<-ceiling(max(z$upper)*100/10)*10
f <- ggplot(z,aes(landmark,100*estimate,colour=outcome,shape=outcome))+
 geom_errorbar(aes(ymin=ifelse(lower==0 & upper==0,NA_real_,100*lower),ymax=ifelse(lower==0 & upper==0,NA_real_,100*upper)),position=position_dodge(.25),width=.12,linewidth=.35,na.rm=TRUE)+
 geom_point(position=position_dodge(.25),size=1.5)+facet_grid(histology~site,switch='y')+
 geom_text(data=flags,aes(x=landmark,y=ytop-2,label='*'),inherit.aes=FALSE,size=3)+
 scale_x_continuous(breaks=c(0,1,3,5))+scale_y_continuous(limits=c(0,ytop),breaks=seq(0,ytop,20))+
 scale_colour_manual(values=c(Cancer='#A33A35',Other='#276789'),labels=c('SEER-attributed cancer death','Other-cause death'))+
 scale_shape_manual(values=c(Cancer=16,Other=17),labels=c('SEER-attributed cancer death','Other-cause death'))+
 labs(x='Years already survived since diagnosis',y='Risk of death within the next 3 years (%)',colour=NULL,shape=NULL,
 caption='Pointwise 95% bootstrap intervals. * Limited support at the indicated landmark.\nAll eligible stages and treatments; each landmark uses its observed survivor population.')+
 theme_classic(base_size=8,base_family='Arial')+
 theme(legend.position='bottom',strip.background=element_blank(),strip.text=element_text(face='bold',size=8),
 strip.placement='outside',plot.caption=element_text(hjust=0,size=7),legend.text=element_text(size=7),panel.spacing=grid::unit(.7,'lines'))
ggsave(file.path(p,'figure1_conditional_profiles.png'),f,width=183,height=195,units='mm',dpi=400)
ggsave(file.path(p,'figure1_conditional_profiles.svg'),f,width=183,height=195,units='mm',device=svglite::svglite)
fwrite(z,file.path(p,'revision/focus/figure1_source.csv'))
# Figure 2 retains the age evidence but makes the subset relationship explicit.
d<-fread(file.path(p,'survivorship/conditional_cause_age_risks.csv'))
z<-d[scheme=='Detailed' & landmark==5 & kind %chin% c('Overall','Attained age') & outcome %chin% c('SEER-attributed cancer','All other attributed','Cardiovascular','Malignant-coded other death')]
z[,group:=factor(group,levels=c('Overall','<65','65-74','75+'))]
z[,panel:=factor(outcome,levels=c('SEER-attributed cancer','All other attributed','Cardiovascular','Malignant-coded other death'),labels=c('A  SEER-attributed cancer death','B  Total other-cause death','C  Cardiovascular death\nComponent of B','D  Malignant-coded other death\nComponent of B'))]
f<-ggplot(z,aes(group,100*estimate))+geom_errorbar(aes(ymin=100*lower,ymax=100*upper),width=.12,linewidth=.45,colour='#276789')+
 geom_point(aes(shape=sparse),size=2.3,colour='#276789')+facet_wrap(~panel,ncol=2)+
 scale_y_continuous(limits=c(0,24),breaks=seq(0,24,5))+
 scale_shape_manual(values=c('FALSE'=16,'TRUE'=1),labels=c('At least 20 events','Fewer than 20 events'))+
 labs(x='Estimated age at the five-year landmark (years)',y='Risk of death within the next 3 years (%)',shape=NULL,
 subtitle='Overall n = 5,980; <65: n = 3,368; 65-74: n = 1,473; 75+: n = 1,139',
 caption='C and D are selected components of B: do not add them to B. Other components are not plotted.\nPointwise 95% bootstrap intervals. Five-year survival does not establish disease-free status.')+
 theme_classic(base_size=9,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=8),legend.position='bottom',plot.subtitle=element_text(size=8),plot.caption=element_text(hjust=0,size=7))
ggsave(file.path(p,'figure3_survivorship_age.png'),f,width=183,height=145,units='mm',dpi=400)
ggsave(file.path(p,'figure3_survivorship_age.svg'),f,width=183,height=145,units='mm',device=svglite::svglite)
fwrite(z,file.path(p,'revision/focus/figure2_source.csv'))
