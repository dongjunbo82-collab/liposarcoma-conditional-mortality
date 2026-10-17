# Figure contract: demonstrate heterogeneity of competing mortality profiles;
# quantitative comparison, R only; no inferred continuous crossover between landmarks.
# PNG + editable SVG; 1000-bootstrap pointwise intervals; no PDF.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
p<-Sys.getenv("LPS_OUTPUT_DIR", unset="outputs/population_study");d<-fread(file.path(p,'primary_conditional_risks.csv'))
group_order<-c('Well-differentiated','Dedifferentiated','Myxoid/round cell','Pleomorphic','Retroperitoneum','Extremity','Trunk/pelvis soft tissue')
profile_order<-c('Overall',group_order[1:4],'Other/NOS',group_order[5:7])
z<-d[horizon==3 & group %chin% profile_order & kind %chin% c('Overall','Subtype','Site') & outcome %chin% c('Cancer','Other')]
z[,group:=factor(group,levels=profile_order)];z[,outcome:=factor(outcome,levels=c('Cancer','Other'))]
f<-ggplot(z,aes(x=landmark,y=100*estimate,colour=outcome,shape=outcome))+
 geom_errorbar(aes(ymin=100*lower,ymax=100*upper),position=position_dodge(width=.26),width=.15,linewidth=.35)+
 geom_point(position=position_dodge(width=.26),size=1.8)+
 facet_wrap(~group,ncol=3)+scale_x_continuous(breaks=c(0,1,3,5))+
 scale_colour_manual(values=c(Cancer='#A33A35',Other='#276789'),labels=c('SEER-attributed cancer death','Other-cause death'))+
 scale_shape_manual(values=c(Cancer=16,Other=17),labels=c('SEER-attributed cancer death','Other-cause death'))+
 labs(x='Years already survived',y='Conditional risk of death within the next 3 years (%)',colour=NULL,shape=NULL,
 caption='Adults with histologically confirmed first/only primary liposarcoma, 2004-2023.\nAalen-Johansen estimates; pointwise 95% intervals from 1000 patient-level bootstrap samples.')+
 theme_classic(base_size=9,base_family='Arial')+theme(legend.position='bottom',strip.background=element_blank(),strip.text=element_text(face='bold',size=8),plot.caption=element_text(hjust=0,size=7))
ggsave(file.path(p,'figure1_conditional_profiles.png'),f,width=183,height=160,units='mm',dpi=400)
ggsave(file.path(p,'figure1_conditional_profiles.svg'),f,width=183,height=160,units='mm',device=svglite::svglite)
z<-d[horizon==3 & landmark==5 & outcome=='Difference' & ((kind=='Overall') | (group %chin% group_order & kind %chin% c('Subtype','Site')))]
z[,label:=paste0(group,' (n=',n,')')];z[,label:=factor(label,levels=rev(label[order(match(group,c('Overall',group_order)))]))]
f2<-ggplot(z,aes(y=label,x=100*estimate))+
 geom_vline(xintercept=0,linetype='dashed',colour='grey50',linewidth=.4)+
 geom_errorbar(aes(xmin=100*lower,xmax=100*upper),orientation='y',width=.18,colour='#276789',linewidth=.6)+
 geom_point(colour='#276789',size=2.4)+
 labs(y=NULL,x='Other-cause minus index-cancer death risk (percentage points)',title='Mortality balance among five-year survivors',subtitle='Risks over the next 3 years',caption='Aalen-Johansen estimates; pointwise 95% percentile intervals from 1000 patient bootstraps.\nUnadjusted subgroup estimates; histology and site groups overlap. Positive values favour higher other-cause risk.')+
 theme_classic(base_size=9,base_family='Arial')+theme(plot.title=element_text(face='bold',size=11),plot.caption=element_text(hjust=0,size=7))
ggsave(file.path(p,'figure2_five_year_risk_difference.png'),f2,width=183,height=125,units='mm',dpi=400)
ggsave(file.path(p,'figure2_five_year_risk_difference.svg'),f2,width=183,height=125,units='mm',device=svglite::svglite)
# Formal marginal contrasts at year 5, all obtained from shared bootstrap draws.
b<-readRDS(file.path(p,'primary_bootstrap_array.rds'));m<-fread(file.path(p,'landmark_support.csv'))
ct<-list()
for(pair in list(c('Subtype','Dedifferentiated','Well-differentiated'),c('Site','Retroperitoneum','Extremity'))){
 ia<-m[kind==pair[1]&group==pair[2]&landmark==5&horizon==3,index];ib<-m[kind==pair[1]&group==pair[3]&landmark==5&horizon==3,index]
 v<-b[ia,4,]-b[ib,4,];v0<-d[index==ia&horizon==3&outcome=='Difference',estimate]-d[index==ib&horizon==3&outcome=='Difference',estimate]
 ct[[length(ct)+1]]<-data.table(contrast=paste(pair[2],'-',pair[3]),estimand='Difference of other-minus-index-cancer mortality differences',estimate=v0,lower=quantile(v,.025),upper=quantile(v,.975))
}
fwrite(rbindlist(ct),file.path(p,'five_year_subgroup_contrasts.csv'))
