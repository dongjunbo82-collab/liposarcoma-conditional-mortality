"""Reconstruct landmark follow-up support, using local data only."""
from pathlib import Path
import os
import pandas as pd
P=Path(os.environ.get('LPS_OUTPUT_DIR','outputs/population_study'))
d=pd.read_csv(P/'primary_analysis_cohort.csv',dtype={'patient_id':str})
support=pd.read_csv(P/'landmark_support.csv')
def subset(kind,g):
 if kind=='Overall':return d
 if kind=='Subtype':return d[d.subtype==g]
 if kind=='Site':return d[d.site==g]
 a,b=g.split(' | ');return d[(d.subtype==a)&(d.site==b)]
def counts(z,L,H):
 z=z[z.time>L*12];t=z.time-L*12;h=H*12
 ans=dict(n=len(z),cancer_events=int(((t<=h)&(z.event==1)).sum()),
 other_events=int(((t<=h)&(z.event==2)).sum()),unknown_events=int(((t<=h)&(z.event==3)).sum()),
 censored_before_end=int(((t<h)&(z.event==0)).sum()),
 censored_at_end=int(((t==h)&(z.event==0)).sum()),
 deaths_at_end=int(((t==h)&(z.event!=0)).sum()),
 observed_beyond_horizon=int((t>h).sum()))
 ans['observed_to_end']=int((t>=h).sum())
 ans['alive_observed_through_end']=ans['censored_at_end']+ans['observed_beyond_horizon']
 assert ans['n']==sum(ans[k] for k in ['cancer_events','other_events','unknown_events','censored_before_end','censored_at_end','observed_beyond_horizon'])
 assert ans['observed_to_end']==ans['alive_observed_through_end']+ans['deaths_at_end']
 return ans
rows=[]
for a in support.to_dict('records'):
 v=counts(subset(a['kind'],a['group']),a['landmark'],a['horizon'])
 for k in ['n','cancer_events','other_events','unknown_events','observed_beyond_horizon']:assert v[k]==a[k],(a,k,v[k])
 rows.append(dict(kind=a['kind'],group=a['group'],landmark=a['landmark'],horizon=a['horizon'],**v))
for age in ['<65','65-74','75+']:
 m=(d.age+5<65) if age=='<65' else ((d.age+5>=65)&(d.age+5<75)) if age=='65-74' else (d.age+5>=75)
 rows.append(dict(kind='Attained age',group=age,landmark=5,horizon=3,**counts(d[m],5,3)))

Q=P/'revision/focus';Q.mkdir(parents=True,exist_ok=True)
pd.DataFrame(rows).to_csv(Q/'followup_accounting.csv',index=False)
