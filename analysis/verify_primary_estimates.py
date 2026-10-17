"""Independent cumulative-incidence verification; no documents or network access."""
from pathlib import Path
import os
import numpy as np
import pandas as pd
p=Path(os.environ.get('LPS_OUTPUT_DIR','outputs/population_study'))
d=pd.read_csv(p/'primary_analysis_cohort.csv',dtype={'patient_id':str})
r=pd.read_csv(p/'primary_conditional_risks.csv');error=0.0
for (kind,group,L,H),z in r.groupby(['kind','group','landmark','horizon']):
    x=d
    if kind=='Subtype':x=x[x.subtype==group]
    elif kind=='Site':x=x[x.site==group]
    elif kind=='Subtype x site':
        a,b=group.split(' | ');x=x[(x.subtype==a)&(x.site==b)]
    x=x[x.time>L*12];n=len(x);survival=1.0;incidence=np.zeros(3)
    for t,events in x.groupby('time').event.agg(list).items():
        if t>(L+H)*12:break
        deaths=np.array([events.count(k) for k in [1,2,3]])
        incidence+=survival*deaths/n;survival*=1-deaths.sum()/n;n-=len(events)
    values=dict(zip(['Cancer','Other','Unknown'],incidence));values['Difference']=incidence[1]-incidence[0]
    for a in z.itertuples():error=max(error,abs(a.estimate-values[a.outcome]))
assert error<1e-12,error
print('Independent Aalen-Johansen checks passed.')
