import os
"""Exploratory external CVD mortality comparison, years 5-8 after diagnosis.

Only public aggregate reference data are external. No individual records uploaded.
Run from project root. Expected counts integrate piecewise-constant population
rates exactly across calendar-year and five-year attained-age boundaries.
"""
from pathlib import Path
from functools import lru_cache
import json, hashlib, subprocess, io
import numpy as np
import pandas as pd
ci_csv=subprocess.check_output([os.environ.get('RSCRIPT', 'Rscript'),'-e',
 'o<-0:1000;write.csv(data.frame(observed=o,lower=ifelse(o==0,0,qchisq(.025,2*o)/2),upper=qchisq(.975,2*(o+1))/2),stdout(),row.names=FALSE)'],text=True)
count_ci=pd.read_csv(io.StringIO(ci_csv)).set_index('observed')

P=Path(os.environ.get('LPS_OUTPUT_DIR', 'outputs/population_study')); O=P/'revision/smr';O.mkdir(exist_ok=True)
rf=P/'revision/wonder_20261016/cdc_cardiovascular_reference_2009_2023.csv'
rates=pd.read_csv(rf)
lookup={(int(r.year),r.sex,85 if r.age_group=='85+ years' else int(r.age_group[:2])):r.rate_per_person_year for r in rates.itertuples()}
x=pd.read_csv(P/'primary_analysis_cohort.csv',dtype={'patient_id':str}).merge(
 pd.read_csv(P/'survivorship/local_patient_event_mapping.csv',dtype={'patient_id':str}),on='patient_id',validate='one_to_one')
assert len(x)==12219 and (x.time>60).sum()==5980
assert x.sex.isin(['Male','Female']).all()
x['cvd']=x.cod_code.isin([50060,50070,50080,50090,50100,50110])
x['unknown']=x.cod_code.isin([41000,99999])
assert ((x.certificate==2)==x.cvd).all()

plan={
 'date':'2026-10-16','status':'Exploratory extension specified after original conditional-risk results; not preregistered',
 'target':'Observed survivors beyond 60 months; follow months >60 to <=96, stop at any death or censoring',
 'endpoint':'Death-certificate cardiovascular recodes 50060,50070,50080,50090,50100,50110, irrespective of SEER attribution',
 'reference':'US national, all races, sex/year/five-year ages 20-84 and 85+, 2009-2023; D76 to 2017, D158 from 2018',
 'main_time':'Original integer survival months; main age at diagnosis = recorded completed age + 0.5; attained age updates continuously',
 'main_calendar':'Average expected counts across equally weighted diagnosis-month midpoints compatible with the full recorded follow-up ending before 2024; calendar time is never silently clamped to 2023',
 'intervals':'Exact two-sided 95% Poisson intervals for observed deaths, with expected counts and person-years treated as fixed',
 'aer':'(observed-expected)/person-years * 10000; interval obtained by replacing observed with its exact Poisson count limits',
 'age_strata':'Age at five-year landmark using recorded integer age + 5: <65,65-74,75+ (fixed across sensitivity analyses)',
 'sensitivities':['Earliest/latest compatible diagnosis month','Age fraction 0 and 11/12 instead of 0.5','Joint calendar/age corners','Diagnosis <=2015, with potential full eight-year follow-up by 2023','Administrative censoring at attained age 85','Exclude diagnosis ages topcoded 90+','Add 0.5 to every survival month and reconstruct landmark eligibility','Assign all unknown-cause deaths in the window to cardiovascular death (extreme bound)'],
 'interpretation':'Comparison with general population, not cancer-free controls or treatment causal effects; no individual linkage, no geographic or race matching; no multiplicity-adjusted subgroup tests',
 'post_result_extension':'After observing the exploratory pleomorphic subgroup SMR, the same sensitivity scenarios were additionally tabulated for that subgroup; this is a post hoc robustness check, not confirmatory testing.',
 'reference_sha256':hashlib.sha256(rf.read_bytes()).hexdigest()}
(O/'analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')

@lru_cache(None)
def integrate(year,age,sex,end,dx,phase):
 """Exact integral on diagnosis-relative years [5,end]."""
 t=5.; ans=0.
 while t < end-1e-10:
  calendar=year+dx+t; attained=age+phase+t
  yr=int(np.floor(calendar+1e-9));ag=min(85,int(np.floor((attained+1e-9)/5))*5)
  assert (yr,sex,ag) in lookup,(yr,sex,ag,calendar,end)
  cy=yr+1-year-dx
  ca=(ag+5-age-phase) if ag<85 else float('inf')
  stop=min(end,cy,ca)
  assert stop>t+1e-10,(t,stop)
  ans+=(stop-t)*lookup[(yr,sex,ag)];t=stop
 return ans

def person_table(calendar='average',phase=.5,month_shift=0.,older=False,censor85=False,exclude90=False):
 z=x.copy();z['follow'] = z.months+month_shift
 z=z[z.follow>60].copy()
 if older:z=z[z.year<=2015]
 if exclude90:z=z[z.age<90]
 out=[]
 for r in z.itertuples():
  end=min(r.follow/12,8.)
  if censor85:end=min(end,85-r.age-phase)
  if end<=5:continue
  candidates=[(m+.5)/12 for m in range(12) if r.year+(m+.5)/12+r.follow/12<=2024+1e-9]
  assert candidates,(r.patient_id,r.year,r.follow)
  offsets=candidates if calendar=='average' else [min(candidates) if calendar=='earliest' else max(candidates)]
  expected=np.mean([integrate(int(r.year),int(r.age),r.sex,end,d,phase) for d in offsets])
  death_in_window=(r.event>0 and r.follow/12<=end+1e-10)
  out.append((r.patient_id,r.year,r.age,r.sex,r.subtype,r.site,
   '<65' if r.age+5<65 else '65-74' if r.age+5<75 else '75+',
   end-5,expected,int(death_in_window and r.cvd),int(death_in_window and r.unknown),
   len(candidates),min(candidates),max(candidates)))
 return pd.DataFrame(out,columns=['patient_id','diagnosis_year','diagnosis_age','sex','subtype','site','landmark_age_group','person_years','expected','observed','unknown_deaths','feasible_months','earliest_offset','latest_offset'])

def estimate(z,scenario,group='Overall',unknown_as_cvd=False):
 n=len(z);o=int(z.observed.sum()+(z.unknown_deaths.sum() if unknown_as_cvd else 0));e=z.expected.sum();py=z.person_years.sum()
 lo=count_ci.loc[o,'lower'];hi=count_ci.loc[o,'upper']
 return dict(scenario=scenario,group=group,n=n,observed=o,expected=e,person_years=py,unknown_deaths=int(z.unknown_deaths.sum()),
  smr=o/e,smr_lower=lo/e,smr_upper=hi/e,observed_rate_per10000=o/py*10000,expected_rate_per10000=e/py*10000,
  aer_per10000=(o-e)/py*10000,aer_lower=(lo-e)/py*10000,aer_upper=(hi-e)/py*10000,sparse=o<20)

main=person_table();assert len(main)==5980 and main.observed.sum()==91
assert abs(main.person_years.sum()-14938.833333333333)<1e-8
main.to_csv(O/'local_person_time_expected.csv',index=False)
results=[estimate(main,'Primary')]
for key in ['landmark_age_group','sex','subtype','site']:
 for name,z in main.groupby(key):results.append(estimate(z,'Primary',key+': '+name))
pd.DataFrame(results).to_csv(O/'primary_and_subgroups.csv',index=False)
sens=[results[0]]
scenarios=[('Earliest compatible diagnosis month',dict(calendar='earliest')),('Latest compatible diagnosis month',dict(calendar='latest')),
 ('Age fraction zero',dict(phase=0.)),('Age fraction eleven twelfths',dict(phase=11/12)),
 ('Early diagnosis and age fraction zero',dict(calendar='earliest',phase=0.)),('Early diagnosis and age fraction eleven twelfths',dict(calendar='earliest',phase=11/12)),
 ('Late diagnosis and age fraction zero',dict(calendar='latest',phase=0.)),('Late diagnosis and age fraction eleven twelfths',dict(calendar='latest',phase=11/12)),
 ('Diagnosis 2004-2015',dict(older=True)),('Censor at attained age 85',dict(censor85=True)),('Exclude diagnosis age 90+',dict(exclude90=True)),
 ('Survival months plus 0.5',dict(month_shift=.5))]
for name,args in scenarios:
 z=person_table(**args)
 for group in ['Overall','<65','65-74','75+','Pleomorphic']:
  zz=z if group=='Overall' else z[z.subtype==group] if group=='Pleomorphic' else z[z.landmark_age_group==group]
  sens.append(estimate(zz,name,group))
for group in ['Overall','<65','65-74','75+','Pleomorphic']:
 zz=main if group=='Overall' else main[main.subtype==group] if group=='Pleomorphic' else main[main.landmark_age_group==group]
 sens.append(estimate(zz,'All unknown deaths assigned to cardiovascular',group,True))
pd.DataFrame(sens).to_csv(O/'sensitivity_results.csv',index=False)
# Audit certificate endpoint against the original attribution-defined component.
window=x[(x.time>60)&(x.time<=96)&(x.event>0)]
window.groupby(['event','certificate']).size().rename('n').reset_index().to_csv(O/'death_definition_audit.csv',index=False)
audit={'cohort_n':len(x),'five_year_n':len(main),'cvd_deaths':int(main.observed.sum()),'unknown_deaths':int(main.unknown_deaths.sum()),
 'person_years':main.person_years.sum(),'expected':main.expected.sum(),'min_feasible_months':int(main.feasible_months.min()),
 'fewer_than_12_feasible_months':int((main.feasible_months<12).sum()),'topcoded_diagnosis_age_n':int((main.diagnosis_age>=90).sum()),
 'calendar_reference_min':2009,'calendar_reference_max':2023,'intervals':'Poisson, conditional on expected counts',
 'rate_validation_file':'../wonder_20261016/retrieval_validation.json'}
(O/'analysis_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(pd.DataFrame(results).head(6).to_string(index=False))
print(pd.DataFrame(sens).query('group=="Overall"').to_string(index=False))
