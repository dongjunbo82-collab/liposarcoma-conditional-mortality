import os
from pathlib import Path
import csv,re,json
P=(Path(os.environ.get('LPS_OUTPUT_DIR', 'outputs/population_study'))/'survivorship')
P.mkdir(exist_ok=True)
def norm(s):return re.sub(r'[^a-z0-9]','',s.lower())
mapping={}
for line in (P/'seer_cod_recode.txt').read_text().splitlines():
 parts=line.split(';')
 if len(parts)>=5:
  code=re.sub(r'[^0-9]','',parts[4])
  if len(code)==5:mapping[norm(parts[0])]=int(code)
mapping.update({norm('Acute Myeloid Leukemia'):35021,norm('Other Infectious and Parasitic Diseases including HIV'):50040,norm('Alive'):0,norm('State DC not available or state DC available but no COD'):41000,norm('Unknown/missing/invalid COD'):99999})
with open(os.environ.get('SEER_EXPORT', 'inputs/seer_export.tsv'),encoding='utf-8-sig') as f:
 rows=list(csv.reader(f,delimiter='\t'))
labels=sorted({r[22] for r in rows[1:]})
missing=[x for x in labels if norm(x) not in mapping]
assert not missing,missing
with (P/'cod_label_mapping.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['cod_label','cod_code','certificate_group'])
 for label in labels:
  c=mapping[norm(label)]
  group=('Alive' if c==0 else 'Malignant cancer' if 20000<=c<=37000 else 'Cardiovascular' if 50060<=c<=50110 else 'Unknown COD' if c in (41000,99999) else 'Residual or ill-defined' if c in (50200,50300) else 'Other specified nonmalignant')
  w.writerow([label,c,group])
(P/'analysis_addendum.json').write_text(json.dumps({'date':'2026-10-07','status':'Secondary exploratory extension after primary results and initial event-count inspection; not preregistered','cohort':'Unchanged main cohort, N=12219','landmarks_years':[0,1,3,5],'future_window_years':3,'age_at_landmark':'diagnosis age + landmark years; grouped <65,65-74,75+; diagnosis ages 90+ topcoded, all remain75+','primary_decomposition':'Retain SEER-attributed cancer death; partition other deaths by COD recode into malignant-coded, cardiovascular, other specified nonmalignant, residual/ill-defined; unknown retained','sensitivity':'Independent death-certificate grouping: all malignant cancer, cardiovascular, other specified nonmalignant, residual/ill-defined, unknown; does not estimate confirmed second primaries','bootstrap':1000,'seed':20261009,'intervals':'pointwise percentile, shared patient resamples across groups; no multiplicity correction','support':'flag fewer than20 events or fewer than20 observed beyond horizon or N<100','limitations':'No raw ICD death codes; cannot recreate individual SEER attribution or adjudicate true cause; no recurrence or completion-of-treatment data'},indent=2),encoding='utf-8')
ap=P/'analysis_addendum.json'
plan=json.loads(ap.read_text(encoding='utf-8'));plan['seed']=20261008;plan['resampling']='Same patient ordering and random draws as original main analysis to preserve shared estimates and intervals'
ap.write_text(json.dumps(plan,indent=2),encoding='utf-8')
print('Mapped',len(labels),'labels; wrote secondary analysis addendum')
