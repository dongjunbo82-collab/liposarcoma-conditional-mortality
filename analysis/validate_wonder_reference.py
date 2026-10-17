import os
from pathlib import Path
from lxml import etree as E
import pandas as pd
import json
p=(Path(os.environ.get("LPS_OUTPUT_DIR", "outputs/population_study"))/"revision/wonder_20261016")
rows=[]
for db,span in [("D76","2009_2017"),("D158","2018_2023")]:
 for suffix in ["","_ten_year"]:
  for row in E.parse(str(p/(db+"_"+span+suffix+"_response.xml"))).findall(".//data-table/r"):
   c=row.findall("c"); assert len(c)==6
   vals=[x.get("l",x.get("v")).strip() for x in c]
   age=vals[2]
   if suffix and age!="85+ years":continue
   if not suffix and age not in [str(a)+"-"+str(a+4)+" years" for a in range(20,85,5)]:continue
   rows.append([db]+vals)
d=pd.DataFrame(rows,columns=["database","year","sex","age_group","deaths","population","rate_per_100000"])
for col in ["year","deaths","population","rate_per_100000"]:
 d[col]=pd.to_numeric(d[col].astype(str).str.replace(",","",regex=False),errors="raise")
assert len(d)==15*2*14
assert not d.duplicated(["year","sex","age_group"]).any()
assert (d.groupby(["year","sex"]).size()==14).all()
assert (d.population>0).all()
assert ((d.deaths/d.population*100000-d.rate_per_100000).abs()<0.000006).all()
d["rate_per_person_year"]=d.deaths/d.population
d.to_csv(p/"cdc_cardiovascular_reference_2009_2023.csv",index=False)
report={"retrieved":"2026-10-16","status":"reference_rates_retrieved_and_validated_SMR_not_yet_calculated","cells":len(d),"years":[2009,2023],"sexes":["Female","Male"],"ages":"20-24 through 80-84 in five-year bands, plus 85+","icd10_definition":"I00-I13, I20-I51, I60-I78; selectable underlying-cause ICD categories","validation":"All 420 cells numeric, unique, complete; counts/population reproduce published crude rates within rounding","limitations":["National all-race reference; no geographic matching","85+ open age band; residual age confounding possible","SEER diagnosis month unavailable in current export; calendar-time assignment needs prespecified approximation and sensitivity analysis","Observed cardiovascular deaths must match death-certificate definition rather than only SEER other-cause subset"],"error_correction":"Prior HTTP 500 responses included invalid request codes/options, not only server outage. D176 was not the correct expanded underlying-cause endpoint; D158 used. HTTP 429 explicitly requires >=15 seconds spacing; script now uses 16 seconds."}
(p/"retrieval_validation.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print(json.dumps(report,indent=2))
