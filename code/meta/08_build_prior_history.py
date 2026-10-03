from pathlib import Path
import os
import pandas as pd
import numpy as np
REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get('BMI_BP_OUTPUT_DIR', REPO_ROOT / 'outputs'))
OUT = ROOT / 'submission_extensions'
OUT.mkdir(exist_ok=True, parents=True)
hrs_root = Path(os.environ['HRS_DATA_DIR'])
elsa_root = Path(os.environ['ELSA_DATA_DIR'])
hrs = hrs_root / 'randhrs1992_2022v1.dta'
cols=['hhidpn']+[f'r{w}hibpe' for w in range(1,13)]+[f'r{w}{s}' for w in range(8,12) for s in ['bpsys','bpdia']]+['r12pmbmi']
h=pd.read_stata(hrs,columns=cols,convert_categoricals=False)
hc=[f'r{w}hibpe' for w in range(1,13)]
h['prior_positive']=h[hc].eq(1).any(axis=1)
h['history_observed']=h[hc].isin([0,1]).any(axis=1)
for w in range(8,12):
    h['prior_positive'] |= h[f'r{w}bpsys'].ge(140)|h[f'r{w}bpdia'].ge(90)
h=h.rename(columns={'hhidpn':'person_id','r12pmbmi':'simple_bmi'})
h[['person_id','prior_positive','history_observed','simple_bmi']].to_csv(OUT/'hrs_history.csv',index=False)
ep = elsa_root
cols=['idauniq']+[f'r{w}{s}' for w in range(1,7) for s in ['hibpe','rxhibp']]
e=pd.read_stata(ep/'gh_elsa_h.dta',columns=cols,convert_categoricals=False)
e['prior_positive']=e[cols[1:]].eq(1).any(axis=1)
e['history_observed']=e[cols[1:]].isin([0,1]).any(axis=1)
for w,name in [(2,'wave_2_nurse_data_v2.dta'),(4,'wave_4_nurse_data.dta')]:
    n=pd.read_stata(ep/name,columns=['idauniq','sysval','diaval'],convert_categoricals=False)
    pos=set(n.loc[n.sysval.ge(140)|n.diaval.ge(90),'idauniq'])
    e['prior_positive'] |= e.idauniq.isin(pos)
e=e.rename(columns={'idauniq':'person_id'})
b=pd.read_csv(ROOT/'elsa_analysis/elsa_baseline_w6.csv')
e=e.merge(b[['person_id','bmi_w6']].rename(columns={'bmi_w6':'simple_bmi'}),on='person_id',how='inner',validate='one_to_one')
e[['person_id','prior_positive','history_observed','simple_bmi']].to_csv(OUT/'elsa_history.csv',index=False)
c=pd.read_csv(ROOT/'chns_analysis/chns_baseline_2006.csv')
c['prior_positive']=False;c['history_observed']=False
for year in [2000,2004]:
    c['prior_positive'] |= c[f'state5_{year}'].isin([3,4,5])|c[f'sbp_{year}'].ge(140)|c[f'dbp_{year}'].ge(90)|c[f'med_current_{year}'].eq(1)
    c['history_observed'] |= c[f'state5_{year}'].notna()
c=c.rename(columns={'bmi_2006':'simple_bmi'})
c[['person_id','prior_positive','history_observed','simple_bmi']].to_csv(OUT/'chns_history.csv',index=False)
print('History tables generated locally; no participant-level data exported externally.')


