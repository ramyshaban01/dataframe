import pandas as pd, numpy as np, pathlib, subprocess, pycountry, gzip, json, os
ROOT=pathlib.Path(".")
inp=ROOT/"inputs"; work=ROOT/"analysis_work"; work.mkdir(exist_ok=True)

bis_zip=next((inp/"bis_artifact").glob("WS_LBS_D_PUB_csv_flat.zip"))
cepii_zip=next((inp/"cepii_artifact").glob("Gravity_csv_V202211.zip"))
ideal_dir=inp/"ideal_artifact"

subprocess.run(["unzip","-o",str(bis_zip),"-d",str(work/"bis")],check=True)
subprocess.run(["unzip","-o",str(cepii_zip),"-d",str(work/"cepii")],check=True)
bis_csv=next((work/"bis").glob("*.csv"))
cepii_csv=next((work/"cepii").glob("Gravity_V202211.csv"))
ideal=next(ideal_dir.glob("IdealPointDyads1946-2025.tab"))

# Extract only Q4 stocks needed from the 6GB BIS flat CSV. Bank-sector positions
# are reconstructed as all-sector minus non-bank positions (Table A6.2 logic).
header=open(bis_csv,encoding="utf-8").readline().rstrip("\n")
target=work/"bis_target.csv"
with open(target,"w",encoding="utf-8") as out:
    out.write(header+"\n")
prefix_template='Q: Quarterly,S: Amounts outstanding / Stocks,{pos},{instr},TO1: All currencies,A: All currencies (=D+F+U),5J: All countries,"A: All reporting banks/institutions (domestic, foreign, consortium and unclassified)",'
for pos in ["C: Total claims","L: Total liabilities"]:
    for instr in ["A: All instruments","G: Loans and deposits"]:
        p=prefix_template.format(pos=pos,instr=instr)
        cmd=f'''LC_ALL=C rg -F {p!r} {str(bis_csv)!r} | LC_ALL=C rg ',N: Cross-border,(19|20)[0-9]{{2}}-Q4,' >> {str(target)!r} || true'''
        subprocess.run(["bash","-lc",cmd],check=True)

use=['L_POSITION:Balance sheet position','L_INSTR:Type of instruments','L_REP_CTY:Reporting country',
     'L_CP_SECTOR:Counterparty sector','L_CP_COUNTRY:Counterparty country',
     'TIME_PERIOD:Time period or range','OBS_VALUE:Observation Value']
d=pd.read_csv(target,usecols=use,low_memory=False)
code=lambda s:s.astype(str).str.split(':').str[0].str.strip()
d['pos']=code(d[use[0]]); d['instr']=code(d[use[1]]); d['rep2']=code(d[use[2]])
d['sec']=code(d[use[3]]); d['cp2']=code(d[use[4]])
d['year']=pd.to_numeric(d[use[5]].astype(str).str[:4],errors='coerce')
d['value']=pd.to_numeric(d[use[6]],errors='coerce')
d=d[d.rep2.str.fullmatch('[A-Z]{2}',na=False)&d.cp2.str.fullmatch('[A-Z]{2}',na=False)&
    d.sec.isin(['A','N'])&d.pos.isin(['C','L'])&d.instr.isin(['A','G'])&d.year.between(1977,2025)].copy()
manual={'TW':'TWN','XK':'XKX','CS':'CSK','YU':'YUG','SU':'SUN','AN':'ANT','DD':'DDR','ZR':'ZAR',
        'BU':'MMR','TP':'TLS','UK':'GBR','PU':None}
def iso3(c):
    if c in manual: return manual[c]
    x=pycountry.countries.get(alpha_2=c)
    return x.alpha_3 if x else None
d['iso3_o']=d.rep2.map(iso3); d['iso3_d']=d.cp2.map(iso3)
d=d[d.iso3_o.notna()&d.iso3_d.notna()&(d.iso3_o!=d.iso3_d)]

bank_sets=[]
diagnostics={}
for instr,label in [('A','all_instruments'),('G','loans_deposits')]:
    x=d[d.instr==instr]
    piv=x.pivot_table(index=['rep2','cp2','iso3_o','iso3_d','year','pos'],columns='sec',values='value',aggfunc='first').reset_index()
    for s in ['A','N']:
        if s not in piv: piv[s]=np.nan
    piv['bank_component']=piv['A']-piv['N']
    nneg=int((piv.bank_component<0).sum())
    piv.loc[piv.bank_component<0,'bank_component']=np.nan
    C=piv[piv.pos=='C'][['iso3_o','iso3_d','year','bank_component']].rename(columns={'bank_component':'reported_claim_musd'})
    L=piv[piv.pos=='L'][['iso3_o','iso3_d','year','bank_component']].rename(
        columns={'iso3_o':'iso3_d','iso3_d':'iso3_o','bank_component':'reverse_liability_musd'})
    m=C.merge(L,on=['iso3_o','iso3_d','year'],how='outer')
    m['bank_musd']=m.reported_claim_musd.combine_first(m.reverse_liability_musd)
    m['data_source']=np.where(m.reported_claim_musd.notna(),'reported_claim',
                       np.where(m.reverse_liability_musd.notna(),'mirrored_liability','missing'))
    m=m[m.bank_musd.notna()].copy(); m['measure']=label
    both=m.dropna(subset=['reported_claim_musd','reverse_liability_musd'])
    diagnostics[label]={'observations':int(len(m)),'directed_pairs':int(m[['iso3_o','iso3_d']].drop_duplicates().shape[0]),
                        'negative_differences_dropped':nneg,'direct_mirror_overlap':int(len(both)),
                        'direct_mirror_corr':float(both[['reported_claim_musd','reverse_liability_musd']].corr().iloc[0,1]) if len(both)>2 else None}
    bank_sets.append(m)
bank=pd.concat(bank_sets,ignore_index=True)
bank.to_csv(work/"bank_interbank_q4.csv",index=False)
json.dump(diagnostics,open(work/"bank_diagnostics.json","w"),indent=2)

# CEPII Gravity / DOTS trade. Origin-reported DOTS is preferred, destination report fills missing.
parts=[]
for ch in pd.read_csv(cepii_csv,usecols=['year','iso3_o','iso3_d','tradeflow_imf_o','tradeflow_imf_d'],chunksize=250000,low_memory=False):
    ch=ch[ch.year.between(1978,2020)&ch.iso3_o.notna()&ch.iso3_d.notna()&(ch.iso3_o!=ch.iso3_d)].copy()
    ch['trade_o']=pd.to_numeric(ch.tradeflow_imf_o,errors='coerce')
    ch['trade_d']=pd.to_numeric(ch.tradeflow_imf_d,errors='coerce')
    ch['trade_musd']=ch.trade_o.combine_first(ch.trade_d)/1000.0
    parts.append(ch[['year','iso3_o','iso3_d','trade_o','trade_d','trade_musd']])
trade=pd.concat(parts,ignore_index=True).drop_duplicates(['year','iso3_o','iso3_d'])

# BSV/Voeten dyadic ideal-point distance; one-year lag relative to economic outcome year.
ip=pd.read_csv(ideal,sep='\t',usecols=['iso3c1','iso3c2','year','AbsIdealDiff'])
ip['year']=pd.to_numeric(ip.year,errors='coerce')
ip=ip.dropna()
ip['a']=ip[['iso3c1','iso3c2']].min(axis=1); ip['b']=ip[['iso3c1','iso3c2']].max(axis=1)
ip['econ_year']=ip.year.astype(int)+1
ip=ip[['a','b','econ_year','AbsIdealDiff']].drop_duplicates(['a','b','econ_year'])

b=bank[(bank.measure=='all_instruments')&bank.year.between(1978,2020)].copy()
common=b.merge(trade,on=['year','iso3_o','iso3_d'],how='left',validate='one_to_one')
common['a']=common[['iso3_o','iso3_d']].min(axis=1);common['b']=common[['iso3_o','iso3_d']].max(axis=1)
common=common.merge(ip,left_on=['a','b','year'],right_on=['a','b','econ_year'],how='left').drop(columns='econ_year')
common=common[common.bank_musd.notna()&common.trade_musd.notna()&common.AbsIdealDiff.notna()&
              (common.bank_musd>=0)&(common.trade_musd>=0)].copy()
common['pair_id']=common.iso3_o+'_'+common.iso3_d
common['op_year']=common.iso3_o+'_'+common.year.astype(str)
common['dp_year']=common.iso3_d+'_'+common.year.astype(str)
common['undirected_pair']=common[['iso3_o','iso3_d']].min(axis=1)+'_'+common[['iso3_o','iso3_d']].max(axis=1)

# Exact matched estimation support: remove fixed-effect groups with no positive banking outcome,
# then use precisely those observations for both bank and trade PPML.
while True:
    bad=np.zeros(len(common),dtype=bool)
    for c in ['pair_id','op_year','dp_year']:
        bad |= common.groupby(c)['bank_musd'].transform('sum').le(0).to_numpy()
    if not bad.any(): break
    common=common.loc[~bad].copy()
mean=common.AbsIdealDiff.mean(); sd=common.AbsIdealDiff.std(ddof=0)
common['ipd_z']=(common.AbsIdealDiff-mean)/sd
common=common.sort_values(['year','iso3_o','iso3_d'])
common.to_csv(work/"common_matched.csv",index=False)
summary={'nobs':int(len(common)),'directed_pairs':int(common.pair_id.nunique()),
         'undirected_pairs':int(common.undirected_pair.nunique()),
         'countries':int(len(set(common.iso3_o)|set(common.iso3_d))),
         'year_min':int(common.year.min()),'year_max':int(common.year.max()),
         'ipd_mean':float(mean),'ipd_sd':float(sd),
         'bank_zero_share':float((common.bank_musd==0).mean()),
         'coverage_by_year':{str(int(y)):int(n) for y,n in common.groupby('year').size().items()}}
json.dump(summary,open(work/"common_summary.json","w"),indent=2)
print(json.dumps(summary,indent=2))
