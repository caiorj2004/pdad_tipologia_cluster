# -*- coding: utf-8 -*-
# Variáveis de apoio à produção, fora do modelo: proporções por grupo e subpopulação, Kruskal-Wallis por segmento e agregados por RA.
# Entrada: bases corrigidas (programa 01), segmentos_completo.csv (programa 03)
# Saída: apoio_prop.csv, apoio_kw.csv, apoio_ra.csv
import os, re, numpy as np, pandas as pd
from scipy.stats import kruskal
ID,RA,SUB,EST='domicilio_id','a01ra','a01subpopulacao','a01estratoV2'
BASE=os.environ.get('TCC_BASE','.')   # pasta com as bases corrigidas (saída do programa 01)
dom=pd.read_excel(os.path.join(BASE,'pdadr_domicilios.xlsx')); mor=pd.read_excel(os.path.join(BASE,'pdadr_moradores.xlsx'))
for df in (dom,mor):
    df.columns=df.columns.astype(str).str.strip()
    for c in (RA,SUB,EST): df[c]=df[c].astype('string').str.replace(r'\s+',' ',regex=True).str.strip().str.replace('Assentamento/ Agrovila','Assentamento/Agrovila',regex=False)
N_H={"Agricultura Empresarial - Leste":5087,"Agricultura Empresarial - Norte":2968,"Agricultura Empresarial - Oeste":2613,"Agricultura Empresarial - Sul":4663,"Agricultura Familiar - Central Adjacente I":595,"Agricultura Familiar - Central Adjacente II":688,"Agricultura Familiar - Leste":1483,"Agricultura Familiar - Norte":6097,"Agricultura Familiar - Oeste":6998,"Agricultura Familiar - Sul":1216,"Assentamento/Agrovila - Central Adjacente I":106,"Assentamento/Agrovila - Leste":1213,"Assentamento/Agrovila - Norte":318,"Assentamento/Agrovila - Oeste":465,"Assentamento/Agrovila - Sul":351,"Características Urbanas - Central Adjacente I":218,"Características Urbanas - Central Adjacente II":5452,"Características Urbanas - Leste":9664,"Características Urbanas - Norte":5556,"Características Urbanas - Oeste":5780,"Características Urbanas - Sul":5009}
n_h=dom.groupby(EST)[ID].count(); dom['w']=dom[EST].map(N_H)/dom[EST].map(n_h)
print('soma pesos',dom.w.sum().round(0))
mor['w']=mor[ID].map(dom.set_index(ID).w)
seg=pd.read_csv('segmentos_completo.csv')[['a01ra','a01subpopulacao','Cluster_K3','ra_norm']]
dom=dom.merge(seg,on=[RA,SUB],how='left'); mor=mor.merge(seg,on=[RA,SUB],how='left')
print('dom fora dos 49 seg:',dom.Cluster_K3.isna().sum(),' mor:',mor.Cluster_K3.isna().sum())
def sn(s):  # 1 sim, 0 não, nan caso contrário
    s=s.astype('string').str.strip()
    return s.map({'Sim':1.0,'Não':0.0}).astype(float)
dom['incentivo']=sn(dom['B22_1_a_sua_propriedade_ja_recebeu_incentivos_a_producao'])
dom['outorga']=sn(dom['B14_9_possui_outorga_para_uso_da_agua'])
q=dom['B22_2_a_sua_propriedade_ja_recebeu_incentivos_a_producao_quais'].astype('string')
inc=dom.incentivo==1
for k,pat in {'Pronaf':'Pronaf','ATER':'ATER','Fomento':'Fomento','FCO':'FCO'}.items():
    dom[k]=np.where(inc,q.str.contains(pat,na=False).astype(float),np.nan)
mor['beneficio']=sn(mor['I24_2__Recebe_ou_já_recebeu_benefícios_de_programas_sociais'])
mor['pediu_auxilio']=sn(mor['F06_7__Pedido_de_auxílio_benefícios_governamentais'])
mor['cnpj']=sn(mor['I16__Possui_CNPJ'])
mor['pix']=sn(mor['F06_6__Utilizou_a_internet_para_transações_financeiras_comerciais_PIX'])
mor['curso_agr']=sn(mor['H15_2__Curso_em_andamento__Profissional_voltado_para_mercado_agrícola'])
DV=['incentivo','Pronaf','ATER','Fomento','FCO','outorga']; MV=['beneficio','pediu_auxilio','cnpj','pix','curso_agr']
lab={'incentivo':'Incentivo à produção (já recebeu)','Pronaf':'Pronaf (entre os que receberam incentivo)','ATER':'ATER (entre os que receberam incentivo)','Fomento':'Fomento (entre os que receberam incentivo)','FCO':'FCO (entre os que receberam incentivo)','outorga':'Outorga de uso da água','beneficio':'Benefício de programa social (recebe ou já recebeu)','pediu_auxilio':'Pediu auxílio/benefício governamental (usuários de internet)','cnpj':'Possui CNPJ (ocupados)','pix':'Usou internet para transações financeiras/PIX','curso_agr':'Curso profissional agrícola em andamento'}
def wilson(p,n):
    if n<=0 or np.isnan(p): return (np.nan,np.nan)
    z=1.96; d=1+z*z/n; c=p+z*z/(2*n); a=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)); return ((c-a)/d,(c+a)/d)
def est(df,v):
    x=df[[v,'w']].dropna()
    if len(x)==0: return dict(p=np.nan,lo=np.nan,hi=np.nan,n=0,k=0)
    p=np.average(x[v],weights=x.w); neff=x.w.sum()**2/(x.w**2).sum(); lo,hi=wilson(p,neff)
    return dict(p=p*100,lo=lo*100,hi=hi*100,n=len(x),k=int(x[v].sum()))
rows=[]
for by in ['a01subpopulacao','Cluster_K3','ALL']:
    for v in DV+MV:
        df=(dom if v in DV else mor).copy(); df['ALL']='Total'
        for g,x in df.groupby(by): rows.append(dict(by=by,grupo=g,var=v,**est(x,v)))
r=pd.DataFrame(rows); r.to_csv('apoio_prop.csv',index=False)
pd.set_option('display.width',220)
for by in ['ALL','a01subpopulacao','Cluster_K3']:
    print(by); print(r[r.by==by].pivot(index='var',columns='grupo',values='p').round(1).loc[DV+MV])
print(r[r.by=='ALL'][['var','p','lo','hi','n','k']].round(1))
print('n por subpop/cluster p/ incentivo'); print(r[(r['var']=='incentivo')&(r.by!='ALL')][['by','grupo','n','k']])
# KW segmento-nível
sg=[]
for v in DV+MV:
    df=dom if v in DV else mor
    t=df.groupby(['a01ra','a01subpopulacao','Cluster_K3']).apply(lambda x: pd.Series(dict(p=np.average(x[v].dropna(),weights=x.w[x[v].notna()]) if x[v].notna().sum()>=5 else np.nan))).reset_index()
    t=t.dropna(); 
    if v in ('Fomento','FCO'): continue
    gs=[g.p.values for _,g in t.groupby('Cluster_K3') if len(g)>0]
    if len(gs)>1:
        H,p=kruskal(*gs); sg.append(dict(var=v,nseg=len(t),H=H,p=p))
print(pd.DataFrame(sg).round(3))
pd.DataFrame(sg).to_csv('apoio_kw.csv',index=False)
# por RA (para mapas)
ra=[]
for v in ['incentivo','beneficio','cnpj','outorga','pix','pediu_auxilio']:
    df=dom if v in DV else mor
    for g,x in df.groupby('ra_norm'):
        e=est(x,v); ra.append(dict(ra=g,var=v,**e))
pd.DataFrame(ra).to_csv('apoio_ra.csv',index=False)
# comp. contagens B22_2
print(dom.loc[inc,['Pronaf','ATER','Fomento','FCO']].sum(), inc.sum())
