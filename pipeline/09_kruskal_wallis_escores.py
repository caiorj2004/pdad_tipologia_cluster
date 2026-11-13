# -*- coding: utf-8 -*-
# Kruskal-Wallis das 34 variáveis e dos 6 componentes entre os 5 grupos, com correção FDR e tamanho de efeito.
# Entrada: segmentos_completo.csv (programa 03)
# Saída: kw_todas.csv, escores_medios_grupo.csv
import numpy as np, pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
d=pd.read_csv('segmentos_completo.csv'); CL='Cluster_K3'
V={'D1':["log_renda_domiciliar_pc","possui_automovel","possui_motocicleta","possui_geladeira","possui_freezer","possui_lavadora","possui_secadora","possui_tv","possui_ar_condicionado","possui_microondas"],
'D2':["escolaridade_adulto","trabalhou_30d","tem_ctps","contribui_previdencia","ivs_df"],
'D3':["comodos_por_pessoa","banheiros_por_pessoa","agua_rede","agua_poco_artesiano","agua_poco_cisterna","captacao_chuva","caixa_agua","solucao_sanitaria_estruturada","acesso_energetico_estruturado","destino_com_coleta_direta","rua_pavimentada","iluminacao_publica","drenagem","ponto_onibus","qtd_estabelecimentos"],
'D4':["internet_domiciliar","internet_banda_larga_fixa","internet_movel","possui_computador"]}
comp=['D1_Bens_Duraveis_MC_PC1','D2_Humano_Trabalho_MC_PC1','D3_Fisico_Habitacional_Territorial_MC_PC1','D3_Fisico_Habitacional_Territorial_MC_PC2','D4_Tecnologico_Conectividade_MC_PC1']
rows=[]
n=len(d)
for dim,vs in V.items():
    for v in vs:
        g=[x[v].dropna().values for _,x in d.groupby(CL)]
        H,p=kruskal(*g); k=len(g)
        rows.append(dict(tipo='var',dim=dim,var=v,H=H,p=p,eps2=max(0,(H-k+1)/(n-k))))
for c in comp:
    g=[x[c].values for _,x in d.groupby(CL)]
    H,p=kruskal(*g); k=len(g)
    rows.append(dict(tipo='comp',dim=c[:2],var=c,H=H,p=p,eps2=max(0,(H-k+1)/(n-k))))
r=pd.DataFrame(rows)
r['q']=multipletests(r.p,method='fdr_bh')[1]
r.to_csv('kw_todas.csv',index=False)
print(r.round(4).to_string())
print('sig p<.05:',(r.p<.05).sum(),'de',len(r),' q<.05:',(r.q<.05).sum(), ' vars:', ((r.tipo=='var')&(r.q<.05)).sum(),'de 34; comps q<.05:',((r.tipo=='comp')&(r.q<.05)).sum())
m=d.groupby(CL)[comp].mean(); m['n']=d.groupby(CL).size(); print(m.round(2)); m.to_csv('escores_medios_grupo.csv')
sd=d[comp].std(); print((d.groupby(CL)[comp].mean()/sd).round(2))
