# -*- coding: utf-8 -*-
# Pareia os rótulos de Ward aos do K-Means (algoritmo húngaro), cruza com a classe núcleo/fronteira e desenha o dendrograma.
# Entrada: X_pca.npy, seg_labels.csv, figuras/coassociacao_segmentos.csv (programa 03)
# Saída: km_ward_labels.csv, figuras/dendrograma_ward.png
import numpy as np, pandas as pd, re, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.optimize import linear_sum_assignment
X=np.load('X_pca.npy'); N=len(X)
seg=pd.read_csv('seg_labels.csv')
seg['sub']=seg.a01subpopulacao.str.replace('Características Urbanas','Moradia e Lazer').str.replace('Agricultura ','Agr. ').str.replace('Assentamento/Agrovila','Assentamento')
seg['rot']=seg.a01ra+' — '+seg['sub']
km=KMeans(5,n_init=50,random_state=42).fit_predict(X)
wd=AgglomerativeClustering(5,linkage='ward').fit_predict(X)
ct=pd.crosstab(km,wd); r,c=linear_sum_assignment(-ct.values)
mp={ct.columns[cc]:ct.index[rr] for rr,cc in zip(r,c)}
wd2=np.array([mp[x] for x in wd])
co=pd.read_csv('figuras/coassociacao_segmentos.csv')
print(co.columns.tolist())
out=seg[['a01ra','sub','rot','amostra_n']].copy(); out['km']=km; out['ward']=wd2
out=out.merge(co[['Segmento','Classe']],left_on='rot',right_on='Segmento',how='left')
out.to_csv('km_ward_labels.csv',index=False)
print(pd.crosstab(out.km,out.ward)); print(out[out.km!=out.ward][['rot','amostra_n','km','ward','Classe']])
print(pd.crosstab(out.Classe,out.km==out.ward))
# dendrogram
Z=linkage(X,'ward')
cols=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e']
fig,ax=plt.subplots(figsize=(8,9))
thr=(Z[-4,2]+Z[-5,2])/2
d=dendrogram(Z,orientation='right',labels=out.rot.values,leaf_font_size=7,link_color_func=lambda k:'#555555',ax=ax)
ax.axvline(thr,ls='--',c='k',lw=.8)
lab=dict(zip(out.rot,out.km))
for t in ax.get_ymajorticklabels():
    t.set_color(cols[lab[t.get_text()]])
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=cols[i],label=f'Grupo {i} (K-Means)') for i in range(5)],loc='lower right',fontsize=8,frameon=False,bbox_to_anchor=(1,0.18))
ax.set_xlabel('Distância de Ward (altura de fusão)'); 
ax.set_title('Dendrograma de Ward (49 segmentos)\nCor do rótulo = grupo do K-Means (K=5); linha = corte em 5 grupos',fontsize=9)
for s in ['top','right']: ax.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig('figuras/dendrograma_ward.png',dpi=200)
h=Z[:,2]; print([ (k,round(h[-(k-1)],2)) for k in range(2,9)])
