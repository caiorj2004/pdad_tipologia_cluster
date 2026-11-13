# -*- coding: utf-8 -*-
# Bateria de critérios (silhueta, CH, DB, gap, cofenética) e estabilidade de K-Means e Ward por subamostragem.
# Entrada: X_pca.npy (programa 03)
# Saída: ward_vs_km_bateria.csv
import numpy as np, pandas as pd, warnings, sys
warnings.filterwarnings('ignore')
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score, adjusted_rand_score
from scipy.cluster.hierarchy import linkage, fcluster, cophenet
from scipy.spatial.distance import pdist
X=np.load('X_pca.npy'); N,p=X.shape
RNG=np.random.default_rng(42)
def km(Z,k,rs=42,ni=10): return KMeans(k,n_init=ni,random_state=rs).fit_predict(Z)
def wd(Z,k,rs=0,ni=0): return AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(Z)
def wss(Z,l): return sum(((Z[l==c]-Z[l==c].mean(0))**2).sum() for c in np.unique(l))
def gap(Z,k,f,nref=50):
    lo,hi=Z.min(0),Z.max(0)
    r=[np.log(wss(U,f(U,k))) for U in (RNG.uniform(lo,hi,size=Z.shape) for _ in range(nref))]
    return np.mean(r)-np.log(wss(Z,f(Z,k)))
cov,mu=np.cov(X,rowvar=False),X.mean(0)
NP=500; NS=200
rows=[]
for name,f,ni in (('K-Means',km,5),('Ward',wd,0)):
    for k in range(2,9):
        lab=km(X,k,42,50) if name=='K-Means' else wd(X,k)
        real=silhouette_score(X,lab)
        nc,nm=[],[]
        for i in range(NP):
            Z=X.copy()
            for c in range(p): RNG.shuffle(Z[:,c])
            nc.append(silhouette_score(Z,f(Z,k,i,ni)))
            Z2=RNG.multivariate_normal(mu,cov,size=N)
            nm.append(silhouette_score(Z2,f(Z2,k,i,ni)))
        nc,nm=np.array(nc),np.array(nm)
        ar=[]
        for b in range(NS):
            idx=RNG.choice(N,int(.8*N),replace=False)
            ar.append(adjusted_rand_score(lab[idx],f(X[idx],k,b,10)))
        rows.append(dict(metodo=name,K=k,silhueta=round(real,3),CH=round(calinski_harabasz_score(X,lab),2),DB=round(davies_bouldin_score(X,lab),2),
            gap=round(gap(X,k,lambda Z,kk: f(Z,kk,42,10)),3),menor=int(np.bincount(lab).min()),
            nulo_emb=round(nc.mean(),3),p_emb=round((np.sum(nc>=real)+1)/(NP+1),3),
            nulo_cov=round(nm.mean(),3),p_cov=round((np.sum(nm>=real)+1)/(NP+1),3),estab=round(np.mean(ar),3)))
        print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv('ward_vs_km_bateria.csv',index=False)
