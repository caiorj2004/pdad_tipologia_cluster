# -*- coding: utf-8 -*-
# Comparação entre K-Means, Ward, GMM (três covariâncias) e DPGMM sobre os 6 componentes, de K=2 a 8.
# Entrada: X_pca.npy (programa 03)
# Saída: comparacao_metodos.csv
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score, adjusted_rand_score
X=np.load('X_pca.npy'); n=len(X)
def km(X,k,seed=42): return KMeans(k,n_init=50,random_state=seed).fit_predict(X)
def ward(X,k): return AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(X)
def gmm(ct):
    return lambda X,k: GaussianMixture(k,covariance_type=ct,n_init=20,random_state=42,reg_covar=1e-3).fit(X).predict(X)
def dp(ct):
    def f(X,k):
        m=BayesianGaussianMixture(n_components=k,weight_concentration_prior_type='dirichlet_process',covariance_type=ct,n_init=20,random_state=42,reg_covar=1e-3,max_iter=500).fit(X)
        return m.predict(X)
    return f
M={'K-Means':km,'Ward':ward,'GMM esférica':gmm('spherical'),'GMM diagonal':gmm('diag'),'GMM completa':gmm('full'),'DPGMM diagonal':dp('diag')}
rows=[];ref={}
for k in range(2,9): ref[k]=km(X,k)
rng=np.random.default_rng(42)
subs=[rng.choice(n,int(round(.8*n)),replace=False) for _ in range(200)]
for name,f in M.items():
    for k in range(2,9):
        try:
            l=f(X,k)
        except Exception as e:
            continue
        nc=len(set(l))
        sil=silhouette_score(X,l) if 1<nc<n else np.nan
        ch=calinski_harabasz_score(X,l) if nc>1 else np.nan
        db=davies_bouldin_score(X,l) if nc>1 else np.nan
        mn=np.bincount(l).min() if True else 0
        mn=min(c for c in np.bincount(l) if c>0)
        ari_ref=adjusted_rand_score(ref[k],l)
        aris=[]
        for s in subs:
            try:
                ls=f(X[s],k); aris.append(adjusted_rand_score(l[s],ls))
            except Exception: pass
        rows.append(dict(metodo=name,K=k,grupos_obs=nc,silhueta=round(sil,3),CH=round(ch,2),DB=round(db,2),menor=mn,ARI_vs_KMeans=round(ari_ref,2),estab_ARI=round(np.mean(aris),3)))
df=pd.DataFrame(rows);df.to_csv('comparacao_metodos.csv',index=False)
pd.set_option('display.width',200)
print(df.to_string(index=False))
