# -*- coding: utf-8 -*-
# Estabilidade de K-Means e Ward por subamostragem de 80% (200 repetições), inclusive cruzada.
# Entrada: X_pca.npy (programa 03)
# Saída: impressão no terminal
import numpy as np, warnings
warnings.filterwarnings('ignore')
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score
X=np.load('X_pca.npy');N=len(X);rng=np.random.default_rng(11)
lab=KMeans(5,n_init=50,random_state=42).fit_predict(X); wl=AgglomerativeClustering(5,linkage='ward').fit_predict(X)
subs=[rng.choice(N,int(.8*N),replace=False) for _ in range(200)]
for ni in (1,10,50):
    a=[adjusted_rand_score(lab[s],KMeans(5,n_init=ni,random_state=b).fit_predict(X[s])) for b,s in enumerate(subs)]
    print('KM n_init',ni,round(np.mean(a),3),'sd',round(np.std(a),3))
a=[adjusted_rand_score(wl[s],AgglomerativeClustering(5,linkage='ward').fit_predict(X[s])) for s in subs]
print('Ward',round(np.mean(a),3),'sd',round(np.std(a),3))
# estabilidade cruzada: partição de referência de um método vs refit do outro na subamostra
a=[adjusted_rand_score(lab[s],AgglomerativeClustering(5,linkage='ward').fit_predict(X[s])) for s in subs]
print('ref KM vs Ward refit',round(np.mean(a),3))
a=[adjusted_rand_score(wl[s],KMeans(5,n_init=50,random_state=b).fit_predict(X[s])) for b,s in enumerate(subs)]
print('ref Ward vs KM refit',round(np.mean(a),3))
