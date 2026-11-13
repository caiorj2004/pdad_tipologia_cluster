# -*- coding: utf-8 -*-
"""S6b — Deixa-uma-RA-fora (e deixa-uma-subpopulação-fora) REFAZENDO O PCA.

Diferença para o teste S6 do programa 21: lá o PCA ficava fixo (calculado com os 49 segmentos) e só o K-Means era refeito.
Aqui, para cada RA (e para cada subpopulação) removida:
  1. remove os segmentos dela da matriz de variáveis (Fase 1 do programa 03, já agregada);
  2. refaz a padronização e o PCA de cada dimensão, e REDECIDE o número de componentes pela análise paralela de Monte Carlo (percentil 95);
  3. refaz o K-Means (K da variante, 50 inicializações);
  4. mede o ARI contra a partição de referência, só nos segmentos que restaram.
Rodar dentro de resultados/<variante>/ (lê config.txt e segmentos_completo.csv). Saída: stress_S6b_*.csv.
"""
import os, sys, numpy as np, pandas as pd
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score as ARI

cfg = dict(kv.split('=', 1) for kv in Path('config.txt').read_text().split())
os.environ.update({'V2': cfg['V2'], 'K_FINAL': cfg['K'], 'WINSOR_P': cfg['WINSOR_P'], 'COMP': cfg['COMP']})
os.environ.setdefault('TCC_BASE', str(Path(__file__).resolve().parent.parent / 'dados' / 'gold'))
os.environ.setdefault('TCC_SHP', str(Path(__file__).resolve().parent.parent / 'dados' / 'regioes_administrativas.zip'))
src = open(Path(__file__).resolve().parent / '03_analise_cluster_4d.py', encoding='utf-8').read()
exec(compile(src[:src.index('def monte_carlo_parallel_analysis')], '03_fase1', 'exec'))   # df_matriz_A, DIMENSOES_REF, K_FINAL, CHAVES, RA, SUB

K = K_FINAL; SEED = 42; N_MC = 300
ref = pd.read_csv('segmentos_completo.csv')
chave = df_matriz_A[CHAVES].merge(ref[CHAVES + [f'Cluster_K{K}']], on=CHAVES, how='left')
assert chave[f'Cluster_K{K}'].notna().all() and len(chave) == len(df_matriz_A)
lab_ref = chave[f'Cluster_K{K}'].astype(int).to_numpy()
ncomp_ref = {d: sum(c.startswith(d + '_MC_PC') for c in ref.columns) for d in DIMENSOES_REF}

def mc_p95(n, p, it=N_MC, seed=SEED):
    rs = np.random.RandomState(seed); sim = np.zeros((it, p))
    for i in range(it): sim[i] = np.linalg.eigvalsh(np.corrcoef(rs.normal(size=(n, p)), rowvar=False))[::-1]
    return np.percentile(sim, 95, axis=0)

def pca_refeito(M):
    """Padroniza, faz o PCA por dimensão e decide o nº de componentes por análise paralela (mesma regra do programa 03)."""
    Z, nc = [], {}
    for d, vs in DIMENSOES_REF.items():
        v = [x for x in vs if x in M.columns and M[x].notna().sum() >= 3 and M[x].dropna().nunique() > 1]
        X = StandardScaler().fit_transform(M[v].apply(pd.to_numeric, errors='coerce').fillna(0))
        eig = np.linalg.eigvalsh(np.corrcoef(X, rowvar=False))[::-1]
        k = max(1, int((eig > mc_p95(len(X), len(v))).sum())); nc[d] = k
        Z.append(PCA(n_components=k).fit_transform(X))
    return np.hstack(Z), nc

rows_ra, rows_sub = [], []
for col, rows in ((RA, rows_ra), (SUB, rows_sub)):
    for g in sorted(df_matriz_A[col].unique()):
        m = (df_matriz_A[col] != g).to_numpy()
        if m.sum() < 4 * K: continue
        Z, nc = pca_refeito(df_matriz_A[m].reset_index(drop=True))
        lab = KMeans(K, n_init=50, random_state=SEED).fit(Z).labels_
        rows.append({col: g, 'n_removidos': int((~m).sum()), 'ARI_vs_ref_restantes': ARI(lab_ref[m], lab),
                     'n_componentes_novo': sum(nc.values()), 'n_componentes_ref': sum(ncomp_ref.values()), 'tamanhos': str(sorted(np.bincount(lab, minlength=K).tolist(), reverse=True))})
a = pd.DataFrame(rows_ra); b = pd.DataFrame(rows_sub)
a.to_csv('stress_S6b_deixa_uma_RA_refaz_pca.csv', index=False); b.to_csv('stress_S6b_deixa_uma_subpop_refaz_pca.csv', index=False)
R = {'S6b_ARI_medio_sem_uma_RA': a.ARI_vs_ref_restantes.mean(), 'S6b_ARI_min_sem_uma_RA': a.ARI_vs_ref_restantes.min(), 'S6b_RA_critica': a.loc[a.ARI_vs_ref_restantes.idxmin(), RA],
     'S6b_pct_RAs_ARI<0.5': float((a.ARI_vs_ref_restantes < 0.5).mean()), 'S6b_n_componentes_varia': bool(a.n_componentes_novo.nunique() > 1),
     'S6b_ARI_medio_sem_uma_subpop': b.ARI_vs_ref_restantes.mean(), 'S6b_ARI_min_sem_uma_subpop': b.ARI_vs_ref_restantes.min()}
pd.DataFrame([R]).to_csv('stress_S6b_resumo.csv', index=False)
print(a.round(3).to_string(index=False)); print(b.round(3).to_string(index=False)); print(R)
