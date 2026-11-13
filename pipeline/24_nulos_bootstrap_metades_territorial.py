# -*- coding: utf-8 -*-
"""Três nulos que faltavam para interpretar a estabilidade e a "homogeneidade territorial" (rodar em resultados/<variante>/ depois do 16).
 (1) NULO DO BOOTSTRAP: o ARI do bootstrap só tem leitura se houver um valor de comparação para dados SEM estrutura de grupos. Gera matrizes de 49 segmentos com a mesma
     covariância entre variáveis da matriz observada (descontada a variância amostral estimada pelo próprio bootstrap), soma o ruído amostral de cada segmento (vetores de desvio das réplicas em boot_matrizes.npy, que preservam a correlação do ruído dentro do segmento), roda a mesma cadeia (PCA por dimensão com nº de componentes fixo + K-Means) e mede o ARI do "bootstrap" paramétrico contra a partição dessa matriz nula.
 (2) REPLICABILIDADE POR METADES: cada segmento tem seus domicílios divididos ao acaso em duas metades; PCA e K-Means em cada metade; ARI entre as duas partições. Comparado ao mesmo cálculo no nulo (ruído com variância dobrada).
 (3) NULO DE ATRIBUIÇÃO DOMICILIAR: dentro de cada estrato amostral, os domicílios são redistribuídos ao acaso entre os segmentos do estrato (mantendo o tamanho de cada segmento). Isso destrói o vínculo domicílio->RA
     mas mantém estrato e subpopulação; mede quanto da AMI(RA) e AMI(subpopulação) observada seria esperada só pelo desenho amostral e pelo tamanho dos segmentos.
Env: N_NULO (padrão 100), N_META (100), N_PERM (100). Saída: nulos_24_resumo.csv, saida em stdout."""
import os, numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score as ARI, adjusted_mutual_info_score as AMI

SEED = 42; N_NULO = int(os.environ.get('N_NULO', 100)); N_META = int(os.environ.get('N_META', 100)); N_PERM = int(os.environ.get('N_PERM', 100)); NI = 10
cfg = dict(kv.split('=', 1) for kv in Path('config.txt').read_text().split())
os.environ.update({'V2': cfg['V2'], 'K_FINAL': cfg['K'], 'WINSOR_P': cfg['WINSOR_P'], 'COMP': cfg['COMP']})
os.environ.setdefault('TCC_BASE', str(Path(__file__).resolve().parent.parent / 'dados' / 'gold')); os.environ.setdefault('TCC_SHP', str(Path(__file__).resolve().parent.parent / 'dados' / 'regioes_administrativas.zip'))
src = open(Path(__file__).resolve().parent / '03_analise_cluster_4d.py', encoding='utf-8').read()
exec(compile(src[:src.index('def monte_carlo_parallel_analysis')], '03_fase1', 'exec'))
K = K_FINAL; ref = pd.read_csv('segmentos_completo.csv')
lab_ref = df_matriz_A[CHAVES].merge(ref[CHAVES + [CL]], on=CHAVES, how='left')[CL].astype(int).to_numpy()
ncomp = {d: sum(c.startswith(d + '_MC_PC') for c in ref.columns) for d in DIMENSOES_REF}
cols = [v for vs in DIMENSOES_REF.values() for v in vs]; pos = {d: [cols.index(v) for v in vs] for d, vs in DIMENSOES_REF.items()}
X_obs = df_matriz_A[cols].apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(float)
def cadeia(X, seed, ni=NI):
    Z = np.hstack([PCA(n_components=ncomp[d]).fit_transform(StandardScaler().fit_transform(X[:, pos[d]])) for d in DIMENSOES_REF])
    return KMeans(K, n_init=ni, random_state=seed).fit(Z).labels_
BM = np.load('boot_matrizes.npy'); assert BM.shape[2] == len(cols)
DEV = np.nan_to_num(BM - BM.mean(axis=0)); sd = np.nan_to_num(BM.std(axis=0, ddof=1)); rng = np.random.default_rng(SEED)
S = np.cov(X_obs, rowvar=False) - np.diag((sd ** 2).mean(axis=0))
w, V = np.linalg.eigh(S); S = (V * np.clip(w, 1e-8, None)) @ V.T; mu = X_obs.mean(axis=0)
def ruido(fator_var=1.0):
    """Ruído amostral de cada segmento: sorteia, por segmento, o vetor de desvios de uma réplica do bootstrap (preserva a correlação entre variáveis dentro do segmento)."""
    b = rng.integers(0, DEV.shape[0], X_obs.shape[0]); return np.sqrt(fator_var) * DEV[b, np.arange(X_obs.shape[0])]
def nulo(fator_var=1.0):
    T = rng.multivariate_normal(mu, S, size=X_obs.shape[0]); return T, fator_var
R = {}
# (1) nulo do bootstrap
obs_boot = np.load('boot_ari.npy'); R['boot_ARI_obs_medio'] = float(obs_boot.mean())
m_null = []
for i in range(N_NULO):
    T, s = nulo(); X0 = T + ruido(s); l0 = cadeia(X0, i)
    m_null.append(np.mean([ARI(l0, cadeia(X0 + ruido(s), 1000 + j)) for j in range(15)]))
m_null = np.array(m_null); R.update(boot_ARI_nulo_medio=float(m_null.mean()), boot_ARI_nulo_p95=float(np.percentile(m_null, 95)), boot_ARI_p_vs_nulo=float((1 + (m_null >= obs_boot.mean()).sum()) / (1 + N_NULO)))
print(f'(1) NULO DO BOOTSTRAP ({N_NULO} matrizes nulas x 15 réplicas): ARI observado médio {obs_boot.mean():.3f} | nulo médio {m_null.mean():.3f} (p95 {np.percentile(m_null,95):.3f}) | p = {R["boot_ARI_p_vs_nulo"]:.3f}')
# (2) metades
seg_id = df_dom.groupby(CHAVES, dropna=False).ngroup().to_numpy(); hh_pos = {h: i for i, h in enumerate(df_dom[ID])}; mor_hh = df_mor[ID].map(hh_pos).to_numpy()
df_dom_s = df_dom[CHAVES + VAR_DOM + ['peso_dom_corrigido']].copy(); df_mor_s = df_mor[CHAVES + VAR_MOR + [PESO_MOR]].copy()
def matriz_de(m):
    d = df_dom_s.copy(); d['w_b'] = d['peso_dom_corrigido'] * m; o = df_mor_s.copy(); o['w_b'] = o[PESO_MOR] * m[mor_hh]
    md = media_ponderada_vetorizada(d, CHAVES, VAR_DOM, 'w_b'); mo = media_ponderada_vetorizada(o, CHAVES, VAR_MOR, 'w_b')
    md['log_renda_domiciliar_pc'] = np.where(md['renda_domiciliar_pc'] >= 0, np.log1p(md['renda_domiciliar_pc']), np.nan)
    M = derivar_variaveis_modelo(df_matriz_A[CHAVES + ['amostra_n']].merge(md, on=CHAVES, how='left').merge(mo, on=CHAVES, how='left'))
    return M
idx_seg = [np.where(seg_id == s)[0] for s in range(seg_id.max() + 1)]
ari_meta = []
for b in range(N_META):
    ma = np.zeros(len(df_dom)); mb = np.zeros(len(df_dom))
    for ix in idx_seg:
        p = rng.permutation(ix); h = len(p) // 2; ma[p[:h]] = 1; mb[p[h:]] = 1
    A = matriz_de(ma)[cols].apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(float); Bm = matriz_de(mb)[cols].apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(float)
    ari_meta.append(ARI(cadeia(A, b), cadeia(Bm, 5000 + b)))
ari_meta = np.array(ari_meta); n_meta = []
for i in range(N_META):
    T, s = nulo(2.0); n_meta.append(ARI(cadeia(T + ruido(s), i), cadeia(T + ruido(s), 9000 + i)))
n_meta = np.array(n_meta)
R.update(metades_ARI_obs_medio=float(ari_meta.mean()), metades_ARI_obs_p5=float(np.percentile(ari_meta, 5)), metades_ARI_nulo_medio=float(n_meta.mean()), metades_ARI_nulo_p95=float(np.percentile(n_meta, 95)))
print(f'(2) METADES ({N_META} divisões): ARI entre metades observado médio {ari_meta.mean():.3f} (p5 {np.percentile(ari_meta,5):.3f}) | nulo médio {n_meta.mean():.3f} (p95 {np.percentile(n_meta,95):.3f})')
# (3) atribuição domiciliar
ra_seg = df_matriz_A[RA].to_numpy(); sub_seg = df_matriz_A[SUB].to_numpy()
seg_idx_matriz = df_matriz_A[CHAVES].reset_index(drop=True).merge(df_dom.groupby(CHAVES, dropna=False).ngroup().rename('sid').to_frame().join(df_dom[CHAVES]).drop_duplicates(), on=CHAVES, how='left')['sid'].to_numpy()
est = df_dom[ESTRATO].to_numpy(); lab_obs = lab_ref
R['AMI_RA_obs'] = float(AMI(ra_seg, lab_obs)); R['AMI_subpop_obs'] = float(AMI(sub_seg, lab_obs))
a_ra, a_sub, a_ari = [], [], []
for b in range(N_PERM):
    novo = seg_id.copy()
    for e in np.unique(est):
        ix = np.where(est == e)[0]; novo[ix] = rng.permutation(seg_id[ix])
    d = df_dom_s.copy(); inv = {s: i for i, s in enumerate(seg_idx_matriz)}
    # reconstrói chaves (RA, subpop) a partir do segmento sorteado
    chaves_seg = df_dom[CHAVES].iloc[[np.where(seg_id == s)[0][0] for s in range(seg_id.max() + 1)]].reset_index(drop=True)
    d[CHAVES] = chaves_seg.iloc[novo].reset_index(drop=True).to_numpy(); o = df_mor_s.copy(); o[CHAVES] = d[CHAVES].to_numpy()[mor_hh]
    md = media_ponderada_vetorizada(d, CHAVES, VAR_DOM, 'peso_dom_corrigido'); mo = media_ponderada_vetorizada(o, CHAVES, VAR_MOR, PESO_MOR)
    md['log_renda_domiciliar_pc'] = np.where(md['renda_domiciliar_pc'] >= 0, np.log1p(md['renda_domiciliar_pc']), np.nan)
    M = derivar_variaveis_modelo(df_matriz_A[CHAVES + ['amostra_n']].merge(md, on=CHAVES, how='left').merge(mo, on=CHAVES, how='left'))
    Xp = M[cols].apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(float); lp = cadeia(Xp, b)
    a_ra.append(AMI(ra_seg, lp)); a_sub.append(AMI(sub_seg, lp)); a_ari.append(ARI(lab_ref, lp))
a_ra, a_sub = np.array(a_ra), np.array(a_sub)
R.update(AMI_RA_nulo_atrib_medio=float(a_ra.mean()), AMI_RA_nulo_atrib_p95=float(np.percentile(a_ra, 95)), AMI_RA_p_atrib=float((1 + (a_ra >= R['AMI_RA_obs']).sum()) / (1 + N_PERM)),
         AMI_subpop_nulo_atrib_medio=float(a_sub.mean()), AMI_subpop_p_atrib=float((1 + (a_sub >= R['AMI_subpop_obs']).sum()) / (1 + N_PERM)))
print(f'(3) ATRIBUIÇÃO DOMICILIAR ALEATÓRIA DENTRO DO ESTRATO ({N_PERM} sorteios): AMI(RA) obs {R["AMI_RA_obs"]:.3f} | nulo médio {a_ra.mean():.3f} (p95 {np.percentile(a_ra,95):.3f}) | p = {R["AMI_RA_p_atrib"]:.3f}')
print(f'    AMI(subpop) obs {R["AMI_subpop_obs"]:.3f} | nulo médio {a_sub.mean():.3f} | p = {R["AMI_subpop_p_atrib"]:.3f}')
pd.DataFrame([R]).to_csv('nulos_24_resumo.csv', index=False)
