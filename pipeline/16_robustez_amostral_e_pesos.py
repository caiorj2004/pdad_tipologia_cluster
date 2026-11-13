# -*- coding: utf-8 -*-
# Robustez à amostragem de domicílios dentro dos segmentos e ao peso de cada segmento.
#   (a) bootstrap por domicílios dentro de cada segmento: refaz médias ponderadas, variáveis derivadas (ilr, renda), PCA por dimensão e K-Means (K=K_FINAL), com perfis e ICs dos grupos (alinhamento húngaro);
#   (c) ARI de estabilidade contra o nulo de covariância; (d) encolhimento empírico-Bayes das médias dos segmentos pequenos; (e) regra de 'não classificável'
#   (b) K-Means ponderado pelo número de domicílios expandidos de cada segmento (peso integral e raiz quadrada).
# Entrada: bases corrigidas (programa 01), Base_dados-IVS.xlsx e segmentos_completo.csv (programa 03; usa a Fase 1 do 03 para reconstruir as variáveis)
# Saída: robustez_amostral_segmentos.csv, robustez_pesos_segmentos.csv, saida_16.txt
import os, re, sys, numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score, silhouette_score

B = int(os.environ.get('N_BOOT', 500)); SEED = 42
NINIT = int(os.environ.get('NINIT_BOOT', 10))   # inicializações do K-Means em cada réplica (10 = rápido; 50 = igual à referência, separa ótimo local de incerteza amostral)
src = open(os.environ.get('TCC_03', '03_analise_cluster_4d.py'), encoding='utf-8').read()
# reaproveita configuração e Fase 1 do programa 03 (variáveis, pesos e matriz de 49 segmentos), sem rodar as etapas seguintes
exec(compile(src[:src.index('def monte_carlo_parallel_analysis')], '03_fase1', 'exec'))

ref = pd.read_csv('segmentos_completo.csv')
chave = df_matriz_A[CHAVES].merge(ref[CHAVES + [CL]], on=CHAVES, how='left')
assert chave[CL].notna().all() and len(chave) == 49
lab_ref = chave[CL].astype(int).to_numpy()
K = K_FINAL
ncomp = {d: sum(c.startswith(d + '_MC_PC') for c in ref.columns) for d in DIMENSOES_REF}
print('componentes por dimensão:', ncomp, '| total', sum(ncomp.values()))

def pca_matriz(M):
    """M: DataFrame com as 49 linhas e as variáveis; devolve a matriz de escores (mesma regra do programa 03)."""
    Z = []
    for d, vs in DIMENSOES_REF.items():
        X = StandardScaler().fit_transform(M[vs].apply(pd.to_numeric, errors='coerce').fillna(0))
        Z.append(PCA(n_components=ncomp[d]).fit_transform(X))
    return np.hstack(Z)

# confere a reprodução da partição de referência
M0 = df_matriz_A.copy()
Z0 = pca_matriz(M0)
lab0 = KMeans(K, n_init=50, random_state=SEED).fit(Z0).labels_
print('ARI reprodução da referência:', round(adjusted_rand_score(lab0, lab_ref), 3))

# ---------------- (a) bootstrap por domicílios dentro de cada segmento
rng = np.random.default_rng(SEED)
seg_id = df_dom.groupby(CHAVES, dropna=False).ngroup().to_numpy()
hh_pos = {h: i for i, h in enumerate(df_dom[ID])}
mor_hh = df_mor[ID].map(hh_pos).to_numpy()
idx_por_seg = [np.where(seg_id == s)[0] for s in range(seg_id.max() + 1)]
df_dom_s = df_dom[CHAVES + VAR_DOM + ['peso_dom_corrigido']].copy(); df_mor_s = df_mor[CHAVES + VAR_MOR + [PESO_MOR]].copy()   # só as colunas usadas
coassoc = np.zeros((49, 49)); aris = []; sil = []; Zs = []
PERF = ['renda_domiciliar_pc', 'escolaridade_adulto', 'internet_banda_larga_fixa', 'possui_computador', 'agua_rede', 'solucao_sanitaria_estruturada', 'rua_pavimentada', 'trabalhou_30d', 'contribui_previdencia', 'possui_automovel', 'ivs_df']
PERF = [p for p in PERF if p in df_matriz_A.columns]
def perfil_grupos(M, lab, w, K):
    out = np.full((K, len(PERF)), np.nan)
    for g in range(K):
        m = lab == g
        if m.any(): out[g] = [np.average(M.loc[m, p].fillna(0), weights=w[m]) for p in PERF]
    return out
perfs = []; ARI_BOOT_N = []
for b in range(B):
    m = np.zeros(len(df_dom))
    for idx in idx_por_seg:
        m[idx] = np.bincount(rng.integers(0, len(idx), len(idx)), minlength=len(idx))
    d = df_dom_s.copy(); d['w_b'] = d['peso_dom_corrigido'] * m
    o = df_mor_s.copy(); o['w_b'] = o[PESO_MOR] * m[mor_hh]
    md = media_ponderada_vetorizada(d, CHAVES, VAR_DOM, 'w_b'); mo = media_ponderada_vetorizada(o, CHAVES, VAR_MOR, 'w_b')
    md['log_renda_domiciliar_pc'] = np.where(md['renda_domiciliar_pc'] >= 0, np.log1p(md['renda_domiciliar_pc']), np.nan)
    M = df_matriz_A[CHAVES + ['amostra_n']].merge(md, on=CHAVES, how='left').merge(mo, on=CHAVES, how='left')
    M = derivar_variaveis_modelo(M)
    Z = pca_matriz(M); Zs.append(Z)
    lab = KMeans(K, n_init=NINIT, random_state=b).fit(Z).labels_
    ct = pd.crosstab(lab_ref, lab).reindex(index=range(K), columns=range(K), fill_value=0).to_numpy()
    r_, c_ = linear_sum_assignment(-ct); mp = {c: r for r, c in zip(r_, c_)}
    lab_al = np.array([mp[l] for l in lab])                      # alinhamento húngaro ao grupo de referência
    for pn in PERF:
        if pn not in M.columns: M[pn] = df_matriz_A[pn].to_numpy()   # variáveis só descritivas (ex. ivs_df) vêm da matriz original
    wb = d.groupby(CHAVES, sort=False)['w_b'].sum().reindex(pd.MultiIndex.from_frame(df_matriz_A[CHAVES])).to_numpy()
    perfs.append(perfil_grupos(M, lab_al, np.where(np.isnan(wb), 0, wb) + 1e-9, K))
    aris.append(adjusted_rand_score(lab_ref, lab)); sil.append(silhouette_score(Z, lab))
    coassoc += (lab[:, None] == lab[None, :])
coassoc /= B
res = []
for i in range(49):
    mates = [j for j in range(49) if lab_ref[j] == lab_ref[i] and j != i]
    outros = [j for j in range(49) if lab_ref[j] != lab_ref[i]]
    res.append({'segmento': f"{df_matriz_A.loc[i, RA]} | {df_matriz_A.loc[i, SUB]}", 'grupo_ref': lab_ref[i], 'amostra_n': int(df_matriz_A.loc[i, 'amostra_n']),
                'junto_do_proprio_grupo': coassoc[i, mates].mean(), 'junto_dos_outros_grupos': coassoc[i, outros].mean()})
res = pd.DataFrame(res)
_cf = Path(os.environ.get('TCC_FIG', 'figuras')) / 'coassociacao_segmentos.csv'     # núcleo/fronteira (programa 04 ou 07)
if _cf.exists():
    _mp = {'Agricultura Familiar': 'Agr. Familiar', 'Agricultura Empresarial': 'Agr. Empresarial', 'Assentamento/Agrovila': 'Assentamento', 'Características Urbanas': 'Moradia e Lazer'}
    res['Segmento'] = res['segmento'].map(lambda t: f"{t.split(' | ')[0]} — {_mp[t.split(' | ')[1]]}")
    res = res.merge(pd.read_csv(_cf)[['Segmento', 'Classe', 'Escore_subamostras']], on='Segmento', how='left')
    print('\nColigação no bootstrap por classe (núcleo/fronteira):'); print(res.groupby('Classe')['junto_do_proprio_grupo'].agg(['mean', 'median', 'min', 'max', 'count']).round(2).to_string())
    print('Spearman com o escore de co-associação por subamostras:', round(res[['junto_do_proprio_grupo', 'Escore_subamostras']].corr(method='spearman').iloc[0, 1], 2))
res.to_csv('robustez_amostral_segmentos.csv', index=False)
aris = np.array(aris)
print(f'\n(a) BOOTSTRAP POR DOMICÍLIOS ({B} réplicas, K={K}, n_init={NINIT})')
print(f'ARI vs. partição de referência: média {aris.mean():.3f} | mediana {np.median(aris):.3f} | p5 {np.percentile(aris,5):.3f} | p95 {np.percentile(aris,95):.3f} | min {aris.min():.3f}')
print(f'Silhueta nas réplicas: média {np.mean(sil):.3f} (referência {silhouette_score(Z0, lab_ref):.3f})')
print('\nColigação média com o próprio grupo, por grupo (e por tamanho de amostra):')
print(res.groupby('grupo_ref')['junto_do_proprio_grupo'].agg(['mean', 'min', 'max']).round(2).to_string())
res['n_cat'] = pd.cut(res.amostra_n, [0, 9, 14, 99], labels=['5-9', '10-14', '15+'])
print(res.groupby('n_cat')['junto_do_proprio_grupo'].agg(['mean', 'min', 'count']).round(2).to_string())
print('\nSegmentos com coligação < 0,60:')
print(res[res.junto_do_proprio_grupo < 0.60].sort_values('junto_do_proprio_grupo')[['segmento', 'grupo_ref', 'amostra_n', 'junto_do_proprio_grupo']].round(2).to_string(index=False))

# estabilidade por K nas mesmas réplicas (contexto para o valor de K)
print('\nARI médio do bootstrap por K (referência de cada K = partição de 49 segmentos com n_init=50):')
for k in (2, 3, 4, 5, 6):
    ref_k = KMeans(k, n_init=50, random_state=SEED).fit(Z0).labels_
    a_k = [adjusted_rand_score(ref_k, KMeans(k, n_init=NINIT, random_state=b).fit(Z).labels_) for b, Z in enumerate(Zs[:200])]
    print(f'  K={k}: ARI médio {np.mean(a_k):.3f} (200 réplicas)')

# ---------------- (b) K-Means ponderado por domicílios expandidos
wdom = df_matriz_A['domicilios_expandidos'].to_numpy(float)
linhas = []
print('\n(b) K-MEANS PONDERADO (mesmos 6 escores, K={K}, n_init=50)')
for nome, w in [('sem peso (referência)', np.ones(49)), ('raiz do peso', np.sqrt(wdom)), ('peso integral', wdom)]:
    lab = KMeans(K, n_init=50, random_state=SEED).fit(Z0, sample_weight=w / w.mean()).labels_
    ari = adjusted_rand_score(lab_ref, lab)
    tam = sorted(np.bincount(lab).tolist(), reverse=True)
    # grupo ponderado correspondente a cada grupo de referência (maior sobreposição)
    ct = pd.crosstab(lab_ref, lab); corresp = ct.idxmax(axis=1)
    muda = [i for i in range(49) if corresp[lab_ref[i]] != lab[i]]
    print(f'- {nome}: ARI {ari:.3f}; tamanhos {tam}; segmentos que mudam de grupo: {len(muda)}')
    for i in muda: print('    ', df_matriz_A.loc[i, RA], '|', df_matriz_A.loc[i, SUB], f'(ref G{lab_ref[i]}, domicílios exp. {wdom[i]:.0f})')
    linhas.append({'variante': nome, 'ARI': ari, 'tamanhos': str(tam), 'n_mudam': len(muda)})
    if nome != 'sem peso (referência)':
        pd.DataFrame({'segmento': res.segmento, 'grupo_ref': lab_ref, 'grupo_pond': lab}).to_csv(f"robustez_pesos_{'raiz' if 'raiz' in nome else 'integral'}.csv", index=False)
pd.DataFrame(linhas).to_csv('robustez_pesos_segmentos.csv', index=False)
print('\nDomicílios expandidos por grupo de referência (%):')
print((pd.Series(wdom).groupby(lab_ref).sum() / wdom.sum() * 100).round(1).to_string())


# ---------------- perfis dos grupos com IC 95% (percentis do bootstrap, grupos alinhados)
P = np.array(perfs); w0 = df_matriz_A['domicilios_expandidos'].to_numpy(float)
P0 = perfil_grupos(df_matriz_A, lab_ref, w0, K)
rows = []
for g in range(K):
    for j, pn in enumerate(PERF):
        rows.append({'grupo': g, 'variavel': pn, 'estimativa': P0[g, j], 'ic_inf': np.nanpercentile(P[:, g, j], 2.5), 'ic_sup': np.nanpercentile(P[:, g, j], 97.5)})
pd.DataFrame(rows).to_csv('perfis_grupos_ic.csv', index=False)
print('\nPERFIS DOS GRUPOS COM IC 95% (bootstrap de domicílios; média ponderada por domicílios expandidos)')
print(pd.DataFrame(rows).pivot(index='variavel', columns='grupo', values=['estimativa', 'ic_inf', 'ic_sup']).round(2).to_string())

# ---------------- (c) estabilidade (subamostras de 80%) contra o nulo de covariância
def ari_subamostra(Z, k, n_rep, seed):
    r = np.random.default_rng(seed); full = KMeans(k, n_init=20, random_state=SEED).fit(Z).labels_; v = []
    for _ in range(n_rep):
        idx = r.choice(len(Z), int(0.8 * len(Z)), replace=False)
        v.append(adjusted_rand_score(full[idx], KMeans(k, n_init=5, random_state=int(r.integers(1e6))).fit(Z[idx]).labels_))
    return np.mean(v)
obs = {k: ari_subamostra(Z0, k, 100, 1) for k in (3, 4, 5)}
cov = np.cov(Z0, rowvar=False); rn = np.random.default_rng(7); nulo = {k: [] for k in (3, 4, 5)}
for r_ in range(100):
    Zn = rn.multivariate_normal(np.zeros(Z0.shape[1]), cov, size=len(Z0))
    for k in nulo: nulo[k].append(ari_subamostra(Zn, k, 20, r_))
print('\n(c) ESTABILIDADE (ARI por subamostras de 80%) vs. NULO DE COVARIANCIA (100 conjuntos nulos)')
lin = []
for k in nulo:
    n_ = np.array(nulo[k]); p = (np.sum(n_ >= obs[k]) + 1) / (len(n_) + 1)
    print(f'  K={k}: observado {obs[k]:.3f} | nulo média {n_.mean():.3f} (p95 {np.percentile(n_,95):.3f}) | p = {p:.3f}'); lin.append({'K': k, 'ARI_obs': obs[k], 'ARI_nulo_medio': n_.mean(), 'ARI_nulo_p95': np.percentile(n_, 95), 'p': p})
pd.DataFrame(lin).to_csv('estabilidade_vs_nulo_cov.csv', index=False)

# ---------------- (d) encolhimento empírico-Bayes das médias dos segmentos (n pequeno)
M1 = df_matriz_A.copy(); vars_enc = [v for v in VAR_DOM + VAR_MOR if v in M1.columns and v not in ('renda_domiciliar_pc',)]
n_h = M1['amostra_n'].to_numpy(float)
def encolher(col, base):
    g = base.groupby(CHAVES)[col].var(ddof=1).reindex(pd.MultiIndex.from_frame(M1[CHAVES])).to_numpy()
    s2 = np.nanmean(g); se2 = np.where(np.isnan(g), s2, g) / n_h
    tau2 = max(np.var(M1[col].to_numpy(float), ddof=1) - np.mean(se2), 1e-9)
    lam = tau2 / (tau2 + se2); mu = np.average(M1[col], weights=n_h)
    return lam * M1[col].to_numpy(float) + (1 - lam) * mu, lam
lams = []
for v in vars_enc:
    base = df_dom if v in df_dom.columns else df_mor
    M1[v], lam = encolher(v, base); lams.append(lam)
M1 = derivar_variaveis_modelo(M1)
Zs_ = pca_matriz(M1); lab_s = KMeans(K, n_init=50, random_state=SEED).fit(Zs_).labels_
ct = pd.crosstab(lab_ref, lab_s).reindex(index=range(K), columns=range(K), fill_value=0).to_numpy(); r_, c_ = linear_sum_assignment(-ct); mp = {c: r for r, c in zip(r_, c_)}
lab_s_al = np.array([mp[l] for l in lab_s]); mud = np.where(lab_s_al != lab_ref)[0]
print(f'\n(d) ENCOLHIMENTO EMPÍRICO-BAYES: fator médio lambda = {np.mean(lams):.2f} (segmentos com n<10: {np.mean(np.array(lams)[:, n_h < 10]):.2f}); ARI vs. referência {adjusted_rand_score(lab_ref, lab_s):.3f}; mudam: {len(mud)}')
for i in mud: print('    ', df_matriz_A.loc[i, RA], '|', df_matriz_A.loc[i, SUB], f'n={int(n_h[i])}: G{lab_ref[i]} -> G{lab_s_al[i]}')

# ---------------- (e) regra principal: n<10 = 'não classificável'; coligação no bootstrap < 0,60 = 'fronteira'
res['nao_classificavel'] = res['amostra_n'] < 10
res['fronteira'] = (~res['nao_classificavel']) & (res['junto_do_proprio_grupo'] < 0.60)
res['grupo_encolhido'] = lab_s_al
res.to_csv('robustez_amostral_segmentos.csv', index=False)
print(f"\n(e) NÃO CLASSIFICÁVEIS (n<10): {int(res.nao_classificavel.sum())} de 49 segmentos; FRONTEIRA (coligação<0,60): {int(res.fronteira.sum())}")
print('Domicílios expandidos nos não classificáveis (%):', round(w0[res.nao_classificavel.to_numpy()].sum() / w0.sum() * 100, 1), '| em fronteira (%):', round(w0[res.fronteira.to_numpy()].sum() / w0.sum() * 100, 1))
print(res[res.nao_classificavel][['segmento', 'grupo_ref', 'amostra_n', 'junto_do_proprio_grupo']].round(2).to_string(index=False))
