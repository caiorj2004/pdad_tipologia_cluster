# -*- coding: utf-8 -*-
"""Testes de estresse sobre a partição de UMA variante do modelo (rodar dentro de resultados/<variante>/).

Entrada : segmentos_completo.csv (programa 03) e variáveis de ambiente K_FINAL.
Saída   : stress_resumo.csv (1 linha, métricas comparáveis entre variantes), stress_*.csv (detalhes), saida_21.txt

Testes (cada um responde a um ponto da avaliação da orientadora):
  S1 sementes        : o K-Means converge sempre ao mesmo ótimo? (instabilidade algorítmica)
  S2 tamanho x ponta : segmentos pequenos são sistematicamente mais extremos? (C1; Spearman e Mann-Whitney)
  S3 Wilson          : ICs de Wilson dos segmentos pequenos (banda larga, água da rede)
  S4 aparado         : K-Means aparado (trimmed), alpha 6% e 10%; quais segmentos ficam sem grupo? (alternativas de agrupamento)
  S5 arquétipos      : análise de arquétipos k=2..5; pureza dos segmentos (alternativa ao pertencimento rígido)
  S6 deixa-um-fora   : deixar uma RA (ou uma subpopulação) fora e reagrupar (a partição depende de poucas unidades?)
  S7 algoritmos      : Ward, ligação média, GMM, K-Means ponderado: concordância (ARI) com o K-Means
  S8 corte de n      : excluir segmentos com n<10 ou n<20 e reagrupar
  S9 puxão da RA     : pares de segmentos da mesma RA caem no mesmo grupo mais do que o acaso? (efeito mecânico da RA, C2/I1)
  S10 eixos          : quanto cada componente separa os grupos (eta^2) e quanto do ordenamento é um eixo único
"""
import os, json, numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, norm
from scipy.optimize import nnls
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score as ARI, silhouette_score

RA, SUB = 'a01ra', 'a01subpopulacao'
K = int(os.environ.get('K_FINAL', 3)); SEED = 42
seg = pd.read_csv('segmentos_completo.csv')
zc = [c for c in seg.columns if '_MC_PC' in c]
Z = seg[zc].to_numpy(float); N = len(seg)
lab = seg[f'Cluster_K{K}'].to_numpy(int)
n = seg['amostra_n'].to_numpy(float); w = seg['domicilios_expandidos'].to_numpy(float)
nome = (seg[RA] + ' | ' + seg[SUB]).to_numpy()
R = {'K': K, 'n_segmentos': N, 'n_componentes': len(zc)}
rng = np.random.default_rng(SEED)

def km(Zm, k, seed=SEED, n_init=50, sw=None):
    return KMeans(k, n_init=n_init, random_state=seed).fit(Zm, sample_weight=sw)

# ---------------- S1 sementes
best = km(Z, K); ine = []; aris = []
for s in range(200):
    m = KMeans(K, n_init=1, init='random', random_state=s).fit(Z); ine.append(m.inertia_); aris.append(ARI(best.labels_, m.labels_))
ine = np.array(ine)
R['S1_pct_sementes_no_otimo'] = float(np.mean(ine <= best.inertia_ * (1 + 1e-6)))
R['S1_ARI_medio_semente_unica'] = float(np.mean(aris)); R['S1_ARI_min_semente_unica'] = float(np.min(aris))
R['S1_ARI_ref_vs_otimo'] = float(ARI(lab, best.labels_))
# com n_init=10 (valor usado nas réplicas do bootstrap): quantas sementes chegam ao ótimo de n_init=50?
i10 = np.array([KMeans(K, n_init=10, random_state=s).fit(Z).inertia_ for s in range(100)])
R['S1_pct_ninit10_no_otimo'] = float(np.mean(i10 <= best.inertia_ * (1 + 1e-6)))

# ---------------- S2 tamanho x ponta (C1)
def dist_proprio_centro(Zm, labels, centers):
    return np.linalg.norm(Zm - centers[labels], axis=1)
linhas = []
for kk in sorted({2, 3, 4, 5, K}):
    m = km(Z, kk); d = dist_proprio_centro(Z, m.labels_, m.cluster_centers_)
    rho, p = spearmanr(n, d); peq = n < 20
    u = mannwhitneyu(d[peq], d[~peq], alternative='greater')
    linhas.append({'K': kk, 'spearman_n_dist': rho, 'p': p, 'dist_mediana_n<20': np.median(d[peq]), 'dist_mediana_n>=20': np.median(d[~peq]), 'p_mw_n<20_maior': u.pvalue, 'n_pequenos': int(peq.sum())})
d_glob = np.linalg.norm(Z - Z.mean(0), axis=1); rho_g, p_g = spearmanr(n, d_glob)
pd.DataFrame(linhas).to_csv('stress_S2_tamanho_vs_ponta.csv', index=False)
kf = [l for l in linhas if l['K'] == K][0]
R.update({'S2_spearman_n_dist_Kfinal': kf['spearman_n_dist'], 'S2_p_Kfinal': kf['p'], 'S2_spearman_n_dist_global': rho_g, 'S2_p_global': p_g,
          'S2_dist_med_peq': kf['dist_mediana_n<20'], 'S2_dist_med_grande': kf['dist_mediana_n>=20'], 'S2_p_mw': kf['p_mw_n<20_maior']})

# ---------------- S3 Wilson (segmentos pequenos)
def wilson(p, nn, z=1.96):
    c = (p + z * z / (2 * nn)) / (1 + z * z / nn); h = z * np.sqrt(p * (1 - p) / nn + z * z / (4 * nn * nn)) / (1 + z * z / nn)
    return max(0, c - h), min(1, c + h)
rows = []
for v in ['internet_banda_larga_fixa', 'agua_rede', 'agua_nivel_rede']:
    if v not in seg.columns: continue
    for i in np.where(n < 12)[0]:
        lo, hi = wilson(float(np.clip(seg.loc[i, v], 0, 1)), n[i]); rows.append({'segmento': nome[i], 'n': int(n[i]), 'variavel': v, 'p': seg.loc[i, v], 'ic_inf': lo, 'ic_sup': hi, 'largura': hi - lo, 'grupo': lab[i]})
sw_ = pd.DataFrame(rows); sw_.to_csv('stress_S3_wilson_pequenos.csv', index=False)
R['S3_largura_media_IC_n<12'] = float(sw_['largura'].mean()) if len(sw_) else np.nan

# ---------------- S4 K-Means aparado
def trimmed_kmeans(Zm, k, alpha, n_start=60, seed=SEED):
    r = np.random.default_rng(seed); Nn = len(Zm); h = int(np.ceil((1 - alpha) * Nn)); best = (np.inf, None, None)
    for _ in range(n_start):
        C = Zm[r.choice(Nn, k, replace=False)]
        for _it in range(100):
            D = ((Zm[:, None, :] - C[None]) ** 2).sum(2); a = D.argmin(1); dm = D.min(1)
            keep = np.argsort(dm)[:h]; Cn = C.copy()
            for j in range(k):
                mj = keep[a[keep] == j]
                if len(mj): Cn[j] = Zm[mj].mean(0)
            if np.allclose(Cn, C): break
            C = Cn
        D = ((Zm[:, None, :] - C[None]) ** 2).sum(2); a = D.argmin(1); dm = D.min(1); keep = np.argsort(dm)[:h]; obj = dm[keep].sum()
        if obj < best[0]: best = (obj, a.copy(), keep.copy())
    a = best[1]; mask = np.zeros(Nn, bool); mask[best[2]] = True
    return a, mask
rows = []
for alpha in (0.06, 0.10):
    a, keep = trimmed_kmeans(Z, K, alpha); aparados = np.where(~keep)[0]
    rows.append({'alpha': alpha, 'n_aparados': len(aparados), 'ARI_vs_kmeans_nos_mantidos': ARI(lab[keep], a[keep]),
                 'aparados_n<10': int((n[aparados] < 10).sum()), 'aparados_n<20': int((n[aparados] < 20).sum()),
                 'aparados': '; '.join(f"{nome[i]} (n={int(n[i])})" for i in aparados)})
st4 = pd.DataFrame(rows); st4.to_csv('stress_S4_kmeans_aparado.csv', index=False)
for _, r_ in st4.iterrows():
    tag = f"{int(r_['alpha'] * 100)}"; R[f'S4_aparado{tag}_ARI'] = r_['ARI_vs_kmeans_nos_mantidos']; R[f'S4_aparado{tag}_pct_n<20'] = r_['aparados_n<20'] / max(1, r_['n_aparados'])

# ---------------- S5 arquétipos (Cutler & Breiman), via NNLS alternado
def simplex_nnls(A, b, M=1e3):
    """min ||A x - b||, x >= 0, sum x = 1 (restrição imposta por linha com peso M)."""
    A2 = np.vstack([A, M * np.ones((1, A.shape[1]))]); b2 = np.append(b, M); x, _ = nnls(A2, b2); return x
def arquetipos(X, k, n_start=20, iters=60, seed=SEED):
    r = np.random.default_rng(seed); Nn = len(X); tss = ((X - X.mean(0)) ** 2).sum(); best = (np.inf, None, None)
    for _ in range(n_start):
        idx = list(r.choice(Nn, k, replace=False))
        P = X[idx].copy()
        for _it in range(iters):
            A = np.array([simplex_nnls(P.T, X[i]) for i in range(Nn)])         # pesos de cada segmento
            Pn = np.linalg.lstsq(A, X, rcond=None)[0]                           # arquétipos livres
            B = np.array([simplex_nnls(X.T, Pn[j]) for j in range(k)])          # arquétipos como mistura convexa de segmentos
            Pn2 = B @ X
            if np.allclose(Pn2, P, atol=1e-6): P = Pn2; break
            P = Pn2
        A = np.array([simplex_nnls(P.T, X[i]) for i in range(Nn)]); rss = ((X - A @ P) ** 2).sum()
        if rss < best[0]: best = (rss, A, P)
    return best[0] / tss, best[1], best[2]
rows = []; Aw = None
for k in range(2, 6):
    rs, A, P = arquetipos(Z, k)
    rows.append({'k': k, 'RSS_relativo': rs, 'var_explicada': 1 - rs, 'pureza_media(max peso)': float(A.max(1).mean()), 'pct_segmentos_peso>0.8': float((A.max(1) > 0.8).mean()),
                 'ARI_argmax_vs_kmeans': ARI(A.argmax(1), lab) if k == K else np.nan})
    if k == K: Aw = A
st5 = pd.DataFrame(rows); st5.to_csv('stress_S5_arquetipos.csv', index=False)
if Aw is not None:
    pd.DataFrame(Aw, columns=[f'peso_arq{j}' for j in range(K)]).assign(segmento=nome, n=n, grupo_kmeans=lab).to_csv('stress_S5_arquetipos_pesos.csv', index=False)
    R['S5_pureza_media_Kfinal'] = float(Aw.max(1).mean()); R['S5_pct_peso>0.8_Kfinal'] = float((Aw.max(1) > 0.8).mean()); R['S5_ARI_argmax_vs_kmeans'] = float(ARI(Aw.argmax(1), lab))
    R['S5_var_explicada_Kfinal'] = float(st5.loc[st5.k == K, 'var_explicada'].iloc[0])

# ---------------- S6 deixa-um-fora
def deixa_fora(col):
    out = []
    for g in sorted(seg[col].unique()):
        m = (seg[col] != g).to_numpy()
        if m.sum() < 3 * K: continue
        l2 = km(Z[m], K, n_init=20).labels_; out.append({col: g, 'n_removidos': int((~m).sum()), 'ARI_vs_ref_restantes': ARI(lab[m], l2)})
    return pd.DataFrame(out)
lo_ra = deixa_fora(RA); lo_sub = deixa_fora(SUB)
lo_ra.to_csv('stress_S6_deixa_uma_RA.csv', index=False); lo_sub.to_csv('stress_S6_deixa_uma_subpop.csv', index=False)
R.update({'S6_ARI_medio_sem_uma_RA': lo_ra.ARI_vs_ref_restantes.mean(), 'S6_ARI_min_sem_uma_RA': lo_ra.ARI_vs_ref_restantes.min(),
          'S6_RA_critica': lo_ra.loc[lo_ra.ARI_vs_ref_restantes.idxmin(), RA], 'S6_ARI_medio_sem_uma_subpop': lo_sub.ARI_vs_ref_restantes.mean(), 'S6_ARI_min_sem_uma_subpop': lo_sub.ARI_vs_ref_restantes.min()})

# ---------------- S7 algoritmos
alg = {'Ward': AgglomerativeClustering(K, linkage='ward').fit_predict(Z), 'ligacao_media': AgglomerativeClustering(K, linkage='average').fit_predict(Z),
       'GMM_completo': GaussianMixture(K, covariance_type='full', n_init=10, random_state=SEED, reg_covar=1e-3).fit_predict(Z),
       'GMM_diagonal': GaussianMixture(K, covariance_type='diag', n_init=10, random_state=SEED).fit_predict(Z),
       'KMeans_peso_raiz': km(Z, K, sw=np.sqrt(w)).labels_, 'KMeans_peso_integral': km(Z, K, sw=w).labels_, 'KMeans_peso_n_amostral': km(Z, K, sw=n).labels_}
st7 = pd.DataFrame([{'algoritmo': a, 'ARI_vs_kmeans': ARI(lab, l), 'tamanhos': str(sorted(np.bincount(l, minlength=K).tolist(), reverse=True))} for a, l in alg.items()])
st7.to_csv('stress_S7_algoritmos.csv', index=False)
R['S7_ARI_medio_outros_algoritmos'] = float(st7.ARI_vs_kmeans.mean()); R['S7_ARI_min_outros_algoritmos'] = float(st7.ARI_vs_kmeans.min())
for a in ('Ward', 'GMM_completo', 'KMeans_peso_integral'): R[f'S7_ARI_{a}'] = float(st7.loc[st7.algoritmo == a, 'ARI_vs_kmeans'].iloc[0])

# ---------------- S8 corte de n
rows = []
for corte in (10, 20):
    m = n >= corte; mm = km(Z[m], K, n_init=50); rows.append({'corte_n': corte, 'segmentos_mantidos': int(m.sum()), 'ARI_vs_ref_nos_mantidos': ARI(lab[m], mm.labels_), 'silhueta': silhouette_score(Z[m], mm.labels_), 'tamanhos': str(sorted(np.bincount(mm.labels_, minlength=K).tolist(), reverse=True))})
st8 = pd.DataFrame(rows); st8.to_csv('stress_S8_corte_n.csv', index=False)
for _, r_ in st8.iterrows(): R[f"S8_ARI_n>={int(r_['corte_n'])}"] = r_['ARI_vs_ref_nos_mantidos']

# ---------------- S9 puxão da RA (efeito mecânico)
def frac_mesmo(col, labels):
    g = seg[col].to_numpy(); tot = same = 0
    for a in range(N):
        for b in range(a + 1, N):
            if g[a] == g[b]: tot += 1; same += labels[a] == labels[b]
    return same / tot if tot else np.nan, tot
obs_ra, npar_ra = frac_mesmo(RA, lab); obs_sub, npar_sub = frac_mesmo(SUB, lab)
perm_ra = []; perm_sub = []
for _ in range(1000):
    pl = rng.permutation(lab); perm_ra.append(frac_mesmo(RA, pl)[0]); perm_sub.append(frac_mesmo(SUB, pl)[0])
R.update({'S9_pares_mesma_RA': npar_ra, 'S9_frac_mesmo_grupo_RA': obs_ra, 'S9_frac_esperada_RA': float(np.mean(perm_ra)), 'S9_p_RA': float((np.sum(np.array(perm_ra) >= obs_ra) + 1) / 1001),
          'S9_frac_mesmo_grupo_subpop': obs_sub, 'S9_frac_esperada_subpop': float(np.mean(perm_sub)), 'S9_p_subpop': float((np.sum(np.array(perm_sub) >= obs_sub) + 1) / 1001)})

# ---------------- S10 eixos
rows = []
for j, c in enumerate(zc):
    x = Z[:, j]; ssb = sum(((x[lab == g].mean() - x.mean()) ** 2) * (lab == g).sum() for g in np.unique(lab)); rows.append({'componente': c, 'eta2_grupos': ssb / ((x - x.mean()) ** 2).sum()})
st10 = pd.DataFrame(rows); st10.to_csv('stress_S10_eixos.csv', index=False)
R['S10_eta2_max'] = float(st10.eta2_grupos.max()); R['S10_eixo_dominante'] = st10.loc[st10.eta2_grupos.idxmax(), 'componente']
ev = np.linalg.eigvalsh(np.cov(Z.T))[::-1]; R['S10_pct_var_1o_eixo_da_matriz'] = float(ev[0] / ev.sum())

pd.DataFrame([R]).to_csv('stress_resumo.csv', index=False)
with open('saida_21.txt', 'w', encoding='utf-8') as f:
    for k_, v_ in R.items(): f.write(f'{k_}: {v_}\n')
print(json.dumps({k_: (round(v_, 3) if isinstance(v_, float) else v_) for k_, v_ in R.items()}, ensure_ascii=False, indent=1))
