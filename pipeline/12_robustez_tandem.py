# -*- coding: utf-8 -*-
"""Robustez da sequência PCA -> K-Means (tandem analysis) no TCC PDAD Rural 2022.

Compara a partição de referência (K-Means, K=3, sobre os 6 componentes principais por dimensão)
com partições obtidas SEM a etapa de PCA ou com redução + agrupamento no mesmo passo.
Entrada: segmentos_completo.csv (49 segmentos x variáveis, gerado pelo script 4D).
Saída: tabelas em CSV e impressão no terminal.
"""
import numpy as np, pandas as pd, warnings
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score as ARI, silhouette_score
warnings.filterwarnings("ignore")

SEED, K = 42, 3
CSV = "segmentos_completo.csv"
d = pd.read_csv(CSV)

D1 = ["log_renda_media", "possui_automovel", "possui_motocicleta", "possui_freezer", "possui_lavadora", "possui_tv", "possui_ar_condicionado", "possui_microondas"]
D2 = ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia"]
D3 = ["comodos_por_pessoa", "banheiros_por_pessoa", "agua_ilr_adequacao", "agua_ilr_tipo", "captacao_chuva", "caixa_agua",
      "esgoto_ilr_adequacao", "esgoto_ilr_tipo", "acesso_energetico_estruturado", "destino_com_coleta_direta", "rua_pavimentada",
      "iluminacao_publica", "drenagem", "ponto_onibus"]
D4 = ["internet_domiciliar", "internet_banda_larga_fixa", "internet_movel", "possui_computador"]
DIMS = {"D1": D1, "D2": D2, "D3": D3, "D4": D4}
VARS = D1 + D2 + D3 + D4
assert len(VARS) == 30
X = d[VARS].apply(pd.to_numeric, errors="coerce").fillna(0)
assert (X.nunique() > 1).all()
Z = StandardScaler().fit_transform(X)                    # 49 x 34, z-score (igual ao script principal)
cols = {v: i for i, v in enumerate(VARS)}

pcs = d.filter(like="_MC_PC")
Xpca = pcs.values
ref = KMeans(K, n_init=50, random_state=SEED).fit(Xpca).labels_
assert ARI(ref, d["Cluster_K3"].values) == 1.0, "partição de referência não reproduzida"
print("Referência reproduzida (ARI=1). Tamanhos:", np.bincount(ref).tolist())

rows = []
def reg(nome, lab, esp=None, nota=""):
    sil = silhouette_score(esp, lab) if esp is not None and len(set(lab)) > 1 else np.nan
    rows.append(dict(variante=nome, ARI_vs_ref=round(ARI(ref, lab), 3), n_discordantes=int(None) if False else None,
                     tamanhos=np.sort(np.bincount(lab, minlength=K))[::-1].tolist(), silhueta_no_proprio_espaco=round(sil, 3), nota=nota))
    return lab

def km(M, n_init=200, seed=SEED): return KMeans(K, n_init=n_init, random_state=seed).fit(M).labels_

def discord(lab):
    """segmentos fora do pareamento ótimo de rótulos (húngaro) com a referência"""
    from scipy.optimize import linear_sum_assignment
    C = np.zeros((K, K), int)
    for a, b in zip(ref, lab): C[a, b] += 1
    r, c = linear_sum_assignment(-C)
    m = dict(zip(c, r))
    return [i for i, b in enumerate(lab) if m[b] != ref[i]]

out_disc = {}

# 1) Referência
reg("0. Referência: PCA por dimensão (6 PC) + K-Means", ref, Xpca)

# 2) K-Means direto nas 30 variáveis padronizadas (peso igual por variável)
l = km(Z); reg("1. K-Means nas 30 variáveis padronizadas", l, Z, "peso igual por variável; D3 (14 var.) pesa mais"); out_disc["1"] = discord(l)

# 3) 30 variáveis com peso igual por dimensão (cada bloco com variância total 1)
Zb = Z.copy()
for dim, vs in DIMS.items():
    idx = [cols[v] for v in vs]; Zb[:, idx] = Zb[:, idx] / np.sqrt(len(vs))
l = km(Zb); reg("2. 30 variáveis, peso igual por dimensão", l, Zb, "cada dimensão com variância total 1"); out_disc["2"] = discord(l)

# 4) Mesmas 30 variáveis, Ward
l = AgglomerativeClustering(K, linkage="ward").fit_predict(Z); reg("3. Ward nas 30 variáveis padronizadas", l, Z); out_disc["3"] = discord(l)

# 5) PCA global (sobre as 30 variáveis), nº de componentes pela análise paralela (percentil 95, 1.000 simulações)
rng = np.random.default_rng(SEED)
sim = np.array([np.sort(np.linalg.eigvalsh(np.corrcoef(rng.standard_normal(Z.shape), rowvar=False)))[::-1] for _ in range(1000)])
p95 = np.percentile(sim, 95, axis=0)
eig = np.sort(np.linalg.eigvalsh(np.corrcoef(Z, rowvar=False)))[::-1]
kg = int(np.sum(eig > p95)); print(f"PCA global: análise paralela retém {kg} componentes; variância acumulada = {eig[:kg].sum()/30*100:.1f}%")
Pg = PCA(kg).fit_transform(Z)
l = km(Pg); reg(f"4. PCA global ({kg} PC) + K-Means", l, Pg, "análise paralela sobre as 30 variáveis"); out_disc["4"] = discord(l)

# 6) Sensibilidade ao nº de componentes globais (3 a 12)
for q in (3, 6, 8, 10, 12):
    P = PCA(q).fit_transform(Z); l = km(P)
    reg(f"4b. PCA global ({q} PC) + K-Means", l, P, f"{PCA(q).fit(Z).explained_variance_ratio_.sum()*100:.0f}% da variância")

# 7) Redução + agrupamento no mesmo passo: Reduced K-means (De Soete & Carroll, 1994) e Factorial K-means (Vichi & Kiers, 2001)
def red_fact(Zm, q, kind, n_starts=100, iters=200, seed=SEED):
    """kind='rkm': min ||Z - U Ybar A'||^2 (A = autovetores de maior variância ENTRE grupos)
       kind='fkm': min ||Z A - U Ybar||^2 (A = autovetores de menor variância DENTRO dos grupos)"""
    best = (np.inf, None)
    rs = np.random.default_rng(seed)
    for s in range(n_starts):
        A = np.linalg.qr(rs.standard_normal((Zm.shape[1], q)))[0] if s else PCA(q).fit(Zm).components_.T
        lab = None
        for it in range(iters):
            Y = Zm @ A
            new = KMeans(K, n_init=1, init="k-means++", random_state=int(rs.integers(1e9))).fit(Y).labels_ if lab is None else \
                  KMeans(K, n_init=1, init=np.array([Y[lab == k].mean(0) if (lab == k).any() else Y[rs.integers(len(Y))] for k in range(K)]), random_state=0).fit(Y).labels_
            U = np.eye(K)[new]; n = U.sum(0); n[n == 0] = 1
            M = (U.T @ Zm) / n[:, None]                       # centroides no espaço original (K x p)
            if kind == "rkm":
                B = (M.T * n) @ M                              # dispersão entre grupos (p x p)
                w, V = np.linalg.eigh(B); A_new = V[:, ::-1][:, :q]
            else:
                W = (Zm - U @ M).T @ (Zm - U @ M)              # dispersão dentro dos grupos
                w, V = np.linalg.eigh(W); A_new = V[:, :q]
            if lab is not None and (new == lab).all() and np.allclose(A_new @ A_new.T, A @ A.T, atol=1e-8): break
            A, lab = A_new, new
        Y = Zm @ A; U = np.eye(K)[lab]; n = U.sum(0); n[n == 0] = 1
        if kind == "rkm":
            M = (U.T @ Zm) / n[:, None]; loss = ((Zm - U @ M @ A @ A.T) ** 2).sum()
        else:
            Yb = (U.T @ Y) / n[:, None]; loss = ((Y - U @ Yb) ** 2).sum()
        if loss < best[0] and len(set(lab)) == K: best = (loss, lab.copy(), A.copy())
    return best

for kind, nome in (("rkm", "Reduced K-means"), ("fkm", "Factorial K-means")):
    for q in (6, 8):
        loss, lab, A = red_fact(Z, q, kind)
        reg(f"5. {nome} (q={q})", lab, Z @ A, "redução + agrupamento no mesmo passo");
        if q == 6: out_disc[f"5{kind}"] = discord(lab)

res = pd.DataFrame(rows)
res["n_discordantes"] = res["variante"].map(lambda v: len(out_disc[v[0]]) if v[0] in out_disc and v[:2] in ("1.", "2.", "3.", "4.") else None)
pd.set_option("display.width", 250, "display.max_colwidth", 70)
print(res.to_string(index=False))
res.to_csv("robustez_tandem.csv", index=False)

# quem muda: segmentos discordantes e sua classe (núcleo/fronteira)
cl = pd.read_csv("km_ward_labels.csv") if False else None
seg = (d["a01ra"] + " — " + d["a01subpopulacao"]).values
print("\nSegmentos discordantes da referência (pareamento ótimo de rótulos):")
for k, v in out_disc.items():
    print(f" variante {k}: {len(v)} -> ", "; ".join(f"{seg[i]} (grupo ref {ref[i]})" for i in v))

# variância preservada pelos 6 PC da referência
tot = 34
ret = 0
for dim, vs in DIMS.items():
    sub = StandardScaler().fit_transform(X[vs]); ev = np.sort(np.linalg.eigvalsh(np.corrcoef(sub, rowvar=False)))[::-1]
    kd = sum(c.startswith(f"{dim}_") for c in pcs.columns) if False else sum(c.startswith({"D1": "D1_", "D2": "D2_", "D3": "D3_", "D4": "D4_"}[dim]) for c in pcs.columns)
    ret += ev[:kd].sum(); print(f"  {dim}: {kd} PC retidos, {ev[:kd].sum():.2f} de {len(vs)} unidades de variância ({ev[:kd].sum()/len(vs)*100:.0f}%)")
print(f"Variância total preservada pelos 6 PC (dentro das dimensões): {ret/tot*100:.1f}%")
