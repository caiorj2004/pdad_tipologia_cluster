# -*- coding: utf-8 -*-
"""Testa se a vantagem de validação externa/territorial de um modelo vem do CONTEXTO DA RA (constante dentro da RA) e não do perfil dos domicílios.
(1) Validação externa condicional à RA: Kruskal-Wallis dos desfechos (média por segmento, validacao_externa_segmentos.csv) entre grupos, com p de PERMUTAÇÃO DENTRO DA RA
    (rótulos embaralhados só entre segmentos da mesma RA; preserva a composição por RA). RAs com um só segmento não contribuem. Pergunta: a associação sobrevive depois de
    retirar o que é explicado por 'estar na RA X'?
(2) eta² da RA nos escores de cada dimensão (ANOVA de um fator, 49 segmentos, RA como fator), por variante, bruto e ajustado pelos graus de liberdade (o bruto já vale ≈ 0,46 sob
    independência com 23 RAs e 49 segmentos). Mede quanto da matriz é 'RA' por construção.
Env: N_PERM29 (padrão 5000). Saídas: validacao_condicional_RA_29.csv, eta2_RA_escores_29.csv, validacao_condicional_RA_29_resumo.csv (chaves cond29_*, eta29_*)."""
import os, numpy as np, pandas as pd
from scipy.stats import kruskal
VARS = ['ebia_insegur_qualquer', 'ebia_mod_grave', 'beneficio_social', 'cadunico', 'escritura_registrada', 'doc_formal', 'outorga_agua']
NP = int(os.environ.get('N_PERM29', 5000)); rng = np.random.default_rng(42)
cfg = dict(kv.split('=', 1) for kv in open('config.txt').read().split()); K = int(cfg['K'])
seg = pd.read_csv('segmentos_completo.csv'); ve = pd.read_csv('validacao_externa_segmentos.csv'); CH = ['a01ra', 'a01subpopulacao']
d = seg[CH + [f'Cluster_K{K}']].merge(ve[CH + VARS], on=CH, how='left'); lab = d[f'Cluster_K{K}'].to_numpy(int); ra = d['a01ra'].to_numpy()
blocos = [np.where(ra == r)[0] for r in np.unique(ra) if (ra == r).sum() > 1]
def H(l, y):
    ok = ~np.isnan(y); g = [y[ok & (l == k)] for k in range(K) if (ok & (l == k)).sum() > 1]
    return kruskal(*g)[0] if len(g) > 1 else np.nan
rows = []
for v in VARS:
    y = d[v].to_numpy(float); h0 = H(lab, y); ge = 0; tot = 0
    for _ in range(NP):
        l = lab.copy()
        for b in blocos: l[b] = rng.permutation(lab[b])
        h = H(l, y)
        if np.isfinite(h): tot += 1; ge += h >= h0
    rows.append(dict(desfecho=v, H_obs=h0, p_perm_dentro_RA=(1 + ge) / (1 + tot), n_segmentos_em_RA_multi=int(sum(len(b) for b in blocos))))
out = pd.DataFrame(rows); out.to_csv('validacao_condicional_RA_29.csv', index=False)
pcs = [c for c in seg.columns if '_MC_PC' in c]; er = []
for c in pcs:
    x = seg[c].to_numpy(float); tot = ((x - x.mean()) ** 2).sum(); ent = sum(((x[ra == r] - x[ra == r].mean()) ** 2).sum() for r in np.unique(ra))
    e2 = 1 - ent / tot; g_ = len(np.unique(ra)); n_ = len(x); er.append(dict(componente=c, eta2_RA=e2, eta2_RA_ajustado=1 - (1 - e2) * (n_ - 1) / (n_ - g_)))
er = pd.DataFrame(er); er.to_csv('eta2_RA_escores_29.csv', index=False)
R = dict(cond29_n_desfechos_p_menor_005=int((out.p_perm_dentro_RA < 0.05).sum()), cond29_p_min=float(out.p_perm_dentro_RA.min()), eta29_RA_medio=float(er.eta2_RA.mean()), eta29_RA_max=float(er.eta2_RA.max()), eta29_RA_ajust_medio=float(er.eta2_RA_ajustado.mean()))
for _, r_ in out.iterrows(): R[f'cond29_p_{r_.desfecho}'] = float(r_.p_perm_dentro_RA)
pd.DataFrame([R]).to_csv('validacao_condicional_RA_29_resumo.csv', index=False)
print(out.round(4).to_string(index=False)); print(er.round(3).to_string(index=False)); print(R)
