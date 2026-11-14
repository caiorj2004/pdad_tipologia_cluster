# -*- coding: utf-8 -*-
"""Validação externa sob incerteza da partição (rodar em resultados/<variante>/ depois de 16 e 19).
A validação do programa 19 testa os desfechos externos contra a partição de REFERÊNCIA. Aqui o desfecho de cada segmento (validacao_externa_segmentos.csv)
fica fixo e os grupos mudam: usa-se o rótulo de cada réplica do bootstrap de domicílios (boot_labels_alinhados.npy, já realinhados à referência) e refaz-se o
Kruskal-Wallis entre grupos. Resposta: o resultado de 19 persiste quando a partição é refeita com a incerteza amostral, ou depende da partição única?
Saídas: validacao_externa_incerteza_28.csv (por desfecho) e validacao_externa_incerteza_28_resumo.csv (chaves ext28_*).
Limites: o desfecho de cada segmento não é re-amostrado; só a partição varia. Sem correção para comparações múltiplas. 49 segmentos."""
import numpy as np, pandas as pd
from scipy.stats import kruskal
VARS = ['ebia_insegur_qualquer', 'ebia_mod_grave', 'beneficio_social', 'cadunico', 'escritura_registrada', 'doc_formal', 'outorga_agua']
cfg = dict(kv.split('=', 1) for kv in open('config.txt').read().split()); K = int(cfg['K'])
seg = pd.read_csv('segmentos_completo.csv'); ve = pd.read_csv('validacao_externa_segmentos.csv')
CH = ['a01ra', 'a01subpopulacao']
ve = seg[CH].merge(ve[CH + VARS], on=CH, how='left')          # desfechos na ordem de segmentos_completo
L = np.load('boot_labels_alinhados.npy'); ref = seg[f'Cluster_K{K}'].values
assert L.shape[1] == len(ref) and (L == ref).mean() > 0.5, 'ordem dos segmentos incompatível com a referência'
def kw(lab, y):
    ok = ~np.isnan(y); g = [y[ok & (lab == k)] for k in range(K) if (ok & (lab == k)).sum() > 1]
    if len(g) < 2: return np.nan, np.nan
    try: H, p = kruskal(*g); return H, p
    except ValueError: return np.nan, np.nan
def eps2(H, n, k): return (H - k + 1) / (n - k) if np.isfinite(H) else np.nan
rows = []
for v in VARS:
    y = ve[v].to_numpy(float); n = int((~np.isnan(y)).sum())
    H0, p0 = kw(ref, y)
    res = np.array([kw(L[b], y) for b in range(L.shape[0])]); p = res[:, 1]; ok = np.isfinite(p)
    rows.append(dict(desfecho=v, p_referencia=p0, epsilon2_referencia=eps2(H0, n, K), p_mediana_replicas=float(np.nanmedian(p)), p_q25=float(np.nanpercentile(p, 25)), p_q75=float(np.nanpercentile(p, 75)),
                     pct_replicas_p_menor_005=float(100 * np.mean(p[ok] < 0.05)), epsilon2_mediana=float(np.nanmedian([eps2(h, n, K) for h in res[:, 0]])), n_replicas=int(ok.sum())))
out = pd.DataFrame(rows); out.to_csv('validacao_externa_incerteza_28.csv', index=False)
R = dict(ext28_n_desfechos_p_ref_menor_005=int((out.p_referencia < 0.05).sum()), ext28_n_desfechos_mediana_p_menor_005=int((out.p_mediana_replicas < 0.05).sum()),
         ext28_pct_replicas_p005_media=float(out.pct_replicas_p_menor_005.mean()), ext28_epsilon2_mediano_medio=float(out.epsilon2_mediana.mean()))
pd.DataFrame([R]).to_csv('validacao_externa_incerteza_28_resumo.csv', index=False)
print('VALIDAÇÃO EXTERNA SOB INCERTEZA DA PARTIÇÃO (Kruskal-Wallis por segmento; %d réplicas)' % L.shape[0]); print(out.round(3).to_string(index=False)); print(R)
