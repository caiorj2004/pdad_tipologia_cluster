# -*- coding: utf-8 -*-
# Sensibilidade da partição aos ajustes da matriz (Etapa A da revisão): liga/desliga cada ajuste e compara com a partição anterior.
import os, sys, io, contextlib, numpy as np, pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
src = open(os.environ.get('TCC_03', '03_analise_cluster_4d.py'), encoding='utf-8').read()
corte = src[:src.index('# ETAPA 5 — BENCHMARK DE K')]
ref = pd.read_csv(os.environ.get('REF_CSV', 'segmentos_completo.csv'))   # partição anterior (K=5)
CFGS = [('base (versão anterior)', ''), ('só ctx', 'ctx'), ('só agua', 'agua'), ('só esgoto', 'esgoto'), ('só renda', 'renda'), ('só dur', 'dur'),
        ('agua+esgoto', 'agua,esgoto'), ('todos', 'ctx,agua,esgoto,renda,dur')]
linhas = []
for nome, flags in CFGS:
    os.environ['V2'] = flags
    g = {'__name__': 'cfg'}
    with contextlib.redirect_stdout(io.StringIO()) as buf: exec(compile(corte, '03', 'exec'), g)
    M = g['df_A_mc']; comp = [c for c in M.columns if '_MC_PC' in c]
    Z = M[comp].to_numpy(); m = M[g['CHAVES']].merge(ref[g['CHAVES'] + ['Cluster_K5']], on=g['CHAVES'], how='left')
    row = {'config': nome, 'n_comp': len(comp), 'comps_por_dim': ','.join(str(sum(c.startswith(d + '_MC') for c in comp)) for d in g['DIMENSOES_REF'])}
    for K in (3, 4, 5):
        lab = KMeans(K, n_init=50, random_state=42).fit(Z).labels_
        row[f'sil_K{K}'] = round(silhouette_score(Z, lab), 3)
        if K == 5: row['ARI_vs_anterior'] = round(adjusted_rand_score(m['Cluster_K5'], lab), 3); row['tamanhos_K5'] = sorted(np.bincount(lab).tolist(), reverse=True); lab5 = lab
        if K == 4: row['tam_K4'] = sorted(np.bincount(lab).tolist(), reverse=True)
    linhas.append(row); print(row, flush=True)
    pd.DataFrame({'segmento': M[g['RA']] + ' | ' + M[g['SUB']], 'n': M['amostra_n'], 'K5': lab5}).to_csv(f"ablacao_labels_{flags.replace(',', '_') or 'base'}.csv", index=False)
pd.DataFrame(linhas).to_csv('ablacao_matriz.csv', index=False)
