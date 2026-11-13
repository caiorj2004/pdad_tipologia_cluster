# -*- coding: utf-8 -*-
"""Decomposição multinível, Moran ajustado por composição e validação externa ordenada (rodar em resultados/<variante>/ depois de 03, 18 e 19).
 (1) MULTINÍVEL: para água em rede, banda larga fixa, escolaridade do adulto e ln da renda (médias dos 49 segmentos): MixedLM com subpopulação como efeito fixo e RA como intercepto aleatório.
     Reporta a fração da variância total devida à subpopulação (R2 marginal), à RA (variância do intercepto, dado o efeito fixo) e ao resíduo, e p de permutação do componente de RA
     (rótulos de RA permutados DENTRO de cada subpopulação). Sem ponderação por tamanho; segmentos pequenos têm ruído amostral maior (ver docs).
 (2) MORAN AJUSTADO POR COMPOSIÇÃO: Moran's I (vizinhança Queen entre RAs, 999 permutações) do valor médio da RA, antes e depois de remover a média da subpopulação de cada segmento. Se o I some
     depois do ajuste, a "autocorrelação" vinha da mistura de subpopulações de cada RA.
 (3) VALIDAÇÃO EXTERNA ORDENADA: os grupos são ordenados pela média de ln(renda) dos segmentos; Spearman entre a posição do grupo e cada indicador externo (49 segmentos), p por permutação dos rótulos dos grupos.
Saída: multinivel_25.csv, moran_ajustado_25.csv, validacao_ordenada_25.csv e resumo em multinivel_moran_resumo.csv"""
import json, warnings, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
warnings.filterwarnings('ignore')
cfg = dict(kv.split('=', 1) for kv in Path('config.txt').read_text().split()); K = int(cfg['K']); CL = f'Cluster_K{K}'
seg = pd.read_csv('segmentos_completo.csv'); rng = np.random.default_rng(42)
seg['ra'] = seg['ra_norm']; seg['sub'] = seg['a01subpopulacao']
SAIDAS = {'agua_rede': 'água em rede', 'internet_banda_larga_fixa': 'banda larga fixa', 'escolaridade_adulto': 'escolaridade do adulto', 'log_renda_domiciliar_pc': 'ln(renda pc)'}
SAIDAS = {k: v for k, v in SAIDAS.items() if k in seg.columns}; N_PERM = 100; R = {}
def ajusta(d, y):
    for met in ('powell', 'bfgs', 'nm'):          # lbfgs falha com matriz singular quando a variância de RA tende a zero
        try: m = smf.mixedlm(f'{y} ~ C(sub)', d, groups=d['ra']).fit(reml=True, method=met); break
        except Exception: m = None
    if m is None: raise RuntimeError('MixedLM não convergiu')
    vre = float(m.cov_re.iloc[0, 0]); res = float(m.scale); fe = np.var(m.model.exog @ m.fe_params.values)
    return fe, vre, res
linhas = []
for y, nome in SAIDAS.items():
    d = seg[['ra', 'sub', y]].dropna().copy(); d[y] = (d[y] - d[y].mean()) / d[y].std()
    fe, vre, res = ajusta(d, y); tot = fe + vre + res; obs = vre / (vre + res)
    nulo = []
    for _ in range(N_PERM):
        dp = d.copy()
        for s_, ix in d.groupby('sub').groups.items(): dp.loc[ix, 'ra'] = rng.permutation(d.loc[ix, 'ra'].to_numpy())
        try: f2, v2, r2 = ajusta(dp, y); nulo.append(v2 / (v2 + r2))
        except Exception: pass
    p = (1 + np.sum(np.array(nulo) >= obs)) / (1 + len(nulo))
    linhas.append({'indicador': nome, 'fracao_subpop(fixo)': fe / tot, 'fracao_RA(aleatorio)': vre / tot, 'fracao_residuo': res / tot, 'CCI_RA_dado_subpop': obs, 'p_perm_RA': p})
ml = pd.DataFrame(linhas); ml.to_csv('multinivel_25.csv', index=False); print('(1) MULTINÍVEL (subpopulação fixa, RA aleatória; indicadores padronizados):'); print(ml.round(3).to_string(index=False))
R['multinivel_CCI_RA_medio'] = float(ml['CCI_RA_dado_subpop'].mean()); R['multinivel_fracao_subpop_media'] = float(ml['fracao_subpop(fixo)'].mean()); R['multinivel_n_p<0.05_RA'] = int((ml.p_perm_RA < 0.05).sum())
# (2) Moran
viz = json.load(open('vizinhos_queen.json')); linhas = []
def moran(v, pares, n_perm=999):
    ras = list(v.index); x = v.to_numpy(float); x = x - x.mean(); ix = {r: i for i, r in enumerate(ras)}
    a = np.array([[ix[p], ix[q]] for p, q in pares]); W = len(a) * 2
    def I(z): return (len(z) / W) * 2 * np.sum(z[a[:, 0]] * z[a[:, 1]]) / np.sum(z ** 2)
    obs = I(x); nul = np.array([I(rng.permutation(x)) for _ in range(n_perm)]); return obs, (1 + np.sum(nul >= obs)) / (1 + n_perm)
for y, nome in SAIDAS.items():
    d = seg[['ra', 'sub', y, 'amostra_n']].dropna().copy(); d['adj'] = d[y] - d.groupby('sub')[y].transform('mean')
    for tipo, col in (('bruto', y), ('ajustado por subpopulação', 'adj')):
        v = d.groupby('ra').apply(lambda g: np.average(g[col], weights=g['amostra_n']), include_groups=False); v = v[[r for r in v.index if r in viz]]
        pares = {tuple(sorted((p, q))) for p in v.index for q in viz[p] if q in v.index}; I, p = moran(v, pares)
        linhas.append({'indicador': nome, 'versao': tipo, 'RAs': len(v), 'vinculos': len(pares), 'I_Moran': I, 'p_perm': p})
mo = pd.DataFrame(linhas); mo.to_csv('moran_ajustado_25.csv', index=False); print('\n(2) MORAN I POR RA (bruto x ajustado pela média da subpopulação):'); print(mo.round(3).to_string(index=False))
b = mo[mo.versao == 'bruto']; a = mo[mo.versao != 'bruto']; R['moran_I_bruto_medio'] = float(b.I_Moran.mean()); R['moran_I_ajustado_medio'] = float(a.I_Moran.mean()); R['moran_ajustado_n_p<0.05'] = int((a.p_perm < 0.05).sum())
# (3) validação ordenada
ve = pd.read_csv('validacao_externa_segmentos.csv'); ve = ve.merge(seg[['a01ra', 'a01subpopulacao', 'log_renda_domiciliar_pc']], on=['a01ra', 'a01subpopulacao'], how='left') if 'log_renda_domiciliar_pc' in seg.columns else ve
ordem = ve.groupby(CL)['log_renda_domiciliar_pc'].mean().sort_values().index.tolist(); pos = ve[CL].map({g: i for i, g in enumerate(ordem)}).to_numpy(float); linhas = []
for v in ['ebia_insegur_qualquer', 'ebia_mod_grave', 'beneficio_social', 'cadunico', 'escritura_registrada', 'doc_formal', 'outorga_agua']:
    m = ve[v].notna().to_numpy(); r0 = spearmanr(pos[m], ve.loc[m, v]).statistic
    nul = np.array([spearmanr(rng.permutation(pos[m]), ve.loc[m, v]).statistic for _ in range(999)]); linhas.append({'indicador': v, 'ordem_grupos(pobre->rico em renda)': str(ordem), 'rho_Spearman': r0, 'p_perm': (1 + np.sum(np.abs(nul) >= abs(r0))) / 1000})
vo = pd.DataFrame(linhas); vo.to_csv('validacao_ordenada_25.csv', index=False); print('\n(3) VALIDAÇÃO EXTERNA ORDENADA (posição do grupo na escala de renda x indicador, 49 segmentos):'); print(vo.drop(columns='ordem_grupos(pobre->rico em renda)').round(3).to_string(index=False))
R['ordenada_n_p<0.05'] = int((vo.p_perm < 0.05).sum()); R['ordenada_rho_medio_abs'] = float(vo.rho_Spearman.abs().mean())
pd.DataFrame([R]).to_csv('multinivel_moran_resumo.csv', index=False)
