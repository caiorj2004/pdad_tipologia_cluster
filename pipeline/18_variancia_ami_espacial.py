# -*- coding: utf-8 -*-
# Decomposição de variância (RA x subpopulação), AMI dos grupos com RA e subpopulação, e join count espacial dos grupos.
# Substitui o V de Cramér do texto anterior. Entrada: segmentos_completo.csv e vizinhos_queen.json (programa 03) e as bases corrigidas (Fase 1 do 03).
# Saída: variancia_ra_sub.csv, ami_rotulos.csv, joincount.csv, saida_18.txt
import os, json, numpy as np, pandas as pd
from sklearn.metrics import adjusted_mutual_info_score as AMI
import statsmodels.formula.api as smf
SEED = 42
src = open(os.environ.get('TCC_03', '03_analise_cluster_4d.py'), encoding='utf-8').read()
exec(compile(src[:src.index('def monte_carlo_parallel_analysis')], '03_fase1', 'exec'))
CLc = f"Cluster_K{K_FINAL}"
seg = pd.read_csv('segmentos_completo.csv'); seg['ra_norm'] = seg[RA].map(limpar_nome_ra_base) if RA in seg.columns else seg['ra_norm']
rng = np.random.default_rng(SEED)

# ---- 1. Decomposição de variância dos escores (nível segmento) e dos indicadores (nível domicílio)
def eta(df, y, w=None):
    """Parcelas de variância de y: R2 de RA, de subpopulação, conjunto (aditivo) e interação (modelo saturado)."""
    f = lambda t: smf.ols(f'{y} ~ {t}', df).fit().rsquared
    r_ra, r_sub, r_ad = f('C(ra)'), f('C(sub)'), f('C(ra) + C(sub)')
    r_full = f('C(ra) * C(sub)') if False else np.nan
    return r_ra, r_sub, r_ad
seg2 = seg.rename(columns={RA: 'ra', SUB: 'sub'})
seg2['ra'] = seg2['ra_norm']
comp = [c for c in seg.columns if '_MC_PC' in c]
rows = []
for y in comp:
    r_ra, r_sub, r_ad = eta(seg2, y)
    m = smf.mixedlm(f'{y} ~ C(sub)', seg2, groups=seg2['ra']).fit(reml=True)
    icc = float(m.cov_re.iloc[0, 0]) / (float(m.cov_re.iloc[0, 0]) + m.scale)
    rows.append({'nivel': 'segmento', 'variavel': y, 'R2_RA': r_ra, 'R2_sub': r_sub, 'R2_RA+sub': r_ad, 'ICC_RA_dado_sub': icc})
# domicílio: variância entre RA, entre subpopulação e entre segmentos (RA x sub)
dd = df_dom.copy(); dd['ra'] = dd['ra_norm'] if 'ra_norm' in dd.columns else dd[RA]; dd['sub'] = dd[SUB]; dd['segm'] = dd['ra'].astype(str) + '|' + dd['sub'].astype(str)
ind = [v for v in ['internet_domiciliar', 'internet_banda_larga_fixa', 'possui_computador', 'agua_rede', 'solucao_sanitaria_estruturada', 'rua_pavimentada', 'log_renda_hh', 'possui_automovel', 'comodos_por_pessoa'] if v in dd.columns]
for y in ind:
    d1 = dd[[y, 'ra', 'sub', 'segm']].dropna()
    tot = d1[y].var(ddof=0) * len(d1)
    ss = lambda g: ((d1.groupby(g)[y].transform('mean') - d1[y].mean()) ** 2).sum()
    r_ra, r_sub, r_seg = ss('ra') / tot, ss('sub') / tot, ss('segm') / tot
    r_ad = smf.ols(f'{y} ~ C(ra) + C(sub)', d1).fit().rsquared
    m = smf.mixedlm(f'{y} ~ 1', d1, groups=d1['ra'], vc_formula={'sub': '0 + C(sub)'}).fit(reml=True) if False else None
    rows.append({'nivel': 'domicilio', 'variavel': y, 'R2_RA': r_ra, 'R2_sub': r_sub, 'R2_RA+sub': r_ad, 'R2_segmento(RAxsub)': r_seg, 'residual': 1 - r_seg})
dv = pd.DataFrame(rows); dv.to_csv('variancia_ra_sub.csv', index=False)
print('1. DECOMPOSIÇÃO DE VARIANCIA (R2 de ANOVA; ICC de modelo misto com intercepto aleatório de RA, subpopulação como efeito fixo)')
print(dv.round(3).to_string(index=False))
print('Obs.: com 49 segmentos e 24 RAs o R2 de RA é inflado (muitos parâmetros); o R2 ajustado e o ICC devem ser lidos junto.')

# ---- 2. AMI dos grupos com RA e subpopulação (segmento e domicílio), com p de permutação
lab = seg[CLc].to_numpy(); ra = seg['ra_norm'].to_numpy(); sub = seg[SUB].to_numpy()
def ami_p(a, b, n=999):
    o = AMI(a, b); nulo = [AMI(a, rng.permutation(b)) for _ in range(n)]
    return o, (np.sum(np.array(nulo) >= o) + 1) / (n + 1), float(np.mean(nulo))
out = []
for nome, v in (('RA', ra), ('subpopulação', sub)):
    o, p, mn = ami_p(lab, v); out.append({'nivel': 'segmento', 'rotulo': nome, 'AMI': o, 'AMI_nulo_medio': mn, 'p_perm': p})
mapa = seg.set_index(CHAVES)[CLc]
dd['grupo'] = pd.MultiIndex.from_frame(dd[CHAVES]).map(mapa.to_dict()) if False else dd.merge(seg[CHAVES + [CLc]], on=CHAVES, how='left')[CLc].to_numpy()
d2 = dd.dropna(subset=['grupo'])
for nome, v in (('RA', d2['ra']), ('subpopulação', d2['sub'])):
    o = AMI(d2['grupo'], v); nulo = [AMI(d2['grupo'], rng.permutation(v.to_numpy())) for _ in range(199)]
    out.append({'nivel': 'domicilio', 'rotulo': nome, 'AMI': o, 'AMI_nulo_medio': np.mean(nulo), 'p_perm': (np.sum(np.array(nulo) >= o) + 1) / 200})
am = pd.DataFrame(out); am.to_csv('ami_rotulos.csv', index=False)
print('\n2. AMI DOS GRUPOS COM RA E SUBPOPULAÇÃO'); print(am.round(3).to_string(index=False))

# ---- 3. Join count dos grupos no espaço (RA = grupo dominante ponderado por domicílios expandidos)
viz = json.load(open('vizinhos_queen.json'))
w = seg.groupby(['ra_norm', CLc])['domicilios_expandidos'].sum().reset_index()
dom = w.loc[w.groupby('ra_norm')['domicilios_expandidos'].idxmax()].set_index('ra_norm')[CLc]
ras = [r for r in dom.index if r in viz]
pares = {tuple(sorted((a, b))) for a in ras for b in viz[a] if b in dom.index}
def mesmos(l): return sum(l[a] == l[b] for a, b in pares)
obs = mesmos(dom.to_dict()); nul = []
for _ in range(9999):
    pl = dict(zip(dom.index, rng.permutation(dom.to_numpy()))); nul.append(mesmos(pl))
p_jc = (np.sum(np.array(nul) >= obs) + 1) / 10000
jc = pd.DataFrame([{'RAs': len(dom), 'vinculos_queen': len(pares), 'vinculos_mesmo_grupo': obs, 'esperado_nulo': np.mean(nul), 'p_perm': p_jc}])
jc.to_csv('joincount.csv', index=False)
print('\n3. JOIN COUNT (grupo dominante por RA; Queen; 9999 permutações de rótulos)'); print(jc.round(3).to_string(index=False))
# coesão dos grupos nos segmentos (todas as RAs, grupo de cada segmento): vizinhos RA com segmentos do mesmo grupo
print('Grupo dominante por RA:'); print(dom.to_string())
