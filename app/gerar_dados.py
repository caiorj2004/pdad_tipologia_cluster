#!/usr/bin/env python3
"""Gera app/dados.json a partir dos CSVs versionados em resultados/ e do shapefile das RAs.

Uso (na raiz do repositório):  python app/gerar_dados.py
O JSON é embutido em docs/index.html por app/montar_site.py.
"""
import json, re, unicodedata
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
from sklearn.metrics import adjusted_rand_score as ari

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / 'resultados'
M2 = RES / 'M2_semEscadas_K3'

def norm(s):
    s = unicodedata.normalize('NFKD', str(s)).encode('ascii', 'ignore').decode().lower()
    return re.sub(r'[^a-z0-9]+', ' ', s).strip()

GRUPOS = {0: ('A', 'Urbanização formal consolidada'), 1: ('B', 'Saneamento e vias precários'), 2: ('C', 'Urbanização parcial')}
SUB_CURTO = {'Agricultura Empresarial': 'Agr. Empresarial', 'Agricultura Familiar': 'Agr. Familiar',
             'Assentamento/Agrovila': 'Assentamento', 'Características Urbanas': 'Moradia e Lazer'}

# ---------------------------------------------------------------- segmentos (M2, K=3)
seg = pd.read_csv(M2 / 'segmentos_completo.csv')
rob = pd.read_csv(M2 / 'robustez_amostral_segmentos.csv')
rob['chave'] = rob['segmento']
seg['chave'] = seg['a01ra'] + ' | ' + seg['a01subpopulacao']
seg = seg.merge(rob[['chave', 'junto_do_proprio_grupo', 'grupo_encolhido']], on='chave', how='left')

def situacao(r):
    if r['amostra_n'] < 10: return 'NC'
    return 'F' if r['junto_do_proprio_grupo'] < 0.60 else 'E'
seg['situacao'] = seg.apply(situacao, axis=1)

segmentos = []
for _, r in seg.sort_values(['a01ra', 'a01subpopulacao']).iterrows():
    g = int(r['Cluster_K3'])
    segmentos.append(dict(
        ra=r['a01ra'], sub=SUB_CURTO[r['a01subpopulacao']], grupo=GRUPOS[g][0], n=int(r['amostra_n']), dom=round(float(r['domicilios_expandidos'])),
        renda=round(float(r['renda_domiciliar_pc'])), esc=round(float(r['escolaridade_adulto']), 1), agua=round(float(r['agua_rede']) * 100),
        esg=round(float(r['solucao_sanitaria_estruturada']) * 100), rua=round(float(r['rua_pavimentada']) * 100),
        banda=round(float(r['internet_banda_larga_fixa']) * 100), colig=round(float(r['junto_do_proprio_grupo']), 2), sit=r['situacao']))

# ---------------------------------------------------------------- geometria das RAs
g = gpd.read_file(RAIZ / 'dados' / 'regioes_administrativas.zip')
def conserta(s):
    try: return s.encode('latin-1').decode('utf-8')
    except Exception: return s
g['nome'] = g['ra_nome'].map(conserta)
g['chave'] = g['nome'].map(norm)
alias = {'sol nascente e por do sol': 'sol nascente por do sol', 'nucleo bandeirante': 'nucleo bandeirante'}
ras_dados = {norm(r): r for r in seg['a01ra'].unique()}
ras_dados['sol nascente por do sol'] = 'Sol Nascente/Pôr do Sol'
g['chave'] = g['chave'].replace(alias)
g['rd'] = g['chave'].map(ras_dados)
faltam = set(ras_dados.values()) - set(g['rd'].dropna())
assert not faltam, faltam

g = g.to_crs(31983)
g['geometry'] = g.geometry.simplify(120, preserve_topology=True)
minx, miny, maxx, maxy = g.total_bounds
W = 1000.0; esc_ = W / (maxx - minx); H = (maxy - miny) * esc_

def path(geom):
    polys = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
    d = []
    for p in polys:
        for ring in [p.exterior] + list(p.interiors):
            pts = [f'{(x - minx) * esc_:.1f},{(maxy - y) * esc_:.1f}' for x, y in ring.coords]
            d.append('M' + 'L'.join(pts) + 'Z')
    return ''.join(d)

# agregados por RA
por_ra = {}
for s in segmentos: por_ra.setdefault(s['ra'], []).append(s)
ras = []
for _, r in g.iterrows():
    nome_rd = r['rd'] if isinstance(r['rd'], str) else None
    item = dict(nome=nome_rd or r['nome'].title(), d=path(r.geometry), pdad=nome_rd is not None)
    if nome_rd:
        ss = por_ra[nome_rd]; tot = {}
        for s in ss: tot[s['grupo']] = tot.get(s['grupo'], 0) + s['dom']
        item['pred'] = max(tot, key=tot.get); item['grupos'] = sorted(tot); item['segmentos'] = ss
    ras.append(item)

# ---------------------------------------------------------------- grupos (perfil, validação, escores)
pf = pd.read_csv(M2 / 'perfis_grupos_ic.csv'); ve = pd.read_csv(M2 / 'validacao_externa.csv'); es = pd.read_csv(M2 / 'escores_medios_grupo.csv')
tot_dom = seg['domicilios_expandidos'].sum()
grupos = []
for k in (0, 1, 2):
    sg = seg[seg['Cluster_K3'] == k]
    grupos.append(dict(id=GRUPOS[k][0], nome=GRUPOS[k][1], n_seg=int(len(sg)), pct_dom=round(100 * sg['domicilios_expandidos'].sum() / tot_dom, 1),
                       dom=round(float(sg['domicilios_expandidos'].sum())),
                       perfil={r['variavel']: [round(r['estimativa'], 3), round(r['ic_inf'], 3), round(r['ic_sup'], 3)] for _, r in pf[pf['grupo'] == k].iterrows()},
                       valid={r['variavel']: [round(r['estimativa_pct'], 1), round(r['ic_inf'], 1), round(r['ic_sup'], 1)] for _, r in ve[ve['grupo'] == k].iterrows()},
                       escores={c: round(float(es.loc[es['Cluster_K3'] == k, c].iloc[0]), 2) for c in es.columns if c.startswith('D')}))
kw = pd.read_csv(M2 / 'kw_todas.csv')

# ---------------------------------------------------------------- comparação dos modelos
cv = pd.read_csv(RES / 'comparacao_variantes.csv').set_index('variante')
def ve_sig(r):
    return int(sum(r[c] < 0.05 for c in cv.columns if c.startswith('VE_KW_p_')))
def linha(v, rotulo, nvars):
    r = cv.loc[v]
    return dict(id=v, rotulo=rotulo, vars=nvars, comp=int(r['n_componentes']), tamanhos=str(r['tamanhos_grupos(seg)']), menor=round(float(r['pct_dom_menor_grupo']), 1),
                silh=round(float(r['silhueta']), 3), p_cov=round(float(r['p_covariancia']), 3), ari_sub=round(float(r['ARI_subamostra80']), 2),
                boot=round(float(r['boot_ARI_obs_medio']), 2), boot_nulo=round(float(r['boot_ARI_nulo_medio']), 2), boot_p=round(float(r['boot_ARI_p_vs_nulo']), 2),
                front=round(float(r['pct_dom_fronteira']), 1), sem_ra=round(float(r['S6b_ARI_medio_sem_uma_RA']), 2), sem_sub=round(float(r['S6b_ARI_medio_sem_uma_subpop']), 2),
                spear=round(float(r['S2_spearman_n_dist_Kfinal']), 2), ve=ve_sig(r), ord=int(r['ordenada_n_p<0.05']), cond=int(r['cond29_n_desfechos_p_menor_005']))
modelos = [linha('M0_antigo_K3', 'M0, K=3', 34), linha('M0_antigo_K5', 'M0, K=5', 34), linha('M1_revisado_K3', 'M1, K=3', 30),
           linha('M2_semEscadas_K3', 'M2, K=3', 30), linha('M2b_aguaILR_K3', 'M2b, K=3', 30)]
extras = [linha(f'{m}_K{k}', f'{m.split("_")[0]}, K={k}', 34 if m.startswith('M0') else 30)
          for m in ('M0_antigo', 'M1_revisado', 'M2_semEscadas', 'M2b_aguaILR') for k in (3, 4, 5)]

# ARI entre partições
def rot(v):
    k = int(re.search(r'K(\d)$', v).group(1)); d = pd.read_csv(RES / v / 'segmentos_completo.csv')
    d['chave'] = d['a01ra'] + ' | ' + d['a01subpopulacao']; return d.set_index('chave')[f'Cluster_K{k}']
cols = ['M0_antigo_K3', 'M0_antigo_K5', 'M1_revisado_K3', 'M2_semEscadas_K3', 'M2b_aguaILR_K3']
R = {c: rot(c) for c in cols}
ari_m = [[round(float(ari(R[a].sort_index(), R[b].sort_index())), 2) for b in cols] for a in cols]

# ---------------------------------------------------------------- estabilidade
nulo = pd.read_csv(M2 / 'nulos_24_resumo.csv').iloc[0]
cov = pd.read_csv(M2 / 'estabilidade_vs_nulo_cov.csv')
km = pd.read_csv(M2 / 'comparacao_metodos.csv'); km = km[km['metodo'] == 'K-Means']
s6 = pd.read_csv(M2 / 'stress_S6b_deixa_uma_RA_refaz_pca.csv').sort_values('ARI_vs_ref_restantes')
sub6 = pd.read_csv(M2 / 'stress_S6b_deixa_uma_subpop_refaz_pca.csv')
estab = dict(
    boot=dict(obs=round(float(nulo['boot_ARI_obs_medio']), 3), nulo=round(float(nulo['boot_ARI_nulo_medio']), 3), nulo_p95=round(float(nulo['boot_ARI_nulo_p95']), 3), p=round(float(nulo['boot_ARI_p_vs_nulo']), 2)),
    metades=dict(obs=round(float(nulo['metades_ARI_obs_medio']), 2), nulo=round(float(nulo['metades_ARI_nulo_medio']), 2), nulo_p95=round(float(nulo['metades_ARI_nulo_p95']), 2)),
    cov=[dict(K=int(r['K']), obs=round(r['ARI_obs'], 2), nulo=round(r['ARI_nulo_medio'], 2), p95=round(r['ARI_nulo_p95'], 2), p=round(r['p'], 2)) for _, r in cov.iterrows()],
    kmeans=[dict(K=int(r['K']), silh=round(r['silhueta'], 3), est=round(r['estab_ARI'], 2), menor=int(r['menor'])) for _, r in km.iterrows()],
    sem_ra=[dict(ra=r['a01ra'], ari=round(r['ARI_vs_ref_restantes'], 2)) for _, r in s6.iterrows()],
    sem_sub=[dict(sub=SUB_CURTO.get(r.iloc[0], r.iloc[0]), ari=round(r['ARI_vs_ref_restantes'], 2)) for _, r in sub6.iterrows()],
    colig=[dict(seg=f"{s['ra']} — {s['sub']}", c=s['colig'], g=s['grupo'], sit=s['sit'], n=s['n']) for s in sorted(segmentos, key=lambda x: x['colig'])])

out = dict(grupos=grupos, segmentos=segmentos, ras=ras, mapa=dict(w=W, h=round(H, 1)), modelos=modelos, extras=extras, ari=dict(cols=cols, m=ari_m), estab=estab,
           totais=dict(dom=round(float(tot_dom)), amostra=int(seg['amostra_n'].sum()), n_seg=int(len(seg)), n_ra=int(seg['a01ra'].nunique())))
(RAIZ / 'app').mkdir(exist_ok=True)
(RAIZ / 'app' / 'dados.json').write_text(json.dumps(out, ensure_ascii=False, separators=(',', ':')), encoding='utf-8')
print('ok:', len(ras), 'RAs,', len(segmentos), 'segmentos,', round((RAIZ / 'app' / 'dados.json').stat().st_size / 1024), 'KB')
