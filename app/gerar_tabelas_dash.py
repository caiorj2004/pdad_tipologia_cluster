#!/usr/bin/env python3
"""Gera as tabelas CSV usadas pelas páginas /tipologia e /pdad-rural do dashboard AgroTec.

Duas fontes, as mesmas do app estático e da análise exploratória do TCC:
  1. app/dados.json (resultados do modelo M2, K=3, e da robustez), que vem de app/gerar_dados.py;
  2. a análise exploratória ponderada (pipeline/02_analise_exploratoria_ponderada.py), executada
     sobre os microdados em dados/gold/.

Uso (na raiz do repositório):
    python app/gerar_dados.py                      # se app/dados.json ainda não existir
    python app/gerar_tabelas_dash.py --saida CAMINHO/dados_tcc

Todas as tabelas estão em nível agregado (segmento, RA, grupo ou subpopulação). Nenhuma linha
identifica domicílio ou morador.
"""
import argparse, contextlib, io, json, os, runpy, sys, tempfile
from pathlib import Path
import numpy as np, pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
ap = argparse.ArgumentParser()
ap.add_argument('--saida', required=True, help='pasta de saída (será criada)')
args = ap.parse_args()
OUT = Path(args.saida); OUT.mkdir(parents=True, exist_ok=True)

def salva(df, nome, **kw):
    df.to_csv(OUT / nome, index=False, encoding='utf-8', **kw)
    print(f'  {nome}: {len(df)} linhas')

# ------------------------------------------------------------------ 1. tipologia (M2, K=3)
D = json.load(open(RAIZ / 'app' / 'dados.json', encoding='utf-8'))
print('Tipologia (app/dados.json)')

g_rows, p_rows = [], []
for g in D['grupos']:
    g_rows.append(dict(grupo=g['id'], nome=g['nome'], n_seg=g['n_seg'], pct_dom=g['pct_dom'], dom=g['dom']))
    for bloco in ('perfil', 'valid'):
        for var, (est, lo, hi) in g[bloco].items():
            p_rows.append(dict(grupo=g['id'], bloco=bloco, variavel=var, estimativa=est, ic_inf=lo, ic_sup=hi))
    for var, v in g['escores'].items():
        p_rows.append(dict(grupo=g['id'], bloco='escore', variavel=var, estimativa=v, ic_inf=np.nan, ic_sup=np.nan))
salva(pd.DataFrame(g_rows), 'tipologia_grupos.csv')
salva(pd.DataFrame(p_rows), 'tipologia_perfil_grupos.csv')

seg = pd.DataFrame(D['segmentos']).rename(columns={'sub': 'subpopulacao', 'dom': 'dom_expandidos', 'esc': 'escolaridade',
                                                   'renda': 'renda_pc', 'esg': 'saneamento', 'rua': 'pavimentacao',
                                                   'banda': 'banda_larga', 'colig': 'coligacao', 'sit': 'situacao', 'n': 'amostra_n'})
salva(seg, 'tipologia_segmentos.csv')

# polígonos das RAs como caminhos SVG (x, y em unidades do mapa; origem no canto superior esquerdo)
ras = pd.DataFrame([dict(ra=r['nome'], caminho=r['d'], pdad=bool(r.get('pdad')), grupo_predominante=r.get('pred') or '',
                         grupos=','.join(r.get('grupos') or [])) for r in D['ras']])
salva(ras, 'tipologia_mapa_ras.csv')
pd.DataFrame([dict(largura=D['mapa']['w'], altura=D['mapa']['h'])]).to_csv(OUT / 'tipologia_mapa_dim.csv', index=False)

mod = pd.DataFrame(D['modelos'] + [e for e in D['extras'] if e['id'] not in {m['id'] for m in D['modelos']}])
salva(mod, 'tipologia_modelos.csv')
ari = pd.DataFrame(D['ari']['m'], index=D['ari']['cols'], columns=D['ari']['cols']).reset_index(names='modelo')
salva(ari, 'tipologia_ari_modelos.csv')

E = D['estab']
est = [dict(teste='bootstrap_domicilios', **E['boot']), dict(teste='metades_territoriais', **E['metades'])]
salva(pd.DataFrame(est), 'tipologia_estabilidade_resumo.csv')
salva(pd.DataFrame(E['cov']), 'tipologia_estabilidade_nulo_cov.csv')
salva(pd.DataFrame(E['kmeans']), 'tipologia_escolha_k.csv')
salva(pd.DataFrame(E['sem_ra']), 'tipologia_sem_uma_ra.csv')
salva(pd.DataFrame(E['sem_sub']), 'tipologia_sem_uma_subpopulacao.csv')
salva(pd.DataFrame(E['colig']), 'tipologia_coligacao_segmentos.csv')
salva(pd.DataFrame([D['totais']]), 'tipologia_totais.csv')

# ------------------------------------------------------------------ 2. análise exploratória (PDAD Rural)
print('Análise exploratória ponderada (pipeline/02)')
tmp = tempfile.mkdtemp()
os.environ.setdefault('PDADR_DOM', str(RAIZ / 'dados' / 'gold' / 'pdadr_domicilios.xlsx'))
os.environ.setdefault('PDADR_MOR', str(RAIZ / 'dados' / 'gold' / 'pdadr_moradores.xlsx'))
os.environ['PDADR_FIG_DIR'] = tmp
import matplotlib; matplotlib.use('Agg')
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(RAIZ / 'pipeline' / '02_analise_exploratoria_ponderada.py'))
dom, mor, W, WM, SUB, RA = ns['dom'], ns['mor'], ns['W'], ns['WM'], ns['SUB'], ns['RA']
wmean, wq = ns['wmean'], ns['wquantile']

renomeia = lambda t, c: t.reset_index().rename(columns={t.index.name or 'index': c})
salva(renomeia(ns['t_sub'], 'subpopulacao'), 'pdad_subpopulacao.csv')
salva(renomeia(ns['t_ra'], 'ra'), 'pdad_ra.csv')
salva(renomeia(ns['t_fx'], 'faixa_etaria'), 'pdad_faixa_etaria.csv')
pir = ns['pt'].reset_index().rename(columns={'faixa': 'faixa_etaria'})
pir.columns.name = None
salva(pir, 'pdad_piramide.csv')

linhas = []
for nome, df, col, wc in (('Renda domiciliar total (R$)', dom, 'renda_domiciliar_total', W), ('Renda domiciliar per capita (R$)', dom, 'renda_domiciliar_pc', W),
                          ('Renda individual (R$)', mor, 'renda_individual', WM)):
    q = wq(df[col], df[wc], [.25, .5, .75, .9])
    linhas.append(dict(variavel=nome, n_validos=int(df[col].notna().sum()), media_pond=wmean(df[col], df[wc]), q1=q[0], mediana=q[1], q3=q[2], p90=q[3]))
salva(pd.DataFrame(linhas), 'pdad_renda_resumo.csv')

def renda_por(df, col, wc):
    d = df.dropna(subset=[col, 'renda_individual'])
    out = []
    for k, x in d.groupby(col, observed=True):
        out.append(dict(categoria=str(k), n=len(x), media=wmean(x['renda_individual'], x[wc]), mediana=wq(x['renda_individual'], x[wc], .5)[0]))
    return pd.DataFrame(out)
re_ = renda_por(mor, 'escolaridade_cat', WM)
re_['ordem'] = re_['categoria'].map({k: i for i, k in enumerate(list(ns['MAPA_ANOS']))})
salva(re_.sort_values('ordem').drop(columns='ordem'), 'pdad_renda_escolaridade.csv')
rb = renda_por(mor, 'recebe_beneficios_sociais', WM)
rb['categoria'] = rb['categoria'].map(lambda v: 'Recebe benefício' if str(v) in ('1', '1.0', 'True') else 'Não recebe benefício')
salva(rb, 'pdad_renda_beneficio.csv')

vd = ['renda_domiciliar_pc', 'agua_rede', 'solucao_sanitaria_estruturada', 'rua_pavimentada', 'internet_domiciliar',
      'internet_banda_larga_fixa', 'possui_computador', 'possui_automovel']
por_sub = pd.DataFrame({s: {v: wmean(d[v], d[W]) for v in vd} for s, d in dom.groupby(SUB)}).T
por_sub.index.name = 'subpopulacao'
vm = ['escolaridade_adulto', 'trabalhou_30d', 'tem_ctps', 'contribui_previdencia', 'recebe_beneficios_sociais', 'internet_uso_3m']
por_sub_m = pd.DataFrame({s: {v: wmean(d[v], d[WM]) for v in vm} for s, d in mor.groupby(SUB)}).T
salva(por_sub.join(por_sub_m).reset_index(), 'pdad_indicadores_subpopulacao.csv')

n_ra = dom[RA].nunique()
tot = dict(domicilios_amostra=len(dom), moradores_amostra=len(mor), domicilios_expandidos=float(dom[W].sum()),
           ras=int(n_ra), idade_media_pond=float(wmean(mor['idade_num'], mor[WM])),
           renda_pc_mediana_pond=float(wq(dom['renda_domiciliar_pc'], dom[W], .5)[0]),
           pct_renda_dom_abaixo_2000=float(dom.loc[dom['renda_domiciliar_total'] < 2000, W].sum() / dom.loc[dom['renda_domiciliar_total'].notna(), W].sum() * 100))
salva(pd.DataFrame([tot]), 'pdad_totais.csv')

seg_p = pd.read_csv(Path(tmp) / 'segmentos_ponderados.csv')
cols = ['a01ra', 'a01subpopulacao', 'amostra_n_dom', 'dom_expandidos', 'renda_domiciliar_pc', 'escolaridade_adulto', 'agua_rede',
        'solucao_sanitaria_estruturada', 'rua_pavimentada', 'internet_banda_larga_fixa', 'possui_computador', 'recebe_beneficios_sociais']
cols = [c for c in cols if c in seg_p.columns]
salva(seg_p[cols].rename(columns={'a01ra': 'ra', 'a01subpopulacao': 'subpopulacao'}), 'pdad_segmentos.csv')

# correlações de Spearman entre os 49 segmentos
vs = [c for c in ['renda_domiciliar_pc', 'escolaridade_adulto', 'agua_rede', 'solucao_sanitaria_estruturada', 'rua_pavimentada',
                  'internet_banda_larga_fixa', 'possui_computador', 'recebe_beneficios_sociais'] if c in seg_p.columns]
sp = seg_p[vs].corr(method='spearman')
sp.index.name = 'variavel'
salva(sp.reset_index(), 'pdad_spearman_segmentos.csv')
print('Concluído:', OUT)
