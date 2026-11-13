# -*- coding: utf-8 -*-
# Validação externa dos grupos com variáveis que NÃO entram na matriz de agrupamento:
# insegurança alimentar (EBIA, D01-D14), benefícios sociais/CadÚnico (I24), posse do imóvel (B04, B07) e outorga de água (B14_9).
# Estimativas ponderadas por grupo com IC 95% por bootstrap de domicílios dentro dos segmentos (grupo de cada segmento fixo na referência).
# Saída: validacao_externa.csv, validacao_externa_segmentos.csv, saida_19.txt
import os, numpy as np, pandas as pd
from scipy.stats import kruskal, spearmanr
src = open(os.environ.get('TCC_03', '03_analise_cluster_4d.py'), encoding='utf-8').read()
exec(compile(src[:src.index('def monte_carlo_parallel_analysis')], '03_fase1', 'exec'))
CLc = f"Cluster_K{K_FINAL}"; B = int(os.environ.get('N_BOOT', 500)); rng = np.random.default_rng(42)
seg = pd.read_csv('segmentos_completo.csv')
dom = pd.read_excel(BASE / 'pdadr_domicilios.xlsx'); mor = pd.read_excel(BASE / 'pdadr_moradores.xlsx')
dcols = [c for c in dom.columns if c.startswith('D') and c[1:3].isdigit() and int(c[1:3]) <= 14 and 'ultimos_tres' in c]
assert len(dcols) == 14, dcols
sim = (dom[dcols] == 'Sim').astype(int); ns = (dom[dcols] == 'Não sabe').any(axis=1)
com_menor = (dom[dcols[8:]] != 'Não se aplica').any(axis=1)
score = sim.sum(axis=1)
corte_mod = np.where(com_menor, 6, 4)       # EBIA: 14 itens (com menor de 18): 6-9 moderada, 10-14 grave; 8 itens: 4-5 moderada, 6-8 grave
d = dom[[ID]].copy()
d['ebia_score'] = score; d['ebia_insegur_qualquer'] = (score >= 1).astype(float); d['ebia_mod_grave'] = (score >= corte_mod).astype(float)
d.loc[ns, ['ebia_insegur_qualquer', 'ebia_mod_grave']] = np.nan     # "não sabe" em algum item: fora do cálculo
d['escritura_registrada'] = dom['B07_a_propriedade_possui_escritura_definitiva_registrada_em_cartorio'].map({'Sim': 1.0, 'Não': 0.0})
d['doc_formal'] = dom['B04_documentacao_do_domicilio'].map(lambda x: 1.0 if x in ('Escritura/Matrícula', 'Contrato de Compra e Venda', 'Concessão de Direito de Uso ou Concessão de Direito Real de Uso') else (np.nan if x in ('Não sabe', 'Não se aplica') else 0.0))
d['outorga_agua'] = dom['B14_9_possui_outorga_para_uso_da_agua'].map({'Sim': 1.0, 'Não': 0.0})
c1 = [c for c in mor.columns if c.startswith('I24_1')][0]; c2 = [c for c in mor.columns if c.startswith('I24_2')][0]
cad = mor.groupby(ID)[c1].apply(lambda s: float((s == 'Sim').any())); ben = mor.groupby(ID)[c2].apply(lambda s: float((s == 'Sim').any()))
d['cadunico'] = d[ID].map(cad); d['beneficio_social'] = d[ID].map(ben)
base = df_dom[[ID] + CHAVES + ['peso_dom_corrigido']].merge(d, on=ID, how='left').merge(seg[CHAVES + [CLc, 'amostra_n']], on=CHAVES, how='left')
base = base.dropna(subset=[CLc]); base[CLc] = base[CLc].astype(int)
VARS = ['ebia_insegur_qualquer', 'ebia_mod_grave', 'beneficio_social', 'cadunico', 'escritura_registrada', 'doc_formal', 'outorga_agua']
K = K_FINAL
def est(df, w):
    return np.array([[np.average(df.loc[(df[CLc] == g) & df[v].notna(), v], weights=df.loc[(df[CLc] == g) & df[v].notna(), w]) if ((df[CLc] == g) & df[v].notna()).any() else np.nan for v in VARS] for g in range(K)])
P0 = est(base, 'peso_dom_corrigido')
seg_id = base.groupby(CHAVES).ngroup().to_numpy(); idxs = [np.where(seg_id == s)[0] for s in range(seg_id.max() + 1)]
boot = []
for _ in range(B):
    m = np.zeros(len(base))
    for ix in idxs: m[ix] = np.bincount(rng.integers(0, len(ix), len(ix)), minlength=len(ix))
    b = base.copy(); b['wb'] = b['peso_dom_corrigido'] * m; boot.append(est(b, 'wb'))
boot = np.array(boot)
rows = [{'grupo': g, 'variavel': v, 'estimativa_pct': 100 * P0[g, j], 'ic_inf': 100 * np.nanpercentile(boot[:, g, j], 2.5), 'ic_sup': 100 * np.nanpercentile(boot[:, g, j], 97.5)} for g in range(K) for j, v in enumerate(VARS)]
out = pd.DataFrame(rows); out.to_csv('validacao_externa.csv', index=False)
print('VALIDAÇÃO EXTERNA — % dos domicílios (ponderado) por grupo, IC 95% por bootstrap')
for j, v in enumerate(VARS):
    print(f'  {v}: ' + ' | '.join(f'G{g} {100*P0[g,j]:.1f} [{100*np.nanpercentile(boot[:,g,j],2.5):.1f}; {100*np.nanpercentile(boot[:,g,j],97.5):.1f}]' for g in range(K)))
print('\nEBIA: prevalência global de insegurança alimentar (qualquer):', round(100 * np.average(base['ebia_insegur_qualquer'].dropna(), weights=base.loc[base['ebia_insegur_qualquer'].notna(), 'peso_dom_corrigido']), 1), '%')
# nível segmento: associação entre escores/grupo e as variáveis externas (Kruskal-Wallis; Spearman com os componentes)
sg = base.groupby(CHAVES).apply(lambda g: pd.Series({v: np.average(g.loc[g[v].notna(), v], weights=g.loc[g[v].notna(), 'peso_dom_corrigido']) if g[v].notna().any() else np.nan for v in VARS})).reset_index()
sg = sg.merge(seg[CHAVES + [CLc, 'amostra_n'] + [c for c in seg.columns if '_MC_PC' in c]], on=CHAVES)
sg.to_csv('validacao_externa_segmentos.csv', index=False)
print('\nNível segmento (n=49): Kruskal-Wallis entre grupos e Spearman com os escores')
cp = [c for c in seg.columns if '_MC_PC' in c]
for v in VARS:
    s = sg.dropna(subset=[v]); H, p = kruskal(*[s.loc[s[CLc] == g, v] for g in range(K) if (s[CLc] == g).sum() > 1])
    print(f'  {v}: KW H={H:.2f} p={p:.3f} | ' + ' | '.join(f"{c.split('_MC_')[0][:2]}{c[-3:]} rho={spearmanr(s[v], s[c])[0]:.2f}" for c in cp))
