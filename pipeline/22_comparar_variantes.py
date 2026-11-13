# -*- coding: utf-8 -*-
"""Consolida as métricas de todas as variantes em resultados/ (sem eleger vencedor).
Saída: resultados/comparacao_variantes.csv (uma linha por variante) e resultados/comparacao_variantes.md (métricas nas linhas, variantes nas colunas).
Todas as métricas vêm dos arquivos gravados por 03, 16, 18, 19 e 21 em cada variante; nada é recalculado aqui, exceto agregações simples."""
import re, sys, glob, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import kruskal

RAIZ = Path(__file__).resolve().parent.parent; RES = RAIZ / 'resultados'
def num(s):
    try: return float(s.replace(',', '.'))
    except Exception: return np.nan

def linha_variante(d):
    cfg = dict(kv.split('=', 1) for kv in (d / 'config.txt').read_text().split())
    K = int(cfg['K']); r = {'variante': cfg['variante'], 'V2': cfg['V2'] or '(nenhuma)', 'K': K, 'WINSOR_P': cfg['WINSOR_P'], 'COMP': cfg['COMP'], 'N_BOOT': cfg['N_BOOT'], 'NINIT_BOOT': cfg['NINIT_BOOT']}
    seg = pd.read_csv(d / 'segmentos_completo.csv'); lab = seg[f'Cluster_K{K}']
    r['n_componentes'] = sum('_MC_PC' in c for c in seg.columns); r['tamanhos_grupos(seg)'] = str(sorted(lab.value_counts().tolist(), reverse=True))
    wt = seg.groupby(f'Cluster_K{K}')['domicilios_expandidos'].sum(); r['pct_dom_menor_grupo'] = float(100 * wt.min() / wt.sum())
    # --- 03: benchmark de K, nulos e estabilidade por subamostra
    t03 = (d / 'saida_03.txt').read_text(encoding='utf-8'); tab = {}; bench = {}
    ini = t03.find('Silhueta contra dois nulos')
    for ln in t03[ini:ini + 2500].splitlines():
        p = ln.split()
        if len(p) == 7 and p[0].isdigit(): tab[int(p[0])] = [num(x) for x in p[1:]]
    ini = t03.find('Benchmark de K')
    for ln in t03[ini:ini + 2500].splitlines():
        p = ln.split()
        if len(p) == 7 and p[0].isdigit(): bench[int(p[0])] = [num(x) for x in p[1:]]
    if K in tab:
        r['silhueta'] = tab[K][0]; r['p_embaralhado'] = tab[K][2]; r['p_covariancia'] = tab[K][4]; r['ARI_subamostra80'] = tab[K][5]
    if K in bench: r['CH'] = bench[K][1]; r['DB'] = bench[K][2]
    if tab:
        r['p_cov_minimo_entre_K'] = min(v[4] for v in tab.values()); r['K_do_p_cov_minimo'] = min(tab, key=lambda k: tab[k][4])
        r['n_K_com_p_cov<0.05'] = int(sum(v[4] < 0.05 for v in tab.values()))
        for k in (2, 3, 5):
            if k in tab: r[f'silhueta_K{k}'] = tab[k][0]; r[f'p_cov_K{k}'] = tab[k][4]; r[f'ARI_sub_K{k}'] = tab[k][5]
    # --- 16: bootstrap de domicílios
    t16 = (d / 'saida_16.txt').read_text(encoding='utf-8')
    m = re.search(r'ARI vs\. partição de referência: média ([\d.]+) \| mediana ([\d.]+) \| p5 ([\d.]+) \| p95 ([\d.]+) \| min ([\d.]+)', t16)
    if m: r['boot_ARI_medio'], r['boot_ARI_mediana'], r['boot_ARI_p5'], r['boot_ARI_p95'], r['boot_ARI_min'] = map(float, m.groups())
    for k, v in re.findall(r'K=(\d): ARI médio ([\d.]+) \(\d+ réplicas\)', t16): r[f'boot_ARI_K{k}'] = float(v)
    m = re.search(r'ENCOLHIMENTO.*?ARI vs\. referência ([\d.]+); mudam: (\d+)', t16)
    if m: r['encolhimento_ARI'], r['encolhimento_mudam'] = float(m.group(1)), int(m.group(2))
    m = re.search(r'K=%d: observado ([\d.]+) \| nulo média ([\d.]+) \(p95 ([\d.]+)\) \| p = ([\d.]+)' % K, t16)
    if m: r['estab_obs'], r['estab_nulo_medio'], r['estab_nulo_p95'], r['estab_p_vs_nulo'] = map(float, m.groups())
    m = re.search(r'FRONTEIRA \(coligação<0,60\): (\d+)', t16)
    if m: r['n_fronteira_boot'] = int(m.group(1))
    m = re.search(r'em fronteira \(%\): ([\d.]+)', t16)
    if m: r['pct_dom_fronteira'] = float(m.group(1))
    for ln in t16.splitlines():
        if ln.strip().startswith('K-MEANS PONDERADO') or 'integral' in ln.lower() and 'ARI' in ln:
            pass
    rb = pd.read_csv(d / 'robustez_amostral_segmentos.csv'); g = rb.groupby('grupo_ref')['junto_do_proprio_grupo'].mean()
    r['coligacao_media_do_grupo_menos_estavel'] = float(g.min()); r['coligacao_media_do_grupo_mais_estavel'] = float(g.max()); r['n_nao_classificaveis'] = int(rb['nao_classificavel'].sum())
    # --- 18: AMI, variância, join count
    am = pd.read_csv(d / 'ami_rotulos.csv'); a = am[am.nivel == 'segmento']
    for _, x in a.iterrows(): r[f"AMI_seg_{'RA' if x.rotulo == 'RA' else 'subpop'}"] = x.AMI; r[f"AMI_seg_p_{'RA' if x.rotulo == 'RA' else 'subpop'}"] = x.p_perm
    jc = pd.read_csv(d / 'joincount.csv'); r['joincount_obs'] = jc.vinculos_mesmo_grupo.iloc[0]; r['joincount_esperado'] = jc.esperado_nulo.iloc[0]; r['joincount_p'] = jc.p_perm.iloc[0]
    # --- 19: validação externa no nível de segmento (Kruskal-Wallis entre grupos)
    ve = pd.read_csv(d / 'validacao_externa_segmentos.csv'); cl = f'Cluster_K{K}'
    for v in ['ebia_insegur_qualquer', 'ebia_mod_grave', 'beneficio_social', 'cadunico', 'escritura_registrada', 'doc_formal', 'outorga_agua']:
        gs = [ve.loc[ve[cl] == g, v].dropna().to_numpy() for g in sorted(ve[cl].unique())]
        try: r[f'VE_KW_p_{v}'] = float(kruskal(*gs).pvalue)
        except Exception: r[f'VE_KW_p_{v}'] = np.nan
    if 'ivs_df' in seg.columns:
        gs = [seg.loc[lab == g, 'ivs_df'].dropna().to_numpy() for g in sorted(lab.unique())]
        r['IVS_KW_p'] = float(kruskal(*gs).pvalue)
    # --- 21: estresse
    st = pd.read_csv(d / 'stress_resumo.csv').iloc[0].to_dict()
    for k_, v_ in st.items():
        if k_ not in ('K', 'n_segmentos'): r[k_] = v_
    f6b = d / 'stress_S6b_resumo.csv'
    if f6b.exists():
        for k_, v_ in pd.read_csv(f6b).iloc[0].to_dict().items(): r[k_] = v_
    return r

linhas = []
for d in sorted(RES.iterdir()):
    if d.is_dir() and (d / 'stress_resumo.csv').exists() and (d / 'config.txt').exists():
        try: linhas.append(linha_variante(d))
        except Exception as e: print('ERRO', d.name, repr(e), file=sys.stderr)
df = pd.DataFrame(linhas)
ordem = ['M0_antigo_K5', 'M0_antigo_K3', 'M1_revisado_K3', 'M1_revisado_K5']
df['__o'] = df.variante.map(lambda v: ordem.index(v) if v in ordem else 99); df = df.sort_values(['__o', 'variante']).drop(columns='__o')
df.to_csv(RES / 'comparacao_variantes.csv', index=False)
def fmt(x):
    if isinstance(x, (float, np.floating)): return '' if np.isnan(x) else (f'{x:.3f}' if abs(x) < 1000 else f'{x:.0f}')
    return str(x)
T = df.set_index('variante').T
with open(RES / 'comparacao_variantes.md', 'w', encoding='utf-8') as f:
    f.write('# Comparação de variantes (métricas, sem eleger vencedor)\n\nGerado por `pipeline/22_comparar_variantes.py`. Colunas = variantes (ver `pipeline/variantes.tsv`); linhas = métricas.\n\n')
    f.write('| métrica | ' + ' | '.join(T.columns) + ' |\n|' + '---|' * (len(T.columns) + 1) + '\n')
    for k_, row in T.iterrows(): f.write(f'| {k_} | ' + ' | '.join(fmt(v) for v in row.tolist()) + ' |\n')
print(f'{len(df)} variantes consolidadas -> resultados/comparacao_variantes.csv e .md')
