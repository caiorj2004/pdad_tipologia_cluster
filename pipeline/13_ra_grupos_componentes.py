# -*- coding: utf-8 -*-
# Grupos dentro das RAs: quadro das RAs com mais de um segmento, mapa de heterogeneidade intra-RA e mapas de calor dos
# escores dos componentes (RAs com mais de um grupo; segmentos de agricultura familiar).
# Entrada: segmentos_completo.csv, km_ward_labels.csv, shapefile das RAs. Saída: Tabela_ra_segmentos.tex, figuras/*.png, ra_componentes.csv
import numpy as np, pandas as pd, geopandas as gpd, unicodedata, re, matplotlib
import matplotlib.patches as mpatches
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os; os.makedirs('figuras', exist_ok=True)
d = pd.read_csv('segmentos_completo.csv'); w = 'domicilios_expandidos'
lab = pd.read_csv('km_ward_labels.csv')[['a01ra', 'sub', 'ward', 'Classe']]
SUBN = {'Agricultura Empresarial': 'Agr. Empresarial', 'Agricultura Familiar': 'Agr. Familiar', 'Assentamento/Agrovila': 'Assentamento', 'Características Urbanas': 'Moradia e Lazer'}
d['sub'] = d.a01subpopulacao.map(SUBN); d = d.merge(lab, left_on=['a01ra', 'sub'], right_on=['a01ra', 'sub'], how='left')
assert len(d) == 49 and d.Classe.notna().all()
d['G'] = d.Cluster_K3
COMP = {'D1_Bens_Duraveis_MC_PC1': 'D1 PC1\nBens duráveis', 'D2_Humano_Trabalho_MC_PC1': 'D2 PC1\nFormalização e\nescolaridade', 'D3_Fisico_Habitacional_Territorial_MC_PC1': 'D3 PC1\nInfraestrutura\npública',
        'D3_Fisico_Habitacional_Territorial_MC_PC2': 'D3 PC2\nCondição interna\nda moradia', 'D4_Tecnologico_Conectividade_MC_PC1': 'D4 PC1\nConectividade'}
NOME_G = {0: 'Moradia e saneamento bons, água da rede rara', 1: 'Infraestrutura intermediária, conectividade abaixo da média', 2: 'Alta escolaridade e conectividade (pequeno)',
          3: 'Melhor infraestrutura pública', 4: 'Maior carência de infraestrutura'}
CORES_CL = {0: "#440154", 1: "#3b528b", 2: "#21918c", 3: "#5ec962", 4: "#fde725"}
NOME_RA = {"aguas claras": "Águas Claras", "brazlandia": "Brazlândia", "ceilandia": "Ceilândia", "fercal": "Fercal", "gama": "Gama", "itapoa": "Itapoã", "lago norte": "Lago Norte",
           "nucleo bandeirante": "Núcleo Bandeirante", "paranoa": "Paranoá", "park way": "Park Way", "planaltina": "Planaltina", "recanto das emas": "Recanto das Emas", "riacho fundo": "Riacho Fundo",
           "riacho fundo ii": "Riacho Fundo II", "samambaia": "Samambaia", "santa maria": "Santa Maria", "sao sebastiao": "São Sebastião", "scia estrutural": "SCIA", "sobradinho": "Sobradinho",
           "sobradinho ii": "Sobradinho II", "sol nascente por do sol": "Sol Nascente/Pôr do Sol", "taguatinga": "Taguatinga", "vicente pires": "Vicente Pires"}
d['RA'] = d.a01ra
cl = list(COMP)
# ---- 1) quadro das RAs com mais de um segmento
nseg = d.groupby('a01ra').size(); ng = d.groupby('a01ra').G.nunique()
multi = nseg[nseg >= 2].index
rows = []
ABREV = {'Agr. Empresarial': 'AE', 'Agr. Familiar': 'AF', 'Assentamento': 'AS', 'Moradia e Lazer': 'ML'}
for ra in sorted(multi, key=lambda r: (-ng[r], -nseg[r], r)):
    x = d[d.a01ra == ra].sort_values('sub')
    seg = '; '.join(f"{ABREV[s]} ({chr(65+int(g))})" for s, g in zip(x['sub'], x.G))
    rows.append((ra, int(nseg[ra]), int(ng[ra]), seg))
tq = pd.DataFrame(rows, columns=['RA', 'Segmentos', 'Grupos', 'Detalhe']); tq.to_csv('ra_segmentos_grupos.csv', index=False)
tex = ["\\begin{table}[htbp]", "\\centering", "\\caption{RAs com mais de um segmento: grupo de cada segmento}", "\\label{tab:ra_segmentos}", "\\small",
       "\\begin{tabular}{l c c p{5.6cm}}", "\\hline", "\\textbf{RA} & \\textbf{Segmentos} & \\textbf{Grupos} & \\textbf{Grupo de cada segmento} \\\\ \\hline"]
prev = None
for r in tq.itertuples():
    if prev is not None and r.Grupos != prev: tex.append("\\hline")
    prev = r.Grupos
    tex.append(f"{r.RA} & {r.Segmentos} & {r.Grupos} & {r.Detalhe} \\\\")
tex += ["\\hline", "\\end{tabular}",
        "\\par\\footnotesize Fonte: elaboração própria a partir da PDAD Rural 2022 (IPEDF). AF: Agricultura Familiar; AE: Agricultura Empresarial; AS: Assentamento/Agrovila; ML: Moradia e Lazer. G0 a G4: grupos do $K$-Means ($K = 5$). As 7 RAs com um único segmento não aparecem.",
        "\\end{table}"]
open('Tabela_ra_segmentos.tex', 'w', encoding='utf-8').write('\n'.join(tex))
# ---- 2) mapas de calor dos escores
def heat(x, fname, titulo, figw=9.2, sep=True):
    x = x.reset_index(drop=True); M = x[cl].values
    fig, ax = plt.subplots(figsize=(figw, 0.42 * len(x) + 2.2))
    im = ax.imshow(M, aspect='auto', cmap='RdBu_r', vmin=-4, vmax=4)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:+.1f}", ha='center', va='center', fontsize=8.5, color='white' if abs(M[i, j]) > 2.6 else 'black')
    ax.set_xticks(range(len(cl))); ax.set_xticklabels([COMP[c] for c in cl], fontsize=8.5); ax.xaxis.tick_top()
    ax.set_yticks(range(len(x))); ax.set_yticklabels([f"{r} — {s} ({chr(65+int(g))}){'*' if n < 10 else ''}" for r, s, g, n in zip(x.RA, x['sub'], x.G, x.amostra_n)], fontsize=9)
    for i, g in enumerate(x.G):
        ax.add_patch(mpatches.Rectangle((-0.62, i - 0.5), 0.12, 1, color=CORES_CL[int(g)], clip_on=False, transform=ax.transData))
    for i in range(1, len(x)):
        if (sep and x.RA[i] != x.RA[i - 1]) or (not sep and x.G[i] != x.G[i - 1]): ax.axhline(i - 0.5, color='black', lw=1.2)
    ax.set_xticks(np.arange(-.5, len(cl)), minor=True); ax.tick_params(which='minor', length=0); ax.tick_params(length=0); ax.tick_params(axis='y', pad=16)
    for sp in ax.spines.values(): sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.035, pad=0.02, aspect=40); cb.set_label('Escore do componente (positivo: acima da média dos 49 segmentos; negativo: abaixo). * = menos de 10 domicílios na amostra', fontsize=8.5); cb.ax.tick_params(labelsize=8)
    plt.savefig(f'figuras/{fname}.png', dpi=200, bbox_inches='tight'); plt.close()
mg = d[d.a01ra.isin(ng[ng >= 2].index)].copy()
mg['ord'] = mg.a01ra.map(lambda r: (-ng[r], -nseg[r], r)); mg = mg.sort_values(['ord', 'sub']).drop(columns='ord')
mg['ordn'] = mg.groupby('a01ra').ngroup()
mg[['a01ra', 'sub', 'G'] + cl].to_csv('ra_componentes_multigrupo.csv', index=False)
heat(mg, 'heatmap_componentes_ras_multigrupo', '')
af = d[d['sub'] == 'Agr. Familiar'].sort_values(['G', 'a01ra']).copy()
af[['a01ra', 'G'] + cl].to_csv('af_componentes.csv', index=False)
heat(af, 'heatmap_componentes_af', '', sep=False)
print(tq.to_string()); print(len(mg), 'segmentos em', mg.a01ra.nunique(), 'RAs')
print(mg[['a01ra', 'sub', 'G'] + cl].round(2).to_string())
print(af[['a01ra', 'G'] + cl].round(2).to_string())
# ---- 3) mapa de heterogeneidade intra-RA (legenda fora do mapa)
g = gpd.read_file('shape_df_extracao/regioes_administrativas.shp')
col = next((c for c in g.columns if c.lower() in ["ra", "nome", "ra_nome", "nome_ra", "regiao", "nm_ra"]), g.columns[0])
def nrm(s):
    s = unicodedata.normalize('NFKD', str(s)).encode('ascii', 'ignore').decode().lower()
    return re.sub(r'\s+', ' ', re.sub(r'[^a-z0-9 ]', ' ', s)).strip()
def fix(x):
    try: x = x.encode('latin1').decode('utf8')
    except Exception: pass
    return x
g['k'] = g[col].map(fix).map(nrm).replace({'scia': 'scia estrutural', 'sol nascente e por do sol': 'sol nascente por do sol'})
d['k'] = d.ra_norm.map(nrm)
het = d.groupby('k').agg(n_grupos=('G', 'nunique'), n_segmentos=('G', 'size')).reset_index()
gg = g.merge(het, on='k', how='inner')
exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'helper_maps.txt'), encoding='utf-8').read())
NUM_RA = numerar_ras(d.k)
gg['ra_norm'] = gg['k']
fig, ax = plt.subplots(figsize=(11, 7))
pal = {1: '#fef0d9', 2: '#fc8d59', 3: '#990000'}
gg.plot(ax=ax, color=[pal[int(v)] for v in gg.n_grupos], edgecolor='black', linewidth=.8)
ax.set_axis_off(); rotular_numeros(ax, gg, NUM_RA, fs=11, dmin=0.045)
ax.legend(handles=[mpatches.Patch(facecolor=pal[i], edgecolor='black', label=f'{i} grupo' + ('s' if i > 1 else '')) for i in (1, 2, 3)], title='Grupos distintos\nentre os segmentos da RA', loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=12, title_fontsize=12, frameon=False)
plt.savefig('figuras/mapa_heterogeneidade_ra.png', dpi=200, bbox_inches='tight'); plt.close()
