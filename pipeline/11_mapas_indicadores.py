# -*- coding: utf-8 -*-
# Mapas dos componentes e de indicadores selecionados por RA.
# Entrada: segmentos_completo.csv, shapefile das RAs
# Saída: indicadores_por_ra.csv, figuras/mapa_componentes_ra.png, figuras/mapa_indicadores_ra.png (RAs numeradas, como no programa 03)
import numpy as np, pandas as pd, geopandas as gpd, unicodedata, re, matplotlib
import matplotlib.patches
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=pd.read_csv('segmentos_completo.csv'); w='domicilios_expandidos'
g=gpd.read_file('shape_df_extracao/regioes_administrativas.shp')   # extrair regioes_administrativas.zip em shape_df_extracao/
print(g.columns.tolist())
col=next((c for c in g.columns if c.lower() in ["ra","nome","ra_nome","nome_ra","regiao","nm_ra"]),g.columns[0])
def nrm(s):
    s=unicodedata.normalize('NFKD',str(s)).encode('ascii','ignore').decode().lower()
    return re.sub(r'\s+',' ',re.sub(r'[^a-z0-9 ]',' ',s)).strip()
def fix(x):
    try: x=x.encode('latin1').decode('utf8')
    except Exception: pass
    return x
g['k']=g[col].map(fix).map(nrm).replace({'scia':'scia estrutural','sol nascente e por do sol':'sol nascente por do sol'}); d['k']=d.ra_norm.map(nrm)
print('sem shape:',sorted(set(d.k)-set(g.k))); 
comp={'D1_Bens_Duraveis_MC_PC1':'D1 PC1 — bens duráveis','D2_Humano_Trabalho_MC_PC1':'D2 PC1 — formalização e escolaridade','D3_Fisico_Habitacional_Territorial_MC_PC1':'D3 PC1 — infraestrutura pública','D3_Fisico_Habitacional_Territorial_MC_PC2':'D3 PC2 — condição da moradia','D4_Tecnologico_Conectividade_MC_PC1':'D4 PC1 — conectividade'}
var={'renda_domiciliar_pc':'Renda per capita (R$)','escolaridade_adulto':'Escolaridade (anos)','agua_rede':'Água da rede (%)','solucao_sanitaria_estruturada':'Solução sanitária (%)','rua_pavimentada':'Rua pavimentada (%)','internet_banda_larga_fixa':'Banda larga fixa (%)','possui_computador':'Computador (%)','ivs_df':'IVS-DF'}
pct={'agua_rede','solucao_sanitaria_estruturada','rua_pavimentada','internet_banda_larga_fixa','possui_computador'}
cols=list(comp)+list(var)
def wm(x,c): 
    m=x[c].notna(); return np.average(x.loc[m,c],weights=x.loc[m,w])
ra=d.groupby('k').apply(lambda x: pd.Series({c:wm(x,c) for c in cols}))
ra['ra']=d.groupby('k').ra_norm.first()
for c in pct: ra[c]*=100
ra.to_csv('indicadores_por_ra.csv')
gg=g.merge(ra.reset_index(),on='k',how='left')

# ==============================================================================
# NUMERAÇÃO DAS RAs NOS MAPAS (os nomes ficam em tabela à parte, no texto)
# ==============================================================================
NOME_RA = {"aguas claras": "Águas Claras", "brazlandia": "Brazlândia", "ceilandia": "Ceilândia", "fercal": "Fercal", "gama": "Gama",
           "itapoa": "Itapoã", "lago norte": "Lago Norte", "nucleo bandeirante": "Núcleo Bandeirante", "paranoa": "Paranoá",
           "park way": "Park Way", "planaltina": "Planaltina", "recanto das emas": "Recanto das Emas", "riacho fundo": "Riacho Fundo",
           "riacho fundo ii": "Riacho Fundo II", "samambaia": "Samambaia", "santa maria": "Santa Maria", "sao sebastiao": "São Sebastião",
           "scia estrutural": "SCIA", "sobradinho": "Sobradinho", "sobradinho ii": "Sobradinho II",
           "sol nascente por do sol": "Sol Nascente/Pôr do Sol", "taguatinga": "Taguatinga", "vicente pires": "Vicente Pires"}

def numerar_ras(ras):
    """Número de cada RA (ordem alfabética do nome sem acento), igual em todos os mapas."""
    return {r: i + 1 for i, r in enumerate(sorted(set(ras)))}

def rotular_numeros(ax, g, num, col="ra_norm", fs=10, dmin=0.055):
    """Desenha o número de cada RA no seu ponto interno. Números que ficariam sobrepostos (RAs pequenas do
    centro do DF) são afastados uns dos outros e ligados ao ponto de origem por uma linha."""
    import numpy as np
    g = g[g[col].isin(num)]
    if g.empty: return
    x0, y0, x1, y1 = g.total_bounds
    W = max(x1 - x0, 1e-9); H = max(y1 - y0, 1e-9)
    ref = max(W, 0.6 * H) if W >= H else H
    nomes = list(g[col]); org = []
    for geom in g.geometry:
        p = geom.representative_point(); org.append((p.x, p.y))
    org = np.array(org); pos = org.copy(); d0 = dmin * ref
    rng = np.random.default_rng(1)
    for _ in range(400):
        mov = np.zeros_like(pos)
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                v = pos[i] - pos[j]; dist = np.hypot(*v)
                if dist < d0:
                    if dist < 1e-12: v = rng.normal(size=2); dist = np.hypot(*v)
                    push = (d0 - dist) / 2 * v / dist; mov[i] += push; mov[j] -= push
        pos += mov + 0.02 * (org - pos)
    for n_, o, p in zip(nomes, org, pos):
        if np.hypot(*(p - o)) > 0.5 * d0:
            ax.plot([o[0], p[0]], [o[1], p[1]], color="black", lw=.6, zorder=5)
            ax.plot(o[0], o[1], "o", color="black", ms=2.2, zorder=5)
        ax.text(p[0], p[1], str(num[n_]), fontsize=fs, ha="center", va="center", fontweight="bold", zorder=6,
                bbox=dict(boxstyle="circle,pad=0.18", fc="white", ec="black", lw=.6))

NUM_RA=numerar_ras(d.k)
gg['ra_norm']=gg['k']
def panel(items,ncol,nrow,fname,cmaps,title=None,figsize=(9.5,11)):
    fig,axes=plt.subplots(nrow,ncol,figsize=figsize)
    for ax,(c,lab) in zip(axes.ravel(),items.items()):
        gg.plot(ax=ax,color='#f0f0f0',edgecolor='#bbbbbb',linewidth=.4)
        has=gg[gg[c].notna()]
        cm='RdBu' if c in comp else 'viridis'
        kw=dict(column=c,cmap=cm,edgecolor='black',linewidth=.5,legend=True,legend_kwds={'shrink':.7})
        if c in comp:
            v=np.abs(has[c]).max(); kw.update(vmin=-v,vmax=v)
        has.plot(ax=ax,**kw); ax.set_axis_off(); ax.set_title(lab,fontsize=12,fontweight='bold')
        rotular_numeros(ax,gg,NUM_RA,col='ra_norm',fs=9,dmin=0.065)
    for ax in axes.ravel()[len(items):]: ax.set_axis_off()
    plt.tight_layout(); plt.savefig(fname,dpi=200,bbox_inches='tight'); plt.close()
panel(comp,2,3,'figuras/mapa_componentes_ra.png',None,figsize=(9.5,10.5))
panel(var,2,4,'figuras/mapa_indicadores_ra.png',None,figsize=(9.5,13.5))
print(ra.round(2).head(30).to_string())
