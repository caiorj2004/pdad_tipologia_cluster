# -*- coding: utf-8 -*-
# Gera as tabelas LaTeX de escores, Kruskal-Wallis e apoio, e as figuras de escores e de apoio por RA.
# Entrada: segmentos_completo.csv, kw_todas.csv, apoio_*.csv
# Saída: tab_escores.tex, Tabela_kw_variaveis.tex, tab_apoio_*.tex, figuras/escores_grupos.png, figuras/mapa_apoio_ra.png
import numpy as np, pandas as pd, geopandas as gpd, unicodedata, re, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
def f1(x): return '--' if pd.isna(x) else f"{x:.1f}".replace('.',',')
def f2(x):
    t=f"{x:.2f}"
    if t=="-0.00": t="0.00"
    return t.replace('.',',')
def fp(p): return '$<0{,}001$' if p<0.001 else f"{p:.3f}".replace('.',',')
# ---------- escores por grupo
d=pd.read_csv('segmentos_completo.csv'); kw=pd.read_csv('kw_todas.csv')
comp=kw[kw.tipo=='comp']
names={'D1_Bens_Duraveis_MC_PC1':'D1 PC1 (bens duráveis)','D2_Humano_Trabalho_MC_PC1':'D2 PC1 (formalização e escolaridade)','D3_Fisico_Habitacional_Territorial_MC_PC1':'D3 PC1 (infraestrutura pública)','D3_Fisico_Habitacional_Territorial_MC_PC2':'D3 PC2 (condição da moradia)','D4_Tecnologico_Conectividade_MC_PC1':'D4 PC1 (conectividade)'}
m=d.groupby('Cluster_K3')[list(names)].mean(); sd=d[list(names)].std()
t=r"""\begin{table}[htbp]
\centering
\caption{Escore médio dos componentes por grupo e teste de Kruskal-Wallis entre grupos}
\label{tab:escores}
\resizebox{\textwidth}{!}{%
\begin{tabular}{l r r r c c c}
\hline
\textbf{Componente} & \textbf{A} & \textbf{B} & \textbf{C} & $H$ & $p$ & $\varepsilon^2$ \\ \hline
"""
for c,n in names.items():
    r=comp[comp['var']==c].iloc[0]
    t+=f"{n} & "+" & ".join(f2(m.loc[g,c]) for g in range(3))+f" & {f1(r.H)} & {fp(r.p)} & {f2(r.eps2)} \\\\\n"
t+=r"""\hline
\end{tabular}}
\par\footnotesize Fonte: elaboração própria. Média dos escores dos segmentos de cada grupo (grupos A a C). Os escores têm média zero e desvio-padrão de 1,2 a 2,2 entre os 49 segmentos, e o sinal de cada componente é arbitrário (a leitura segue as cargas da Seção~\ref{sec:res_pca}). $H$: estatística de Kruskal-Wallis (4 graus de liberdade); $\varepsilon^2$: tamanho de efeito. Os $p$-valores são descritivos, porque os grupos foram formados com os mesmos componentes (Seção~\ref{sec:metodo_robustez}); todos permanecem abaixo de 0,01 após a correção de Benjamini-Hochberg.
\end{table}
"""
open('tab_escores.tex','w').write(t)
# heatmap em desvios-padrão
z=(m/sd); import os as _os
NOMES=(["A — Urbanização consolidada","B — Infraestrutura precária","C — Urbanização parcial"] if _os.environ.get("GRUPOS_NOMES")=="m2" else ["A — Urbanização consolidada","B — Urbanização parcial","C — Infraestrutura precária"]); _n=list(d.groupby("Cluster_K3").size()) if "d" in dir() else [17,24,8]; z.index=[f"{NOMES[i]} (n={n})" for i,n in zip(range(3),_n)]
fig,ax=plt.subplots(figsize=(8,3.6)); im=ax.imshow(z.values,cmap='RdBu',vmin=-2,vmax=2,aspect='auto')
ax.set_xticks(range(z.shape[1])); ax.set_xticklabels([v.split(' (')[0] for v in names.values()],fontsize=9); ax.set_yticks(range(3)); ax.set_yticklabels(z.index,fontsize=9)
for i in range(3):
    for j in range(z.shape[1]): ax.text(j,i,f"{z.values[i,j]:+.2f}".replace('.',','),ha='center',va='center',fontsize=9,color='white' if abs(z.values[i,j])>1.3 else 'black')
plt.colorbar(im,label='Escore médio (em desvios-padrão)',shrink=.8); ax.set_title('Escores médios dos componentes por grupo',fontsize=10,fontweight='bold')
plt.tight_layout(); plt.savefig('figuras/escores_grupos.png',dpi=200); plt.close()
# ---------- anexo KW 34 variáveis
lab={"log_renda_domiciliar_pc":"Renda domiciliar \\textit{per capita} (log)","possui_automovel":"Automóvel","possui_motocicleta":"Motocicleta","possui_geladeira":"Geladeira","possui_freezer":"\\textit{Freezer}","possui_lavadora":"Máquina de lavar","possui_secadora":"Secadora","possui_tv":"Televisão","possui_ar_condicionado":"Ar-condicionado","possui_microondas":"Micro-ondas","escolaridade_adulto":"Escolaridade (25 anos ou mais)","trabalhou_30d":"Trabalhou nos últimos 30 dias","tem_ctps":"Carteira assinada (CTPS)","contribui_previdencia":"Contribui para a previdência","ivs_df":"IVS-DF","comodos_por_pessoa":"Cômodos por pessoa","banheiros_por_pessoa":"Banheiros por pessoa","agua_rede":"Água da rede","agua_poco_artesiano":"Poço artesiano","agua_poco_cisterna":"Poço cisterna","captacao_chuva":"Captação de chuva","caixa_agua":"Caixa d'água","solucao_sanitaria_estruturada":"Solução sanitária estruturada","acesso_energetico_estruturado":"Acesso energético estruturado","destino_com_coleta_direta":"Coleta direta de lixo","rua_pavimentada":"Rua pavimentada","iluminacao_publica":"Iluminação pública","drenagem":"Drenagem","ponto_onibus":"Ponto de ônibus","qtd_estabelecimentos":"Escala da RA (domicílios expandidos)","internet_domiciliar":"Internet domiciliar","internet_banda_larga_fixa":"Banda larga fixa","internet_movel":"Internet móvel","possui_computador":"Computador"}
v=kw[kw.tipo=='var']
t=r"""% Tabela do Anexo: Kruskal-Wallis das variáveis originais entre os três grupos
{\small
\begin{longtable}{@{}l l r r r r c@{}}
\caption{Kruskal-Wallis das variáveis originais (inclusive as que ficaram fora do modelo) entre os três grupos ($N = 49$ segmentos)}
\label{tab:kw_variaveis}\\
\toprule
Dim. & Variável & $H$ & $p$ & $q$ (FDR) & $\varepsilon^2$ & Dif.\ ($q<0{,}05$) \\
\midrule
\endfirsthead
\multicolumn{7}{c}{\tablename~\thetable{} -- continuação}\\
\toprule
Dim. & Variável & $H$ & $p$ & $q$ (FDR) & $\varepsilon^2$ & Dif.\ ($q<0{,}05$) \\
\midrule
\endhead
\midrule
\multicolumn{7}{r}{\textit{continua na página seguinte}}\\
\endfoot
\bottomrule
\multicolumn{7}{@{}p{\linewidth}@{}}{\footnotesize \textbf{Fonte:} elaboração própria. $H$: estatística de Kruskal-Wallis (2 graus de liberdade) sobre as médias dos 49 segmentos, agrupados pelo $K$-Means ($K=3$); $q$: valor $p$ ajustado pela taxa de falsas descobertas de Benjamini-Hochberg, sobre os testes (variáveis originais e componentes); $\varepsilon^2$: tamanho de efeito. Os valores são descritivos, porque os grupos foram formados com as mesmas variáveis.}\\
\endlastfoot
"""
for _,r in v.iterrows():
    t+=f"{r.dim} & {lab[r['var']]} & {f1(r.H)} & {fp(r.p)} & {fp(r.q)} & {f2(r.eps2)} & {'sim' if r.q<.05 else 'não'} \\\\\n"
t+="\\end{longtable}\n}\n"
open('Tabela_kw_variaveis.tex','w').write(t)
# ---------- apoio
r=pd.read_csv('apoio_prop.csv'); r['grupo']=r.grupo.astype(str).str.replace(r'\.0$','',regex=True)
rows=[('incentivo','Incentivo à produção (já recebeu)'),('Pronaf','\\quad Pronaf (entre os beneficiários)'),('ATER','\\quad ATER (entre os beneficiários)'),('Fomento','\\quad Fomento rural (entre os beneficiários)'),('FCO','\\quad FCO (entre os beneficiários)'),('outorga','Outorga de uso da água'),('beneficio','Benefício de programa social'),('pediu_auxilio','Pediu auxílio governamental (via internet)'),('cnpj','Possui CNPJ (ocupados)'),('pix','Transações financeiras pela internet (PIX)'),('curso_agr','Curso profissional agrícola em andamento')]
def tab(by,order,heads,cap,lab_,note):
    q=r[r.by==by]; t=f"\\begin{{table}}[htbp]\n\\centering\n\\caption{{{cap}}}\n\\label{{{lab_}}}\n\\resizebox{{\\textwidth}}{{!}}{{%\n\\begin{{tabular}}{{l {'r '*(len(order)+1)}}}\n\\hline\n\\textbf{{Indicador (\\%)}} & "+" & ".join(heads)+" & \\textbf{Total} \\\\ \\hline\n"
    for v,l in rows:
        cells=[]
        for g in order+['Total']:
            b=r[(r.by=='ALL')] if g=='Total' else q[q.grupo.astype(str)==str(g)]
            x=b[b['var']==v].iloc[0]; s=f1(x.p)
            if v in('Pronaf','ATER','Fomento','FCO') and x.n<30: s+='$^{a}$'
            cells.append(s)
        t+=l+" & "+" & ".join(cells)+" \\\\\n"
    # bases
    for v,l in [('incentivo','Base: domicílios com resposta válida'),('Pronaf','Base: domicílios beneficiários')]:
        cells=[]
        for g in order+['Total']:
            b=r[(r.by=='ALL')] if g=='Total' else q[q.grupo.astype(str)==str(g)]
            cells.append(f"{int(b[b['var']==v].iloc[0].n):,}".replace(',','.'))
        t+="\\hline\n" if v=='incentivo' else ""
        t+=l+" & "+" & ".join(cells)+" \\\\\n"
    t+="\\hline\n\\end{tabular}}\n\\par\\footnotesize "+note+"\n\\end{table}\n"
    return t
nota="Fonte: elaboração própria a partir da PDAD Rural 2022. Proporções ponderadas pelo peso corrigido, entre as respostas válidas (``sim'' ou ``não''); respostas ``não sabe'' e ``não se aplica'' excluídas. Os programas são categorias excludentes entre os domicílios que receberam incentivo (as demais respostas são ``outro'' e ``não sabe''). $^{a}$ Base de menos de 30 domicílios beneficiários: leitura apenas indicativa."
t1=tab('a01subpopulacao',['Agricultura Familiar','Agricultura Empresarial','Assentamento/Agrovila','Características Urbanas'],['Agr. Familiar','Agr. Empresarial','Assentamento','Moradia e Lazer'],'Apoio à produção, acesso a programas e a mercado, por subpopulação','tab:apoio_sub',nota)
t2=tab('Cluster_K3',[0,1,2],['A','B','C'],'Apoio à produção, acesso a programas e a mercado, por grupo da tipologia','tab:apoio_grupo',nota+""+("" if _os.environ.get("GRUPOS_NOMES")=="m2" else " G2 tem apenas 23 domicílios com resposta válida e 2 beneficiários."))
open('tab_apoio_sub.tex','w').write(t1); open('tab_apoio_grupo.tex','w').write(t2)
print(t2)
# ---------- mapas apoio por RA
a=pd.read_csv('apoio_ra.csv')
g=gpd.read_file('shape_df_extracao/regioes_administrativas.shp')   # extrair regioes_administrativas.zip em shape_df_extracao/
def nrm(s):
    s=unicodedata.normalize('NFKD',str(s)).encode('ascii','ignore').decode().lower(); return re.sub(r'\s+',' ',re.sub(r'[^a-z0-9 ]',' ',s)).strip()
def fix(x):
    try: return x.encode('latin1').decode('utf8')
    except Exception: return x
g['k']=g.ra_nome.map(fix).map(nrm).replace({'scia':'scia estrutural','sol nascente e por do sol':'sol nascente por do sol'})
a['k']=a.ra.map(nrm)
items={'incentivo':'Incentivo à produção (%)','outorga':'Outorga de uso da água (%)','beneficio':'Benefício de programa social (%)','cnpj':'CNPJ entre ocupados (%)'}

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

NUM_RA=numerar_ras(d.ra_norm)
fig,axes=plt.subplots(2,2,figsize=(9.5,8.5)); axes=axes.ravel()
for ax,(v,lab_) in zip(axes,items.items()):
    s=a[(a['var']==v)&(a.n>=30)][['k','p']]   # RAs com ao menos 30 respostas
    gg=g.merge(s,on='k',how='left'); gg['ra_norm']=gg['k']
    gg.plot(ax=ax,color='#f0f0f0',edgecolor='#bbbbbb',linewidth=.4)
    gg[gg.p.notna()].plot(ax=ax,column='p',cmap='viridis',edgecolor='black',linewidth=.5,legend=True,legend_kwds={'shrink':.7})
    ax.set_axis_off(); ax.set_title(lab_,fontsize=12,fontweight='bold')
    rotular_numeros(ax,gg,NUM_RA,col='ra_norm',fs=9,dmin=0.065)
plt.tight_layout(); plt.savefig('figuras/mapa_apoio_ra.png',dpi=200,bbox_inches='tight'); plt.close()
