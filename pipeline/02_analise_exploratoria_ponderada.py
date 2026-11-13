# -*- coding: utf-8 -*-
"""
Análise Exploratória dos Dados (AED) da PDAD Rural 2022 — VERSÃO PONDERADA

Todas as estimativas descritivas (frequências, médias, quantis, momentos, tabelas
cruzadas e perfis dos segmentos) são ponderadas pelo peso amostral corrigido.
O que muda em relação à versão anterior:
  1. peso_dom_corrigido (por estrato, idêntico ao script de cluster; N_h do IPEDF) em vez do peso bruto.
  2. peso dos moradores = peso_dom_corrigido do domicílio onde moram. O peso_mor bruto do arquivo
     perdeu o separador decimal e mistura escalas (valores de 97 a 9,9e14); não é utilizável como está.
  3. Parser "Sim, próprio/compartilhado" corrigido (internet_domiciliar). Banda larga e internet
     móvel: "Não se aplica" (domicílio sem internet) = 0.
  4. Renda por pessoa do domicílio (mesma variável do modelo) é a variável principal de renda.
  5. Segmentos com n < 10 domicílios na amostra são sinalizados; rankings são mostrados com e sem eles.
Declaração para a redação: "As estimativas descritivas são ponderadas pelo peso amostral
corrigido (5.075 domicílios; 15.357 moradores), que expande a amostra para 66.540 domicílios rurais."
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, spearmanr

warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None
pd.set_option("display.float_format", lambda x: "%.2f" % x)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50); pd.set_option("display.max_rows", 200)
sns.set_theme(style="whitegrid", palette="muted")
plt.rcParams["axes.unicode_minus"] = False

try:
    from IPython.display import display  # Colab / Jupyter
except Exception:
    def display(x): print(x.to_string() if hasattr(x, "to_string") else x)

# ---- Configuração -------------------------------------------------------------
PATH_DOM = os.environ.get("PDADR_DOM", "/content/drive/MyDrive/pdadr_domicilios.xlsx")
PATH_MOR = os.environ.get("PDADR_MOR", "/content/drive/MyDrive/pdadr_moradores.xlsx")
FIG_DIR = os.environ.get("PDADR_FIG_DIR", "")          # se definido, salva os PNGs nessa pasta
N_MIN_SEG = 5          # critério de inclusão do segmento (igual ao script de cluster)
N_ALERTA = 10          # segmentos com n abaixo disso são sinalizados como frágeis
if FIG_DIR:
    os.makedirs(FIG_DIR, exist_ok=True)
    matplotlib.use("Agg")

def fim_fig(nome):
    plt.tight_layout()
    if FIG_DIR: plt.savefig(os.path.join(FIG_DIR, nome + ".png"), dpi=200, bbox_inches="tight")
    plt.show(); plt.close()

def titulo(t):
    print("\n" + "=" * 90 + f"\n{t}\n" + "=" * 90)

# ==============================================================================
# BLOCO 1: CARGA E PESOS
# ==============================================================================
df_dom = pd.read_excel(PATH_DOM); df_mor = pd.read_excel(PATH_MOR)
df_dom.columns = df_dom.columns.astype(str).str.strip(); df_mor.columns = df_mor.columns.astype(str).str.strip()
ID, RA, SUB, EST = "domicilio_id", "a01ra", "a01subpopulacao", "a01estratoV2"
NPESS = "npessoas_qtd_de_pessoas_no_domicilio"
for d in (df_dom, df_mor):
    for c in (RA, SUB, EST):
        d[c] = (d[c].astype("string").str.replace(r"\s+", " ", regex=True).str.strip()
                .str.replace("Assentamento/ Agrovila", "Assentamento/Agrovila", regex=False))

N_H = {  # universo N_h publicado pelo IPEDF (21 estratos), total = 66.540 domicílios
    "Agricultura Empresarial - Leste": 5087, "Agricultura Empresarial - Norte": 2968,
    "Agricultura Empresarial - Oeste": 2613, "Agricultura Empresarial - Sul": 4663,
    "Agricultura Familiar - Central Adjacente I": 595, "Agricultura Familiar - Central Adjacente II": 688,
    "Agricultura Familiar - Leste": 1483, "Agricultura Familiar - Norte": 6097,
    "Agricultura Familiar - Oeste": 6998, "Agricultura Familiar - Sul": 1216,
    "Assentamento/Agrovila - Central Adjacente I": 106, "Assentamento/Agrovila - Leste": 1213,
    "Assentamento/Agrovila - Norte": 318, "Assentamento/Agrovila - Oeste": 465,
    "Assentamento/Agrovila - Sul": 351, "Características Urbanas - Central Adjacente I": 218,
    "Características Urbanas - Central Adjacente II": 5452, "Características Urbanas - Leste": 9664,
    "Características Urbanas - Norte": 5556, "Características Urbanas - Oeste": 5780,
    "Características Urbanas - Sul": 5009,
}
titulo("0. PESOS AMOSTRAIS")
df_dom["peso_dom"] = pd.to_numeric(df_dom["peso_dom"], errors="coerce")
obs = df_dom.groupby(EST).agg(n_h=(ID, "count"), w_orig=("peso_dom", "first")).reset_index()
obs["N_h"] = obs[EST].map(N_H)
assert obs["N_h"].notna().all(), "Estrato sem N_h"
obs["w_esperado"] = obs["N_h"] / obs["n_h"]
obs["peso_dom_corrigido"] = obs["w_orig"] / 10 ** np.round(np.log10(obs["w_orig"] / obs["w_esperado"]))
assert np.allclose(obs["n_h"] * obs["peso_dom_corrigido"], obs["N_h"], rtol=0, atol=1e-8), "n_h * w != N_h"
df_dom = df_dom.merge(obs[[EST, "peso_dom_corrigido"]], on=EST, how="left")
W = "peso_dom_corrigido"
df_mor["peso_mor_corrigido"] = df_mor[ID].map(df_dom.set_index(ID)[W])
assert df_mor["peso_mor_corrigido"].notna().all(), "Morador sem domicílio correspondente"
WM = "peso_mor_corrigido"
print(f"Domicílios na amostra: {len(df_dom):,} | Moradores na amostra: {len(df_mor):,}")
print(f"Domicílios expandidos (soma dos pesos): {df_dom[W].sum():,.0f}  (esperado: {sum(N_H.values()):,})")
print(f"Moradores com o peso do domicílio (soma): {df_mor[WM].sum():,.0f}  -> usar apenas para PROPORÇÕES e MÉDIAS, não para totais populacionais")
print(f"Efeito de desenho (Kish) dos domicílios: deff = {len(df_dom) * (df_dom[W]**2).sum() / df_dom[W].sum()**2:.2f}"
      f" | n efetivo = {df_dom[W].sum()**2 / (df_dom[W]**2).sum():,.0f}")
print(f"Efeito de desenho (Kish) dos moradores: deff = {len(df_mor) * (df_mor[WM]**2).sum() / df_mor[WM].sum()**2:.2f}"
      f" | n efetivo = {df_mor[WM].sum()**2 / (df_mor[WM]**2).sum():,.0f}")

# ---- Funções ponderadas -----------------------------------------------------------
def wmean(x, w):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float); w = np.asarray(w, float)
    ok = ~np.isnan(x) & ~np.isnan(w) & (w > 0)
    return np.average(x[ok], weights=w[ok]) if ok.sum() else np.nan

def wquantile(x, w, q):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float); w = np.asarray(w, float)
    ok = ~np.isnan(x) & ~np.isnan(w) & (w > 0); x, w = x[ok], w[ok]
    if len(x) == 0: return np.full(len(np.atleast_1d(q)), np.nan)
    i = np.argsort(x); x, w = x[i], w[i]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(np.atleast_1d(q), c, x)

def wmoments(x, w):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float); w = np.asarray(w, float)
    ok = ~np.isnan(x) & ~np.isnan(w) & (w > 0); x, w = x[ok], w[ok]
    if len(x) < 3: return {}
    m = np.average(x, weights=w); v = np.average((x - m) ** 2, weights=w); s = np.sqrt(v)
    q = wquantile(x, w, [.25, .5, .75])
    return {"N_amostra": len(x), "Média pond.": m, "Q1": q[0], "Mediana": q[1], "Q3": q[2], "DP pond.": s,
            "CV_%": s / m * 100 if m else np.nan,
            "Assimetria": np.average(((x - m) / s) ** 3, weights=w) if s > 0 else np.nan,
            "Curtose (exc.)": np.average(((x - m) / s) ** 4, weights=w) - 3 if s > 0 else np.nan}

def freq_pond(df, col, wcol, ordem=None):
    t = df.groupby(col, observed=False).agg(n_amostra=(col, "size"), peso=(wcol, "sum"))
    t["% amostra"] = t["n_amostra"] / t["n_amostra"].sum() * 100
    t["% ponderado"] = t["peso"] / t["peso"].sum() * 100
    t = t.drop(columns="peso")
    return t.sort_values("% ponderado", ascending=False) if ordem is None else t.reindex(ordem)

def tab_pct_pond(df, var, grupo, wcol):
    t = df.dropna(subset=[var, grupo]).groupby([grupo, var], observed=False)[wcol].sum().unstack(fill_value=0)
    return t.div(t.sum(axis=1), axis=0) * 100

# ==============================================================================
# BLOCO 2: SANEAMENTO E FEATURE ENGINEERING
# ==============================================================================
SENT = [9999, 99999, 999999, 9999999, 99999999, 8888, 88888, 888888, 8888888, 7777, 77777, 777777, 7777777,
        "9999", "99999", "999999", "8888", "88888", "7777", "77777", 888891600, 855555554, 785555554]
def limpar(df):
    d = df.copy()
    SENT_NUM = [x for x in SENT if not isinstance(x, str)]; SENT_STR = [x for x in SENT if isinstance(x, str)]
    for c in d.columns:                       # por coluna: evita falha do DataFrame.replace em lote em algumas versões do pandas
        if pd.api.types.is_numeric_dtype(d[c]): d[c] = d[c].mask(d[c].isin(SENT_NUM))
        elif d[c].dtype == object or str(d[c].dtype) == "string": d[c] = d[c].mask(d[c].astype("string").str.strip().isin(SENT_STR))
    for c in d.select_dtypes(include=[np.number]).columns:
        if "renda" in c.lower(): d.loc[d[c] >= 100000, c] = np.nan
    return d
dom = limpar(df_dom); mor = limpar(df_mor)
dom[W] = df_dom[W]; mor[WM] = df_mor[WM]          # pesos nunca passam pela limpeza de sentinelas

def parse_binaria(s, nao_se_aplica_zero=False):
    """'Sim' e 'Sim, ...' = 1; 'Não' = 0; opcionalmente 'Não se aplica' = 0."""
    s = s.astype("string").str.strip().str.lower()
    r = pd.Series(np.nan, index=s.index, dtype=float)
    r[s.str.startswith("sim", na=False) | s.isin(["1", "1.0"])] = 1.0
    zeros = ["não", "nao", "0", "0.0"] + (["não se aplica", "nao se aplica"] if nao_se_aplica_zero else [])
    r[s.isin(zeros)] = 0.0
    return r

def parse_posse(s):
    t = s.astype("string").str.strip().str.lower()
    r = pd.Series(np.nan, index=s.index, dtype=float)
    r[t.isin(["não tem", "nao tem"])] = 0.0
    n = pd.to_numeric(s, errors="coerce"); r[n.notna()] = (n[n.notna()] > 0).astype(float)
    return r

def combine_max(df, cols, parser=parse_posse):
    return pd.concat([parser(df[c]) for c in cols if c in df.columns], axis=1).max(axis=1, skipna=True)

g = lambda c: dom[c]
dom["renda_domiciliar_pc"] = np.maximum(
    pd.to_numeric(dom["renda_domiciliar_per_capita"], errors="coerce") if "renda_domiciliar_per_capita" in dom.columns
    else pd.to_numeric(dom["renda_domiciliar_calculada"] / dom[NPESS], errors="coerce"), 0.0)
dom["renda_domiciliar_total"] = pd.to_numeric(dom["renda_domiciliar_calculada"], errors="coerce")
dom["comodos_por_pessoa"] = pd.to_numeric(g("B11_quantidade_de_comodos"), errors="coerce") / dom[NPESS]
dom["banheiros_por_pessoa"] = pd.to_numeric(g("B13_quantidade_de_banheiros_e_ou_sanitarios"), errors="coerce") / dom[NPESS]
dom["possui_automovel"] = parse_binaria(g("C01_possui_automovel")); dom["possui_motocicleta"] = parse_binaria(g("C02_possui_motocicleta"))
dom["possui_bicicleta"] = parse_binaria(g("C03_possui_bicicleta"))
dom["possui_geladeira"] = combine_max(dom, ["C07_2_qtd_de_geladeiras_de_uma_porta", "C07_3_qtd_de_geladeiras_de_duas_ou_mais_portas"])
dom["possui_freezer"] = parse_posse(g("C07_4_qtd_de_freezers"))
dom["possui_lavadora"] = combine_max(dom, ["C07_5_qtd_de_maquinas_de_lavar_roupas", "C07_6_qtd_de_maquinas_de_lavar_e_secar_roupas"])
dom["possui_secadora"] = parse_posse(g("C07_7_qtd_de_secadoras_de_roupa"))
dom["possui_tv"] = combine_max(dom, ["C07_9_qtd_de_televisores_tubo", "C07_10_qtd_de_televisores_tela_fina_plana"])
dom["possui_ar_condicionado"] = parse_posse(g("C07_14_qtd_de_ar_condicionado")); dom["possui_microondas"] = parse_posse(g("C07_16_qtd_de_fornos_microondas"))
dom["agua_rede"] = parse_binaria(g("B14_1_abast_de_agua_rede_geral_caesb")); dom["agua_poco_artesiano"] = parse_binaria(g("B14_2_abast_de_agua_poco_artesiano"))
dom["agua_poco_cisterna"] = parse_binaria(g("B14_3_abast_de_agua_poco_cisterna"))
dom["captacao_chuva"] = parse_binaria(g("B14_7_abast_de_agua_neste_domicilio_e_feita_captacao_de_agua_da_chuva"))
dom["caixa_agua"] = parse_binaria(g("B14_6_abast_de_agua_neste_domicilio_ha_caixa_d_agua"))
dom["solucao_sanitaria_estruturada"] = combine_max(dom, ["B15_1_esgot_sanitario_rede_geral_caesb", "B15_2_esgot_sanitario_fossa_septica", "B15_5_esgot_sanitario_bacia_de_evapotranspiracao"], parse_binaria)
dom["acesso_energetico_estruturado"] = combine_max(dom, ["B16_1_energia_eletrica_rede_geral_ceb", "B16_3_energia_eletrica_proprio_energia_solar", "B16_4_energia_eletrica_outras_fontes_renovaveis"], parse_binaria)
dom["destino_com_coleta_direta"] = combine_max(dom, ["B17_1_tratamento_do_lixo_coleta_seletiva_direta", "B17_2_tratamento_do_lixo_coleta_convencional_direta_ou_nao_seletiva"], parse_binaria)
dom["rua_pavimentada"] = parse_binaria(g("B18_1_infraest_urbana_a_rua_de_acesso_principal_asfaltada_pavimentada"))
dom["iluminacao_publica"] = parse_binaria(g("B18_5_infraest_urbana_na_rua_tem_iluminacao")); dom["drenagem"] = parse_binaria(g("B18_6_infraest_urbana_na_rua_tem_drenagem_de_agua_de_chuva"))
dom["ponto_onibus"] = parse_binaria(g("B20_8_nas_proximidades_existe_ponto_de_onibus"))
dom["internet_domiciliar"] = parse_binaria(g("C05_o_domicilio_ou_alguem_do_domicilio_possuia_acesso_a_internet_no_ultimo_mes"))
dom["internet_banda_larga_fixa"] = parse_binaria(g("C05_1_1_acesso_por_internet_banda_larga_fixa"), True)
dom["internet_movel"] = parse_binaria(g("C05_1_2_acesso_por_internet_de_rede_de_celular"), True)
dom["possui_computador"] = combine_max(dom, ["C07_12_qtd_de_microcomputadores_desktop", "C07_13_qtd_de_notebooks_laptops"])
dom["incentivo_producao"] = parse_binaria(g("B22_1_a_sua_propriedade_ja_recebeu_incentivos_a_producao")) if "B22_1_a_sua_propriedade_ja_recebeu_incentivos_a_producao" in dom.columns else np.nan

# Rótulos reais da coluna `escolaridade` (7 categorias). O mapa antigo usava "ensino fundamental incompleto" etc.,
# que não casam com "Fundamental Incompleto"; só "Sem Instrução" e "Superior ..." eram mapeados (~77% de nulos em adultos).
MAPA_ANOS = {"sem instrução": 0, "fundamental incompleto": 4.5, "fundamental completo": 9, "médio incompleto": 10.5,
             "médio completo": 12, "superior incompleto": 14, "superior completo": 16}
mor["idade_num"] = pd.to_numeric(mor["E02_idade__Idade_informada"], errors="coerce")
mor["faixa_etaria"] = pd.cut(mor["idade_num"], bins=[0, 11, 21, 31, 41, 51, 61, 71, 120], right=False,
                             labels=["0-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71+"])
mor["anos_estudo"] = mor["escolaridade"].astype("string").str.strip().str.lower().map(MAPA_ANOS)
mor["escolaridade_adulto"] = np.where(mor["idade_num"] >= 25, mor["anos_estudo"], np.nan)
mor["trabalhou_30d"] = parse_binaria(mor["I05__Trabalhou_nos_últimos_30_dias"]); mor["tem_ctps"] = parse_binaria(mor["I17__Possui_CTPS_assinada"])
mor["contribui_previdencia"] = parse_binaria(mor["I18__Contribui_para_alguma_Previdência_Social_Pública"])
mor["recebe_beneficios_sociais"] = parse_binaria(mor["I24_2__Recebe_ou_já_recebeu_benefícios_de_programas_sociais"])
mor["pedido_auxilio_governamental"] = parse_binaria(mor["F06_7__Pedido_de_auxílio_benefícios_governamentais"])
COL_VB = "I22_1__Valor__R___recebido_no_mês_passado_de_benefícios_sociais__Bolsa_família__BPC_LOAS__Bolsa_de_estudo__valor"
mor["valor_beneficios_sociais"] = pd.to_numeric(mor[COL_VB], errors="coerce") if COL_VB in mor.columns else np.nan   # NaN = sem informação (não vira 0)
mor["internet_uso_3m"] = parse_binaria(mor["F04__Acessou_a_Internet_nos_últimos_3_meses"])
mor["internet_computador"] = parse_binaria(mor["F05_1__Utilizou_microcomputador_para_acessar_a_internet"])
mor["sexo"] = mor["E03_sexo__Sexo_informado"] if "E03_sexo__Sexo_informado" in mor.columns else mor.get("E03_sexo")
mor["renda_individual"] = pd.to_numeric(mor["renda_individual_total"], errors="coerce")
print("\nSaneamento e feature engineering concluídos (parser 'Sim, ...' corrigido).")

# ==============================================================================
# BLOCO 3: DADOS FALTANTES (após saneamento, variáveis derivadas)
# ==============================================================================
titulo("1. DADOS FALTANTES NAS VARIÁVEIS DERIVADAS")
vd_dom = ["renda_domiciliar_pc", "comodos_por_pessoa", "banheiros_por_pessoa", "possui_automovel", "possui_geladeira", "agua_rede",
          "solucao_sanitaria_estruturada", "acesso_energetico_estruturado", "destino_com_coleta_direta", "rua_pavimentada", "iluminacao_publica",
          "drenagem", "ponto_onibus", "internet_domiciliar", "internet_banda_larga_fixa", "internet_movel", "possui_computador", "incentivo_producao"]
vd_mor = ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia", "recebe_beneficios_sociais", "valor_beneficios_sociais",
          "internet_uso_3m", "renda_individual", "idade_num"]
for nome, df, vs in (("Domicílios", dom, vd_dom), ("Moradores", mor, vd_mor)):
    m = (df[vs].isna().mean() * 100).sort_values(ascending=False)
    print(f"\n% de nulos — {nome}:"); print(m.round(1).to_string())
    m = m[m > 0]
    if len(m):
        plt.figure(figsize=(10, 5)); ax = sns.barplot(x=m.values, y=m.index, color="#c0392b", edgecolor="black")
        for i, v in enumerate(m.values): ax.text(v + m.values.max() * 0.01, i, f"{v:.1f}%", va="center", fontsize=9)
        plt.title(f"% de nulos após saneamento — {nome}", fontweight="bold"); plt.xlabel("%"); plt.ylabel("")
        fim_fig(f"Missing_{nome}")

# ==============================================================================
# BLOCO 4: UNIVARIADA (PONDERADA)
# ==============================================================================
titulo("2. DISTRIBUIÇÕES UNIVARIADAS PONDERADAS")
ORD_SUB = None
print("\n--- Subpopulação (domicílios): % na amostra x % ponderado ---")
t_sub = freq_pond(dom, SUB, W); display(t_sub.round(1))
print("\n--- Região Administrativa (domicílios): % na amostra x % ponderado ---")
t_ra = freq_pond(dom, RA, W); display(t_ra.round(1))
print("\n--- Faixa etária (moradores): % na amostra x % ponderado ---")
t_fx = freq_pond(mor.dropna(subset=["faixa_etaria"]), "faixa_etaria", WM, ordem=list(mor["faixa_etaria"].cat.categories)); display(t_fx.round(1))
print(f"\nIdade média ponderada: {wmean(mor['idade_num'], mor[WM]):.1f} | mediana ponderada: {wquantile(mor['idade_num'], mor[WM], .5)[0]:.1f}"
      f" | média simples: {mor['idade_num'].mean():.1f}")
print(f"Tamanho médio do domicílio: ponderado {wmean(df_dom[NPESS], df_dom[W]):.2f} | amostra {df_dom[NPESS].mean():.2f}")

for t, nome in ((t_sub, "Subpopulacao_distribuicao"), (t_ra, "RA_distribuicao")):
    d = t.reset_index().melt(id_vars=t.index.name, value_vars=["% amostra", "% ponderado"], var_name="Base", value_name="%")
    plt.figure(figsize=(10, 5 if nome.startswith("Sub") else 8))
    if nome.startswith("Sub"):
        ax = sns.barplot(data=d, x=t.index.name, y="%", hue="Base", edgecolor="black"); plt.xticks(rotation=20, ha="right")
    else:
        ax = sns.barplot(data=d, y=t.index.name, x="%", hue="Base", edgecolor="black")
    plt.title("Distribuição dos domicílios: amostra x expandida (ponderada)", fontweight="bold"); fim_fig(nome)

# Pirâmide etária ponderada (proporção da população rural)
pir = mor.dropna(subset=["idade_num", "sexo"]).copy()
pir["faixa"] = pd.cut(pir["idade_num"], bins=list(range(0, 110, 10)), right=False, labels=[f"{i}-{i+9}" for i in range(0, 100, 10)])
pt = pir.groupby(["faixa", "sexo"], observed=False)[WM].sum().unstack(fill_value=0)
pt = pt / pt.values.sum() * 100
print("\n--- Pirâmide etária (% da população rural, ponderada) ---"); display(pt.round(2))
if {"Masculino", "Feminino"} <= set(pt.columns):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(pt.index, -pt["Masculino"], color="#2b5c8f", label="Masculino"); ax.barh(pt.index, pt["Feminino"], color="#d95f02", label="Feminino")
    ax.set_xticklabels([f"{abs(x):.0f}%" for x in ax.get_xticks()]); ax.set_xlabel("% da população rural (ponderado)"); ax.set_ylabel("Faixa etária")
    ax.set_title("Pirâmide etária — PDAD Rural 2022 (ponderada)", fontweight="bold"); ax.legend(); fim_fig("Piramide_etaria")
    print(f"Razão de sexos (H/M) ponderada: {pt['Masculino'].sum() / pt['Feminino'].sum() * 100:.1f} homens por 100 mulheres")

def plot_num_pond(df, col, wcol, titulo_, limite=None, nome=None):
    d = df[[col, wcol]].dropna(); d = d[d[wcol] > 0]
    if limite: d = d[d[col] <= limite]
    q = wquantile(d[col], d[wcol], [.05, .25, .5, .75, .95])
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    sns.histplot(x=d[col], weights=d[wcol], bins=30, kde=False, ax=ax[0], color="#2b5c8f"); ax[0].set_title("Histograma ponderado", fontweight="bold")
    ax[1].barh(["P5", "P25", "Mediana", "P75", "P95"], q, color="#7fa9c7", edgecolor="black"); ax[1].set_title("Quantis ponderados", fontweight="bold")
    for i, v in enumerate(q): ax[1].text(v, i, f" {v:,.0f}", va="center", fontsize=9)
    sns.ecdfplot(x=d[col], weights=d[wcol], ax=ax[2], color="#d95f02", lw=2); ax[2].axhline(.5, color="gray", ls="--", alpha=.6)
    ax[2].set_title("ECDF ponderada", fontweight="bold")
    plt.suptitle(titulo_, fontsize=14, fontweight="bold", y=1.03); fim_fig(nome or col)

# ---- Renda ------------------------------------------------------------------------
print("\n--- Renda: média e quantis ponderados (e valores sem peso para comparação) ---")
linhas = []
for nome, df, col, wc in (("Renda domiciliar total (R$)", dom, "renda_domiciliar_total", W), ("Renda domiciliar per capita (R$)", dom, "renda_domiciliar_pc", W),
                          ("Renda individual (R$)", mor, "renda_individual", WM)):
    q = wquantile(df[col], df[wc], [.25, .5, .75, .9])
    linhas.append({"Variável": nome, "N válidos": int(df[col].notna().sum()), "Média pond.": wmean(df[col], df[wc]), "Q1 pond.": q[0], "Mediana pond.": q[1],
                   "Q3 pond.": q[2], "P90 pond.": q[3], "Média amostra": df[col].mean(), "Mediana amostra": df[col].median()})
display(pd.DataFrame(linhas).set_index("Variável").round(1))
for nome, df, col, wc, lim in (("renda_domiciliar_total", dom, "renda_domiciliar_total", W, 2000), ("renda_domiciliar_pc", dom, "renda_domiciliar_pc", W, 1000),
                               ("renda_individual", mor, "renda_individual", WM, 2000)):
    d = df[[col, wc]].dropna()
    print(f"% abaixo de R$ {lim:,}  em {nome}: ponderado {d.loc[d[col] < lim, wc].sum() / d[wc].sum() * 100:.1f}% | amostra {(d[col] < lim).mean() * 100:.1f}%")
plot_num_pond(dom, "renda_domiciliar_total", W, "Renda domiciliar total (ponderada)", 30000, "Renda_domiciliar")
plot_num_pond(dom, "renda_domiciliar_pc", W, "Renda domiciliar per capita (ponderada)", 10000, "Renda_domiciliar_pc")
plot_num_pond(mor, "renda_individual", WM, "Renda individual (ponderada)", 15000, "Renda_individual")

# ==============================================================================
# BLOCO 5: MOMENTOS PONDERADOS
# ==============================================================================
titulo("3. MOMENTOS DESCRITIVOS PONDERADOS")
vm_dom = ["renda_domiciliar_pc", "comodos_por_pessoa", "banheiros_por_pessoa"] + [v for v in vd_dom if v not in ("renda_domiciliar_pc", "comodos_por_pessoa", "banheiros_por_pessoa")]
tm = pd.DataFrame({v: wmoments(dom[v], dom[W]) for v in vm_dom}).T
print("\n--- Domicílios (ponderado) ---"); display(tm.round(3))
vm_mor = ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia", "recebe_beneficios_sociais", "internet_uso_3m", "renda_individual", "idade_num"]
tmm = pd.DataFrame({v: wmoments(mor[v], mor[WM]) for v in vm_mor}).T
print("\n--- Moradores (ponderado) ---"); display(tmm.round(3))
print("\n--- Média ponderada por subpopulação (domicílios) ---")
display(pd.DataFrame({s: {v: wmean(d[v], d[W]) for v in ["renda_domiciliar_pc", "agua_rede", "solucao_sanitaria_estruturada", "rua_pavimentada",
                                                       "internet_domiciliar", "internet_banda_larga_fixa", "possui_computador"]}
                      for s, d in dom.groupby(SUB)}).T.round(3))

# ==============================================================================
# BLOCO 6: BIVARIADAS PONDERADAS
# ==============================================================================
titulo("4. ANÁLISES BIVARIADAS PONDERADAS")
def renda_por(df, col, wcol, nome_fig, tit):
    d = df.dropna(subset=[col, "renda_individual"])
    t = d.groupby(col, observed=False).apply(lambda x: pd.Series({"n": len(x), "Média pond. (R$)": wmean(x["renda_individual"], x[wcol]),
                                                                  "Mediana pond. (R$)": wquantile(x["renda_individual"], x[wcol], .5)[0]}))
    print(f"\n{tit}"); display(t.round(1))
    m = t.reset_index().melt(id_vars=col, value_vars=["Média pond. (R$)", "Mediana pond. (R$)"], var_name="Estatística", value_name="R$")
    plt.figure(figsize=(10, 5)); ax = sns.barplot(data=m, x=col, y="R$", hue="Estatística", palette="Blues_r", edgecolor="black")
    plt.xticks(rotation=25 if len(t) > 4 else 0, ha="right"); plt.title(tit, fontweight="bold")
    for p in ax.patches:
        h = p.get_height()
        if h > 0: ax.annotate(f"{h:.0f}", (p.get_x() + p.get_width() / 2, h), ha="center", va="bottom", fontsize=8)
    fim_fig(nome_fig)
ORD_ESC = [k for k in MAPA_ANOS]
mor["escolaridade_cat"] = pd.Categorical(mor["escolaridade"].astype("string").str.strip().str.lower(), categories=ORD_ESC, ordered=True)
renda_por(mor, "escolaridade_cat", WM, "Renda_escolaridade", "Renda individual por escolaridade (ponderada)")
renda_por(mor, "sexo", WM, "Renda_sexo", "Renda individual por sexo (ponderada)")
renda_por(mor, "faixa_etaria", WM, "Renda_faixa_etaria", "Renda individual por faixa etária (ponderada)")
renda_por(mor, "recebe_beneficios_sociais", WM, "Renda_beneficio", "Renda individual por recebimento de benefício (ponderada)")

for v, nome in (("internet_uso_3m", "Usou internet nos últimos 3 meses"), ("recebe_beneficios_sociais", "Recebe/recebeu benefícios sociais"), ("tem_ctps", "CTPS assinada")):
    print(f"\n--- {nome} por subpopulação (% ponderado) ---"); display(tab_pct_pond(mor, v, SUB, WM).round(1))

def cramers_v_pond(df, x, y, wcol):
    d = df.dropna(subset=[x, y, wcol]); t = d.groupby([x, y], observed=False)[wcol].sum().unstack(fill_value=0)
    t = t.loc[t.sum(axis=1) > 0, t.sum(axis=0) > 0]
    if min(t.shape) < 2: return np.nan
    n_ef = d[wcol].sum() ** 2 / (d[wcol] ** 2).sum()          # n efetivo (Kish): evita inflar o qui-quadrado com pesos grandes
    t = t / t.values.sum() * n_ef
    chi2 = chi2_contingency(t, correction=False)[0]; r, k = t.shape
    return np.sqrt(max(0, chi2 / n_ef - (k - 1) * (r - 1) / (n_ef - 1)) / min(k - 1, r - 1))
vv = ["agua_rede", "solucao_sanitaria_estruturada", "rua_pavimentada", "internet_domiciliar", "possui_computador"]
print("\n--- V de Cramér ponderado (n efetivo de Kish): subpopulação x variável ---")
print({v: round(cramers_v_pond(dom, SUB, v, W), 3) for v in vv})
print("--- V de Cramér ponderado: RA x variável ---")
print({v: round(cramers_v_pond(dom, RA, v, W), 3) for v in vv})

# ==============================================================================
# BLOCO 7: MATRIZ DOS SEGMENTOS (RA x SUBPOPULAÇÃO), TODOS PONDERADOS
# ==============================================================================
titulo("5. MATRIZ DOS SEGMENTOS TERRITORIAIS (RA x SUBPOPULAÇÃO)")
CH = [RA, SUB]
VS_DOM = ["renda_domiciliar_pc", "renda_domiciliar_total", "comodos_por_pessoa", "banheiros_por_pessoa", "possui_automovel", "possui_motocicleta", "possui_bicicleta",
          "possui_geladeira", "possui_freezer", "possui_lavadora", "possui_secadora", "possui_tv", "possui_ar_condicionado", "possui_microondas", "agua_rede",
          "agua_poco_artesiano", "agua_poco_cisterna", "captacao_chuva", "caixa_agua", "solucao_sanitaria_estruturada", "acesso_energetico_estruturado",
          "destino_com_coleta_direta", "rua_pavimentada", "iluminacao_publica", "drenagem", "ponto_onibus", "internet_domiciliar", "internet_banda_larga_fixa",
          "internet_movel", "possui_computador", "incentivo_producao"]
VS_MOR = ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia", "internet_uso_3m", "internet_computador", "recebe_beneficios_sociais",
          "valor_beneficios_sociais", "pedido_auxilio_governamental", "renda_individual"]
def agrega(df, vs, wcol):
    out = {}
    for k, gp in df.groupby(CH, observed=True):
        out[k] = {v: wmean(gp[v], gp[wcol]) for v in vs if v in gp.columns}
    t = pd.DataFrame(out).T; t.index = pd.MultiIndex.from_tuples(t.index, names=CH); return t
seg_d = agrega(dom, VS_DOM, W); seg_m = agrega(mor, VS_MOR, WM)
cnt = dom.groupby(CH).agg(amostra_n_dom=(ID, "count"), dom_expandidos=(W, "sum"))
cnt_m = mor.groupby(CH).agg(amostra_n_mor=(ID, "count"))
seg = cnt.join(cnt_m).join(seg_d).join(seg_m).reset_index()
seg = seg[seg["amostra_n_dom"] >= N_MIN_SEG].reset_index(drop=True)
seg["frágil_n<10"] = seg["amostra_n_dom"] < N_ALERTA
print(f"Segmentos com n >= {N_MIN_SEG}: {len(seg)} | com n < {N_ALERTA} (frágeis): {int(seg['frágil_n<10'].sum())} | com n < 20: {int((seg['amostra_n_dom'] < 20).sum())}")
print(f"Domicílios expandidos cobertos pelos {len(seg)} segmentos: {seg['dom_expandidos'].sum():,.0f} de {df_dom[W].sum():,.0f}")
print(f"n por segmento: mínimo {seg['amostra_n_dom'].min()} | mediana {seg['amostra_n_dom'].median():.0f} | máximo {seg['amostra_n_dom'].max()}")
print("\nSegmentos frágeis (n < 10):"); display(seg.loc[seg["frágil_n<10"], CH + ["amostra_n_dom", "dom_expandidos"]])

# ==============================================================================
# BLOCO 8: ANÁLISE MACRO (N = 49)
# ==============================================================================
titulo(f"6. ANÁLISE ENTRE SEGMENTOS (N = {len(seg)})")
cols_desc = ["renda_domiciliar_pc", "escolaridade_adulto", "agua_rede", "solucao_sanitaria_estruturada", "rua_pavimentada", "internet_banda_larga_fixa",
             "possui_computador", "recebe_beneficios_sociais", "tem_ctps"]
print("\n--- Estatísticas descritivas entre segmentos (valores já ponderados dentro de cada segmento) ---")
display(seg[cols_desc].describe().T[["mean", "std", "min", "50%", "max"]])

heat = ["renda_domiciliar_pc", "escolaridade_adulto", "internet_banda_larga_fixa", "possui_computador", "agua_rede", "solucao_sanitaria_estruturada",
        "rua_pavimentada", "tem_ctps", "recebe_beneficios_sociais", "valor_beneficios_sociais"]
corr = seg[heat].corr(method="spearman"); print("\n--- Spearman entre segmentos ---"); display(corr.round(2))
plt.figure(figsize=(11, 8)); ax = sns.heatmap(corr, cmap="coolwarm", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right"); plt.title(f"Correlação de Spearman ({len(seg)} segmentos)", fontsize=14); fim_fig("Heatmap_spearman")
for a, b in (("escolaridade_adulto", "renda_domiciliar_pc"), ("agua_rede", "renda_domiciliar_pc"), ("internet_banda_larga_fixa", "renda_domiciliar_pc"),
             ("recebe_beneficios_sociais", "renda_domiciliar_pc")):
    r, p = spearmanr(seg[a], seg[b], nan_policy="omit"); print(f"Spearman {a} x {b}: rho = {r:.2f} (p = {p:.3f})")
seg_ok = seg[~seg["frágil_n<10"]]
r, p = spearmanr(seg_ok["escolaridade_adulto"], seg_ok["renda_domiciliar_pc"], nan_policy="omit")
print(f"Robustez (sem segmentos com n<10, N={len(seg_ok)}): escolaridade x renda pc rho = {r:.2f} (p = {p:.3f})")

SHOW = CH + ["amostra_n_dom", "dom_expandidos", "renda_domiciliar_pc", "escolaridade_adulto", "agua_rede", "solucao_sanitaria_estruturada",
             "internet_banda_larga_fixa", "recebe_beneficios_sociais"]
print("\n--- Top 10 por renda domiciliar per capita (todos os segmentos) ---"); display(seg.sort_values("renda_domiciliar_pc", ascending=False)[SHOW].head(10))
print("\n--- Top 10 por renda domiciliar per capita (apenas n >= 10) ---"); display(seg_ok.sort_values("renda_domiciliar_pc", ascending=False)[SHOW].head(10))
print("\n--- Bottom 10 por renda domiciliar per capita (apenas n >= 10) ---"); display(seg_ok.sort_values("renda_domiciliar_pc")[SHOW].head(10))
SB = CH + ["amostra_n_dom", "recebe_beneficios_sociais", "valor_beneficios_sociais", "tem_ctps", "renda_domiciliar_pc"]
print("\n--- Top 10 por % de moradores que recebem benefícios sociais (apenas n >= 10) ---")
display(seg_ok.sort_values("recebe_beneficios_sociais", ascending=False)[SB].head(10))

titulo("7. PERFIL COMPLETO DOS SEGMENTOS")
SV = CH + ["amostra_n_dom", "dom_expandidos", "renda_domiciliar_pc", "escolaridade_adulto", "agua_rede", "solucao_sanitaria_estruturada", "rua_pavimentada",
           "internet_banda_larga_fixa", "possui_computador", "recebe_beneficios_sociais", "tem_ctps", "frágil_n<10"]
with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 250):
    display(seg.sort_values(CH)[SV].round(3).reset_index(drop=True))

# ---- Prova real: renda individual de quem não recebe benefício ----------------------
titulo("8. PROVA REAL: RENDA INDIVIDUAL DE QUEM NÃO RECEBE BENEFÍCIO SOCIAL (ponderada)")
d = mor[mor["recebe_beneficios_sociais"] == 0.0].dropna(subset=["renda_individual"])
zero = d.loc[d["renda_individual"] == 0, WM].sum() / d[WM].sum() * 100
print(f"Válidos: {len(d)} | média pond.: R$ {wmean(d['renda_individual'], d[WM]):.2f} | mediana pond.: R$ {wquantile(d['renda_individual'], d[WM], .5)[0]:.2f}"
      f" | renda exatamente zero: {zero:.1f}% (ponderado), {(d['renda_individual'] == 0).mean() * 100:.1f}% (amostra)")

seg.to_csv(os.path.join(FIG_DIR, "segmentos_ponderados.csv"), index=False) if FIG_DIR else None
print("\nAED ponderada concluída. Declarar na redação: estimativas ponderadas pelo peso amostral corrigido.")
