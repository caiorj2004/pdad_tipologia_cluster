# -*- coding: utf-8 -*-
"""Análise de Cluster do PDAD Rural — Modelo de Referência 4D (VERSÃO CORRIGIDA)

Mudanças em relação à versão anterior (marcadas no código com "# [CORREÇÃO ...]"):
  1. internet_domiciliar: a resposta "Sim, próprio"/"Sim, compartilhado..." não era reconhecida pelo parse
     (a variável ficava constante em 0 e era descartada do PCA). Agora "Sim..." = 1.
  2. ivs_df: passa a vir da base oficial do IPEDF (IVS-DF 2018-2021), não de valores digitados à mão
     nem do padrão 0,34 para as demais RAs.
  3. banda larga fixa / internet móvel: "Não se aplica" (domicílio sem internet) passa a contar como 0,
     para a variável deixar de ser calculada apenas entre os domicílios que têm internet.
  4. qtd_estabelecimentos mantida; pct_estabelecimentos removida (é qtd / total, correlação perfeita:
     duplicava a variável e inflava um componente artificial).
  5. Escolha do K: benchmark ampliado (silhueta, CH, DB, Gap, BIC), estabilidade por subamostragem
     e dois nulos (embaralhamento por coluna e normal multivariada com a mesma covariância).
  7. escolaridade_adulto: o mapa de rótulos não casava com a coluna `escolaridade` ("Fundamental Incompleto",
     "Médio Completo"...), então só "Sem Instrução" e "Superior ..." eram convertidos. Mapa corrigido.
  8. peso dos moradores: o peso_mor do arquivo perdeu o separador decimal e mistura escalas. Os moradores passam
     a usar o peso corrigido do domicílio onde vivem (mesma probabilidade de inclusão).
  6. LISA: removido o fillna(0) (mascarava variáveis degeneradas); agora há asserção de ausência de NaN.
"""

# Colab: !pip install libpysal esda statsmodels geopandas
import os, re, glob, zipfile, unicodedata, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import (calinski_harabasz_score, davies_bouldin_score, silhouette_score,
                             adjusted_rand_score)
from statsmodels.stats.multitest import multipletests
import libpysal, esda
from libpysal.weights import KNN
warnings.filterwarnings("ignore")

# ==============================================================================
# CONFIGURAÇÃO
# ==============================================================================
BASE = Path(os.environ.get("TCC_BASE", "/content/drive/MyDrive"))
PATH_DOM = BASE / "pdadr_domicilios.xlsx"
PATH_MOR = BASE / "pdadr_moradores.xlsx"
PATH_IVS = BASE / "Base_dados-IVS.xlsx"          # base do IPEDF (IVS-DF 2018-2021)
PATH_SHP_ZIP = Path(os.environ.get("TCC_SHP", BASE / "regioes_administrativas.zip"))
IVS_ANO = 2021                                    # ano mais próximo da PDAD Rural 2022
NAO_SE_APLICA_ZERO = True                         # [CORREÇÃO 3] False reproduz o cálculo condicional antigo
USAR_PCT_ESTAB = False                            # [CORREÇÃO 4] True recoloca a variável redundante
# [V2] Revisão pós-avaliação: cada ajuste pode ser ligado/desligado (V2="" reproduz a versão anterior)
#   ctx   = IVS-DF e escala da RA saem da matriz de agrupamento (ficam só para descrição)
#   agua  = água em escada de serviço (rede / poço / precária) em razões log-isométricas (ilr)
#   esgoto= esgoto em escada de serviço (rede / fossa séptica ou evapotranspiração / rudimentar ou céu aberto) em ilr
#   renda = média do ln(1+renda) domiciliar, com winsorização no P99
#   dur   = geladeira e secadora saem do agrupamento; D1 passa a ser "bens duráveis"
V2 = set(x for x in os.environ.get("V2", "ctx,agua,esgoto,renda,dur").split(",") if x)
K_FINAL = int(os.environ.get("K_FINAL", "3"))   # K=3: maior silhueta entre K>=3, estabilidade 0,66, menor grupo com 8 segmentos (ver Etapas 5 e 5.2)
WINSOR_P = float(os.environ.get("WINSOR_P", "0.99"))   # percentil de winsorização da renda (1.0 = sem winsorização)
COMP = os.environ.get("COMP", "ilr")                    # transformação composicional de água/esgoto: "ilr" (2 coordenadas) ou "clr" (3 colunas)
K_CANDIDATOS = range(2, 9)
N_PERM_NULO = 500                                 # permutações por K nos testes nulos
N_SUBAMOSTRAS = 200                               # subamostras (80%) por K na estabilidade
SEED = 42
RNG = np.random.default_rng(SEED)
CL = f"Cluster_K{K_FINAL}"
NOMES_GRUPOS = {0: "A — Urbanização formal consolidada", 1: "B — Urbanização parcial", 2: "C — Saneamento e vias precários"} if K_FINAL == 3 else {i: f"Grupo {i}" for i in range(K_FINAL)}   # nomes substantivos (ordem: da infraestrutura mais à menos consolidada)
FIG_DIR = Path(os.environ.get("TCC_FIG", "figuras")); FIG_DIR.mkdir(parents=True, exist_ok=True)   # PNGs e CSVs de saída

def salvar(nome):
    plt.tight_layout(); plt.savefig(FIG_DIR / f"{nome}.png", dpi=200, bbox_inches="tight"); plt.show(); plt.close()

# ==============================================================================
# PADRONIZAÇÃO DE NOMES DE RA
# ==============================================================================
DE_PARA_RAS = {
    "sol nascente": "sol nascente por do sol", "por do sol": "sol nascente por do sol",
    "sol nascente por do sol": "sol nascente por do sol", "sol nascente e por do sol": "sol nascente por do sol",
    "nucleo bandeirante": "nucleo bandeirante", "bandeirante": "nucleo bandeirante",
    "n bandeirante": "nucleo bandeirante", "ra viii nucleo bandeirante": "nucleo bandeirante",
    "scia": "scia estrutural", "estrutural": "scia estrutural", "scia estrutural": "scia estrutural",
    "riacho fundo i": "riacho fundo", "riacho fundo": "riacho fundo", "riacho fundo ii": "riacho fundo ii",
    "varjao": "varjao", "fercal": "fercal", "itapoan": "itapoa", "itapoa": "itapoa",
    "sra": "sia", "sia": "sia", "brasilia": "plano piloto", "plano piloto": "plano piloto",
    "parkway": "park way",  # [CORREÇÃO 2] grafia da base do IVS
}

def limpar_nome_ra_base(texto):
    if pd.isna(texto): return ""
    s = unicodedata.normalize("NFKD", str(texto)).encode("ASCII", "ignore").decode("utf-8").lower()
    s = re.sub(r"^ra\s+[ivxlcdm0-9]+\s*[-_–—]?\s*", "", s)
    s = re.sub(r"^regiao\s+administrativa\s+(de\s+)?", "", s)
    s = re.sub(r"[/_\-\.]", " ", s)
    s = re.sub(r"\be\b", " ", s)
    s = " ".join(s.split())
    return DE_PARA_RAS.get(s, s)

# ==============================================================================
# FASE 1: SANEAMENTO, PESOS, AGREGAÇÃO
# ==============================================================================
print("Fase 1: Saneamento, pesos e agregação territorial...\n")
df_dom = pd.read_excel(PATH_DOM)
df_mor = pd.read_excel(PATH_MOR)
df_dom.columns = df_dom.columns.astype(str).str.strip()
df_mor.columns = df_mor.columns.astype(str).str.strip()

ID, RA, SUB, ESTRATO = "domicilio_id", "a01ra", "a01subpopulacao", "a01estratoV2"
PESO_DOM, PESO_MOR = "peso_dom", "peso_mor"
NPESS, ESC = "npessoas_qtd_de_pessoas_no_domicilio", "escolaridade"
CHAVES = [RA, SUB]
print(f"Domicílios: {len(df_dom):,} | Moradores: {len(df_mor):,}")

for df in (df_dom, df_mor):
    for col in (RA, SUB, ESTRATO):
        df[col] = (df[col].astype("string").str.replace(r"\s+", " ", regex=True).str.strip()
                   .str.replace("Assentamento/ Agrovila", "Assentamento/Agrovila", regex=False))
df_dom[PESO_DOM] = pd.to_numeric(df_dom[PESO_DOM], errors="coerce")
df_mor[PESO_MOR] = pd.to_numeric(df_mor[PESO_MOR], errors="coerce")

N_H = {  # universo N_h publicado pelo IPEDF (21 estratos)
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
obs = df_dom.groupby(ESTRATO).agg(n_h=(ID, "count"), w_orig=(PESO_DOM, "first")).reset_index()
obs["N_h"] = obs[ESTRATO].map(N_H)
obs["w_esperado"] = obs["N_h"] / obs["n_h"]
obs["fator_escala"] = 10 ** np.round(np.log10(obs["w_orig"] / obs["w_esperado"]))
obs["peso_dom_corrigido"] = obs["w_orig"] / obs["fator_escala"]
assert np.allclose(obs["n_h"] * obs["peso_dom_corrigido"], obs["N_h"], rtol=0, atol=1e-8), "Erro na validação de N_h!"
print(f"OK Pesos corrigidos para {len(obs)} estratos (N total = {obs['N_h'].sum():,}).")
df_dom = df_dom.merge(obs[[ESTRATO, "peso_dom_corrigido"]], on=ESTRATO, how="left")
# [CORREÇÃO 8] peso do morador = peso corrigido do seu domicílio (peso_mor bruto está corrompido)
df_mor["peso_mor_corrigido"] = df_mor[ID].map(df_dom.set_index(ID)["peso_dom_corrigido"])
assert df_mor["peso_mor_corrigido"].notna().all()
PESO_MOR = "peso_mor_corrigido"

# ---- Funções de parsing ------------------------------------------------------
def parse_binaria(series):
    s = series.astype("string").str.strip().str.lower()
    res = pd.Series(np.nan, index=series.index, dtype=float)
    res[s.isin(["sim", "1"])] = 1.0
    res[s.isin(["não", "nao", "0"])] = 0.0
    return res

def parse_binaria_sim_prefixo(series, nao_se_aplica_zero=False):
    """[CORREÇÃO 1/3] 'Sim, ...' = 1; 'Não' = 0; opcionalmente 'Não se aplica' = 0 (ex.: sem internet => sem banda larga)."""
    s = series.astype("string").str.strip().str.lower()
    res = pd.Series(np.nan, index=series.index, dtype=float)
    res[s.str.startswith("sim", na=False)] = 1.0
    zeros = ["não", "nao", "0"] + (["não se aplica", "nao se aplica"] if nao_se_aplica_zero else [])
    res[s.isin(zeros)] = 0.0
    return res

def parse_posse(series):
    s = series.astype("string").str.strip().str.lower()
    res = pd.Series(np.nan, index=series.index, dtype=float)
    res[s.isin(["não tem", "nao tem"])] = 0.0
    num = pd.to_numeric(series, errors="coerce")
    res[num.notna()] = (num[num.notna()] > 0).astype(float)
    return res

def combine_max(df, cols, parser=parse_posse):
    parsed = pd.concat([parser(df[c]) for c in cols if c in df.columns], axis=1)
    return parsed.max(axis=1, skipna=True)

def media_ponderada_vetorizada(df, group_cols, var_cols, weight_col):
    w = df[weight_col].to_numpy(dtype=float)
    result_df = pd.DataFrame(index=df.groupby(group_cols, dropna=False).size().index).reset_index()
    for var in var_cols:
        x = df[var].to_numpy(dtype=float)
        valid = ~np.isnan(x) & ~np.isnan(w) & (w >= 0)
        tmp = pd.DataFrame({"xw": np.where(valid, x * w, 0.0), "w_val": np.where(valid, w, 0.0)})
        g = tmp.groupby([df[c] for c in group_cols], dropna=False).sum()
        result_df[var] = np.where(g["w_val"] > 0, g["xw"] / g["w_val"], np.nan)
    return result_df

# ---- IVS-DF [CORREÇÃO 2] -----------------------------------------------------
df_dom["ra_norm"] = df_dom[RA].map(limpar_nome_ra_base)
_ivs = pd.read_excel(PATH_IVS)
_ivs = _ivs[_ivs["Ano"] == IVS_ANO].copy()
_ivs["ra_norm"] = _ivs["Território"].map(limpar_nome_ra_base)
MAPA_IVS_DF = dict(zip(_ivs["ra_norm"], _ivs["IVS Geral"]))
faltam_ivs = sorted(set(df_dom["ra_norm"]) - set(MAPA_IVS_DF))
assert not faltam_ivs, f"RAs da PDAD sem IVS na base do IPEDF: {faltam_ivs}"
df_dom["ivs_df"] = df_dom["ra_norm"].map(MAPA_IVS_DF)
print(f"OK IVS-DF {IVS_ANO} associado às {df_dom['ra_norm'].nunique()} RAs da amostra (fonte: base do IPEDF).")

# ---- Escala territorial (qtd_estabelecimentos) -------------------------------
# ATENÇÃO (redação): esta variável é o total de domicílios rurais expandidos da RA na PDAD Rural,
# usado como proxy de escala/adensamento do território. Não é um censo de estabelecimentos.
volumetria_ra = df_dom.groupby("ra_norm")["peso_dom_corrigido"].sum()
df_dom["qtd_estabelecimentos"] = df_dom["ra_norm"].map(volumetria_ra)
df_dom["pct_estabelecimentos"] = df_dom["ra_norm"].map(volumetria_ra / volumetria_ra.sum())

# ---- Atributos domiciliares --------------------------------------------------
df_dom["renda_domiciliar_pc"] = np.maximum(
    pd.to_numeric(df_dom["renda_domiciliar_per_capita"], errors="coerce") if "renda_domiciliar_per_capita" in df_dom.columns
    else pd.to_numeric(df_dom["renda_domiciliar_calculada"] / df_dom[NPESS], errors="coerce"), 0.0)

df_dom["possui_automovel"] = parse_binaria(df_dom["C01_possui_automovel"])
df_dom["possui_motocicleta"] = parse_binaria(df_dom["C02_possui_motocicleta"])
df_dom["possui_bicicleta"] = parse_binaria(df_dom["C03_possui_bicicleta"])
df_dom["possui_geladeira"] = combine_max(df_dom, ["C07_2_qtd_de_geladeiras_de_uma_porta", "C07_3_qtd_de_geladeiras_de_duas_ou_mais_portas"])
df_dom["possui_freezer"] = parse_posse(df_dom["C07_4_qtd_de_freezers"])
df_dom["possui_lavadora"] = combine_max(df_dom, ["C07_5_qtd_de_maquinas_de_lavar_roupas", "C07_6_qtd_de_maquinas_de_lavar_e_secar_roupas"])
df_dom["possui_secadora"] = parse_posse(df_dom["C07_7_qtd_de_secadoras_de_roupa"])
df_dom["possui_tv"] = combine_max(df_dom, ["C07_9_qtd_de_televisores_tubo", "C07_10_qtd_de_televisores_tela_fina_plana"])
df_dom["possui_ar_condicionado"] = parse_posse(df_dom["C07_14_qtd_de_ar_condicionado"])
df_dom["possui_microondas"] = parse_posse(df_dom["C07_16_qtd_de_fornos_microondas"])

df_dom["comodos_por_pessoa"] = pd.to_numeric(df_dom["B11_quantidade_de_comodos"], errors="coerce") / df_dom[NPESS]
df_dom["banheiros_por_pessoa"] = pd.to_numeric(df_dom["B13_quantidade_de_banheiros_e_ou_sanitarios"], errors="coerce") / df_dom[NPESS]

df_dom["agua_rede"] = parse_binaria(df_dom["B14_1_abast_de_agua_rede_geral_caesb"])
df_dom["agua_poco_artesiano"] = parse_binaria(df_dom["B14_2_abast_de_agua_poco_artesiano"])
df_dom["agua_poco_cisterna"] = parse_binaria(df_dom["B14_3_abast_de_agua_poco_cisterna"])
df_dom["captacao_chuva"] = parse_binaria(df_dom["B14_7_abast_de_agua_neste_domicilio_e_feita_captacao_de_agua_da_chuva"])
df_dom["caixa_agua"] = parse_binaria(df_dom["B14_6_abast_de_agua_neste_domicilio_ha_caixa_d_agua"])
df_dom["solucao_sanitaria_estruturada"] = combine_max(df_dom, ["B15_1_esgot_sanitario_rede_geral_caesb", "B15_2_esgot_sanitario_fossa_septica", "B15_5_esgot_sanitario_bacia_de_evapotranspiracao"], parse_binaria)
# [V2] escadas de serviço (hierarquia: vale a melhor fonte declarada; B14 e B15 aceitam várias respostas)
def _sim(col): return parse_binaria(df_dom[col]).fillna(0.0).astype(bool)
_rede_a, _poco_a = _sim("B14_1_abast_de_agua_rede_geral_caesb"), _sim("B14_2_abast_de_agua_poco_artesiano") | _sim("B14_3_abast_de_agua_poco_cisterna")
_nivel_agua = np.where(_rede_a, "rede", np.where(_poco_a, "ind", "prec"))
_ok_agua = parse_binaria(df_dom["B14_1_abast_de_agua_rede_geral_caesb"]).notna() | parse_binaria(df_dom["B14_2_abast_de_agua_poco_artesiano"]).notna()
for k in ("rede", "ind", "prec"): df_dom[f"agua_nivel_{k}"] = np.where(_ok_agua, (_nivel_agua == k).astype(float), np.nan)
_rede_e = _sim("B15_1_esgot_sanitario_rede_geral_caesb")
_sept_e = _sim("B15_2_esgot_sanitario_fossa_septica") | _sim("B15_5_esgot_sanitario_bacia_de_evapotranspiracao")
_nivel_esg = np.where(_rede_e, "rede", np.where(_sept_e, "ind", "prec"))
_ok_esg = parse_binaria(df_dom["B15_1_esgot_sanitario_rede_geral_caesb"]).notna() | parse_binaria(df_dom["B15_2_esgot_sanitario_fossa_septica"]).notna()
for k in ("rede", "ind", "prec"): df_dom[f"esgoto_nivel_{k}"] = np.where(_ok_esg, (_nivel_esg == k).astype(float), np.nan)
# [V2] renda: ln(1+renda pc) por domicílio, winsorizada no P99 (global)
_p99 = df_dom["renda_domiciliar_pc"].quantile(WINSOR_P)
df_dom["renda_pc_winsor"] = df_dom["renda_domiciliar_pc"].clip(upper=_p99)
df_dom["log_renda_hh"] = np.log1p(df_dom["renda_pc_winsor"])
df_dom["acesso_energetico_estruturado"] = combine_max(df_dom, ["B16_1_energia_eletrica_rede_geral_ceb", "B16_3_energia_eletrica_proprio_energia_solar", "B16_4_energia_eletrica_outras_fontes_renovaveis"], parse_binaria)
df_dom["destino_com_coleta_direta"] = combine_max(df_dom, ["B17_1_tratamento_do_lixo_coleta_seletiva_direta", "B17_2_tratamento_do_lixo_coleta_convencional_direta_ou_nao_seletiva"], parse_binaria)
df_dom["rua_pavimentada"] = parse_binaria(df_dom["B18_1_infraest_urbana_a_rua_de_acesso_principal_asfaltada_pavimentada"])
df_dom["iluminacao_publica"] = parse_binaria(df_dom["B18_5_infraest_urbana_na_rua_tem_iluminacao"])
df_dom["drenagem"] = parse_binaria(df_dom["B18_6_infraest_urbana_na_rua_tem_drenagem_de_agua_de_chuva"])
df_dom["ponto_onibus"] = parse_binaria(df_dom["B20_8_nas_proximidades_existe_ponto_de_onibus"])

# [CORREÇÃO 1] e [CORREÇÃO 3]
df_dom["internet_domiciliar"] = parse_binaria_sim_prefixo(df_dom["C05_o_domicilio_ou_alguem_do_domicilio_possuia_acesso_a_internet_no_ultimo_mes"])
df_dom["internet_banda_larga_fixa"] = parse_binaria_sim_prefixo(df_dom["C05_1_1_acesso_por_internet_banda_larga_fixa"], NAO_SE_APLICA_ZERO)
df_dom["internet_movel"] = parse_binaria_sim_prefixo(df_dom["C05_1_2_acesso_por_internet_de_rede_de_celular"], NAO_SE_APLICA_ZERO)
df_dom["possui_computador"] = combine_max(df_dom, ["C07_12_qtd_de_microcomputadores_desktop", "C07_13_qtd_de_notebooks_laptops"])

# ---- Atributos de moradores --------------------------------------------------
# [CORREÇÃO 7] rótulos reais da coluna `escolaridade`; anos de estudo equivalentes por categoria
MAPA_ANOS_ESTUDO = {
    "sem instrução": 0, "fundamental incompleto": 4.5, "fundamental completo": 9, "médio incompleto": 10.5,
    "médio completo": 12, "superior incompleto": 14, "superior completo": 16,
}
col_idade = "E02_idade__Idade_informada" if "E02_idade__Idade_informada" in df_mor.columns else "E02_idade"
df_mor["idade_num"] = pd.to_numeric(df_mor[col_idade], errors="coerce")
df_mor["anos_estudo"] = df_mor[ESC].astype("string").str.strip().str.lower().map(MAPA_ANOS_ESTUDO)
assert df_mor.loc[df_mor["idade_num"] >= 25, "anos_estudo"].notna().all(), "Rótulo de escolaridade sem mapeamento"
df_mor["escolaridade_adulto"] = np.where(df_mor["idade_num"] >= 25, df_mor["anos_estudo"], np.nan)
df_mor["trabalhou_30d"] = parse_binaria(df_mor["I05__Trabalhou_nos_últimos_30_dias"])
df_mor["tem_ctps"] = parse_binaria(df_mor["I17__Possui_CTPS_assinada"])
df_mor["contribui_previdencia"] = parse_binaria(df_mor["I18__Contribui_para_alguma_Previdência_Social_Pública"])

df_dom["tipologia_original"] = df_dom[SUB].astype(str)

VAR_MOR = ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia"]
VAR_DOM = [
    "renda_domiciliar_pc", "possui_automovel", "possui_motocicleta", "possui_bicicleta", "possui_geladeira",
    "possui_freezer", "possui_lavadora", "possui_secadora", "possui_tv", "possui_ar_condicionado",
    "possui_microondas", "comodos_por_pessoa", "banheiros_por_pessoa", "agua_rede", "agua_poco_artesiano",
    "agua_poco_cisterna", "captacao_chuva", "caixa_agua", "solucao_sanitaria_estruturada",
    "acesso_energetico_estruturado", "destino_com_coleta_direta", "rua_pavimentada", "iluminacao_publica",
    "drenagem", "ponto_onibus", "internet_domiciliar", "internet_banda_larga_fixa", "internet_movel",
    "possui_computador", "ivs_df", "qtd_estabelecimentos", "pct_estabelecimentos",
    "agua_nivel_rede", "agua_nivel_ind", "agua_nivel_prec", "esgoto_nivel_rede", "esgoto_nivel_ind", "esgoto_nivel_prec",
    "log_renda_hh", "renda_pc_winsor",
]

def ilr_3partes(p1, p2, p3, n):
    """Razões log-isométricas de uma composição (rede, individual, precária); zeros tratados por suavização aditiva (+0,5 domicílio)."""
    P = np.vstack([p1, p2, p3]).T * n[:, None]
    P = (P + 0.5) / (n[:, None] + 1.5)
    P = P / P.sum(axis=1, keepdims=True)
    if COMP == "clr":
        L = np.log(P); C = L - L.mean(axis=1, keepdims=True)
        return C[:, 0], C[:, 1], C[:, 2]
    adeq = np.sqrt(2 / 3) * np.log(np.sqrt(P[:, 0] * P[:, 1]) / P[:, 2])     # fontes adequadas vs. precárias
    tipo = np.log(P[:, 0] / P[:, 1]) / np.sqrt(2)                            # rede vs. solução individual (sem sentido normativo)
    return adeq, tipo

def derivar_variaveis_modelo(M):
    """Variáveis do segmento que dependem de composições; reaplicável a réplicas de bootstrap."""
    M = M.copy(); n = M["amostra_n"].to_numpy(float)
    if COMP == "clr":
        M["agua_clr_rede"], M["agua_clr_ind"], M["agua_clr_prec"] = ilr_3partes(M["agua_nivel_rede"], M["agua_nivel_ind"], M["agua_nivel_prec"], n)
        M["esgoto_clr_rede"], M["esgoto_clr_ind"], M["esgoto_clr_prec"] = ilr_3partes(M["esgoto_nivel_rede"], M["esgoto_nivel_ind"], M["esgoto_nivel_prec"], n)
    else:
        M["agua_ilr_adequacao"], M["agua_ilr_tipo"] = ilr_3partes(M["agua_nivel_rede"], M["agua_nivel_ind"], M["agua_nivel_prec"], n)
        M["esgoto_ilr_adequacao"], M["esgoto_ilr_tipo"] = ilr_3partes(M["esgoto_nivel_rede"], M["esgoto_nivel_ind"], M["esgoto_nivel_prec"], n)
    M["log_renda_media"] = M["log_renda_hh"]
    return M
matriz_mor_A = media_ponderada_vetorizada(df_mor, CHAVES, VAR_MOR, PESO_MOR)
matriz_dom_A = media_ponderada_vetorizada(df_dom, CHAVES, VAR_DOM, "peso_dom_corrigido")
matriz_dom_A["log_renda_domiciliar_pc"] = np.where(matriz_dom_A["renda_domiciliar_pc"] >= 0, np.log1p(matriz_dom_A["renda_domiciliar_pc"]), np.nan)
counts_A = df_dom.groupby(CHAVES).agg(amostra_n=(ID, "count"), domicilios_expandidos=("peso_dom_corrigido", "sum")).reset_index()

df_matriz_A = counts_A.merge(matriz_dom_A, on=CHAVES, how="left").merge(matriz_mor_A, on=CHAVES, how="left")
df_matriz_A = df_matriz_A[df_matriz_A["amostra_n"] >= 5].reset_index(drop=True)
df_matriz_A = derivar_variaveis_modelo(df_matriz_A)
sub_dom = df_dom.groupby(CHAVES)[SUB].agg(lambda x: x.mode()[0] if not x.mode().empty else "Desconhecido").reset_index(name="tipologia_original")
df_matriz_A = df_matriz_A.merge(sub_dom, on=CHAVES, how="left")
print(f"\nFASE 1 CONCLUÍDA: {len(df_matriz_A)} segmentos territoriais válidos.")

# Checagem de sanidade das variáveis corrigidas
assert df_matriz_A["internet_domiciliar"].notna().all() and df_matriz_A["internet_domiciliar"].nunique() > 1, "internet_domiciliar degenerada"
print(f"   internet_domiciliar: min={df_matriz_A['internet_domiciliar'].min():.3f} | mediana={df_matriz_A['internet_domiciliar'].median():.3f} | max={df_matriz_A['internet_domiciliar'].max():.3f}")

# ==============================================================================
# ETAPA 4 — PCA POR DIMENSÃO (ANÁLISE PARALELA DE MONTE CARLO)
# ==============================================================================
print("\nEtapa 4 — PCA por dimensão (análise paralela de Monte Carlo)...\n")
D1_VARS = (["log_renda_media" if "renda" in V2 else "log_renda_domiciliar_pc", "possui_automovel", "possui_motocicleta"]
           + ([] if "dur" in V2 else ["possui_geladeira"]) + ["possui_freezer", "possui_lavadora"]
           + ([] if "dur" in V2 else ["possui_secadora"]) + ["possui_tv", "possui_ar_condicionado", "possui_microondas"])
D2_VARS = ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia"] + ([] if "ctx" in V2 else ["ivs_df"])
D3_VARS = (["comodos_por_pessoa", "banheiros_por_pessoa"]
           + ((["agua_clr_rede", "agua_clr_ind", "agua_clr_prec"] if COMP == "clr" else ["agua_ilr_adequacao", "agua_ilr_tipo"]) if "agua" in V2 else ["agua_rede", "agua_poco_artesiano", "agua_poco_cisterna"])
           + ["captacao_chuva", "caixa_agua"]
           + ((["esgoto_clr_rede", "esgoto_clr_ind", "esgoto_clr_prec"] if COMP == "clr" else ["esgoto_ilr_adequacao", "esgoto_ilr_tipo"]) if "esgoto" in V2 else ["solucao_sanitaria_estruturada"])
           + ["acesso_energetico_estruturado", "destino_com_coleta_direta", "rua_pavimentada", "iluminacao_publica", "drenagem", "ponto_onibus"]
           + ([] if "ctx" in V2 else ["qtd_estabelecimentos"]) + (["pct_estabelecimentos"] if USAR_PCT_ESTAB else []))
DIMENSOES_REF = {
    "D1_Financeiro_Ativos": D1_VARS,
    "D2_Humano_Trabalho": D2_VARS,
    "D3_Fisico_Habitacional_Territorial": D3_VARS,
    "D4_Tecnologico_Conectividade": ["internet_domiciliar", "internet_banda_larga_fixa", "internet_movel", "possui_computador"],
}

def monte_carlo_parallel_analysis(n_samples, n_features, n_iterations=1000, random_state=SEED):
    rs = np.random.RandomState(random_state)
    sim = np.zeros((n_iterations, n_features))
    for i in range(n_iterations):
        sim[i, :] = np.linalg.eigvalsh(np.corrcoef(rs.normal(size=(n_samples, n_features)), rowvar=False))[::-1]
    return np.percentile(sim, 95, axis=0)

df_componentes_mc = pd.DataFrame(index=df_matriz_A.index)
df_A_mc = df_matriz_A.copy()
df_A_mc["ra_norm"] = df_A_mc[RA].map(limpar_nome_ra_base)
resumo_pca, cargas = [], []

for dim, vars_list in DIMENSOES_REF.items():
    valid = [v for v in vars_list if v in df_matriz_A.columns and df_matriz_A[v].notna().sum() >= 3 and df_matriz_A[v].dropna().nunique() > 1]
    descartadas = sorted(set(vars_list) - set(valid))
    if descartadas: print(f"{dim}: variáveis descartadas por serem constantes/ausentes: {descartadas}")
    X_std = StandardScaler().fit_transform(df_matriz_A[valid].apply(pd.to_numeric, errors="coerce").fillna(0))
    eig_obs = np.linalg.eigvalsh(np.corrcoef(X_std, rowvar=False))[::-1]
    eig_p95 = monte_carlo_parallel_analysis(*X_std.shape)
    k_mc = max(1, int(np.sum(eig_obs > eig_p95)))
    k_kai = max(1, int(np.sum(eig_obs >= 1.0)))
    pca_mc = PCA(n_components=k_mc).fit(X_std)
    scores = pca_mc.transform(X_std)
    for j in range(k_mc):
        nome = f"{dim}_MC_PC{j+1}"
        df_A_mc[nome] = scores[:, j]
        df_componentes_mc[nome] = scores[:, j]
        for idx, v in enumerate(valid):
            c = pca_mc.components_[j, idx]
            cargas.append({"Dimensão": dim, "Componente": f"PC{j+1}", "Variável": v,
                           "Carga Fatorial (Loading)": round(c, 4), "Influência Absoluta": round(abs(c), 4)})
    pca_full = PCA().fit(X_std)
    for j in range(len(eig_obs)):
        resumo_pca.append({"dimensao": dim, "componente": j + 1, "eigenvalue_obs": round(eig_obs[j], 4),
                           "eigenvalue_mc_p95": round(eig_p95[j], 4), "retido_mc": j < k_mc, "retido_kaiser": j < k_kai,
                           "variancia_explicada_pct": round(pca_full.explained_variance_ratio_[j] * 100, 2),
                           "variancia_acumulada_pct": round(np.cumsum(pca_full.explained_variance_ratio_ * 100)[j], 2)})

df_resumo_pca = pd.DataFrame(resumo_pca)
print("=" * 110, "\nRETENÇÃO DE COMPONENTES (MONTE CARLO vs. KAISER)\n", "=" * 110, sep="")
print(df_resumo_pca.to_string(index=False))
df_cargas = pd.DataFrame(cargas).sort_values(["Dimensão", "Componente", "Influência Absoluta"], ascending=[True, True, False])
print("\n" + "=" * 110, "\nCARGAS FATORIAIS DOS COMPONENTES RETIDOS\n", "=" * 110, sep="")
print(df_cargas.to_string(index=False))

# ==============================================================================
# ETAPA 5 — BENCHMARK DE K
# ==============================================================================
comp_cols = df_componentes_mc.columns.tolist()
X_pca = df_componentes_mc.apply(pd.to_numeric).to_numpy()
N, p = X_pca.shape
print(f"\nEtapa 5 — Benchmark de K. Matriz de agrupamento: {N} segmentos x {p} componentes:")
for c in comp_cols: print("  -", c)

def gap_statistic(X, k, n_refs=50):
    disp = KMeans(k, n_init=10, random_state=SEED).fit(X).inertia_
    lo, hi = X.min(axis=0), X.max(axis=0)
    ref = [KMeans(k, n_init=10, random_state=SEED).fit(RNG.uniform(lo, hi, size=X.shape)).inertia_ for _ in range(n_refs)]
    return np.mean(np.log(ref)) - np.log(disp)

linhas_k = []
for k in K_CANDIDATOS:
    lab = KMeans(k, n_init=50, random_state=SEED).fit(X_pca).labels_
    linhas_k.append({"K": k, "Silhouette": round(silhouette_score(X_pca, lab), 4),
                     "Calinski_Harabasz": round(calinski_harabasz_score(X_pca, lab), 2),
                     "Davies_Bouldin": round(davies_bouldin_score(X_pca, lab), 4),
                     "Gap_Statistic": round(gap_statistic(X_pca, k), 4),
                     "BIC": round(GaussianMixture(k, covariance_type="diag", random_state=SEED).fit(X_pca).bic(X_pca), 2),
                     "Menor_cluster_(segmentos)": int(np.bincount(lab).min())})
df_metrics = pd.DataFrame(linhas_k)
print("\n" + "=" * 110, "\nBENCHMARK DE K\n", "=" * 110, sep="")
print(df_metrics.to_string(index=False))

# ==============================================================================
# ETAPA 5.1 / 5.2 — NULOS E ESTABILIDADE POR K
# ==============================================================================
print("\nEtapa 5.1/5.2 — Silhueta contra dois nulos + estabilidade por subamostragem...\n")
def sil_km(Z, k, n_init, rs):
    return silhouette_score(Z, KMeans(k, n_init=n_init, random_state=rs).fit(Z).labels_)

cov, mu = np.cov(X_pca, rowvar=False), X_pca.mean(axis=0)
linhas_nulo = []
for k in K_CANDIDATOS:
    lab = KMeans(k, n_init=50, random_state=SEED).fit(X_pca).labels_
    real = silhouette_score(X_pca, lab)
    nulo_col, nulo_mv = [], []
    for i in range(N_PERM_NULO):
        Z = X_pca.copy()
        for c in range(p): RNG.shuffle(Z[:, c])                                   # nulo 1: destrói todas as correlações
        nulo_col.append(sil_km(Z, k, 5, i))
        nulo_mv.append(sil_km(RNG.multivariate_normal(mu, cov, size=N), k, 5, i))  # nulo 2: preserva a covariância
    nulo_col, nulo_mv = np.array(nulo_col), np.array(nulo_mv)
    aris = []
    for b in range(N_SUBAMOSTRAS):
        idx = RNG.choice(N, int(0.8 * N), replace=False)
        aris.append(adjusted_rand_score(lab[idx], KMeans(k, n_init=10, random_state=b).fit(X_pca[idx]).labels_))
    linhas_nulo.append({"K": k, "Silhueta": round(real, 4),
                        "Nulo_embaralhado_media": round(nulo_col.mean(), 4), "p_embaralhado": round((np.sum(nulo_col >= real) + 1) / (N_PERM_NULO + 1), 4),
                        "Nulo_covariancia_media": round(nulo_mv.mean(), 4), "p_covariancia": round((np.sum(nulo_mv >= real) + 1) / (N_PERM_NULO + 1), 4),
                        "Estabilidade_ARI": round(float(np.mean(aris)), 3)})
df_nulos = pd.DataFrame(linhas_nulo)
print(df_nulos.to_string(index=False))
print("\nLeitura: o nulo 'embaralhado' costuma rejeitar para qualquer K e NÃO discrimina entre eles;"
      " use p_covariancia e a estabilidade (ARI) para comparar candidatos.")

# ==============================================================================
# ETAPA 6 — CLUSTERIZAÇÃO FINAL, CRUZAMENTO COM TIPOLOGIA E KRUSKAL-WALLIS
# ==============================================================================
print(f"\nEtapa 6 — K-Means final (K={K_FINAL}) e cruzamento com a tipologia original...\n")
kmeans_final = KMeans(n_clusters=K_FINAL, n_init=50, random_state=SEED).fit(X_pca)
df_A_mc[CL] = kmeans_final.labels_

# Exportações lidas pelos programas complementares (04 a 12)
np.save("X_pca.npy", X_pca)                                                         # 49 segmentos x 6 componentes
df_A_mc.to_csv("segmentos_completo.csv", index=False)                               # variáveis, componentes e grupo K-Means
df_A_mc[[RA, SUB, "amostra_n", "domicilios_expandidos"]].to_csv("seg_labels.csv", index=False)
obs[[ESTRATO, "N_h", "n_h", "peso_dom_corrigido"]].rename(columns={"peso_dom_corrigido": "peso_N_h_sobre_n_h"}).to_csv("universo_por_estrato.csv", index=False)  # universo e amostra por estrato
col_peso = "domicilios_expandidos"

resumo_clusters = df_A_mc.groupby(CL).agg(n_segmentos=("ra_norm", "size"), domicilios_expandidos=(col_peso, "sum")).reset_index()
print(resumo_clusters.to_string(index=False))
tabela_cruzada = df_A_mc.pivot_table(index=CL, columns="tipologia_original", values=col_peso, aggfunc="sum", fill_value=0)
print("\nDomicílios expandidos por cluster x tipologia original:")
print(tabela_cruzada.round(1).to_string())

print("\nKruskal-Wallis (caráter descritivo, n=49 segmentos):")
variaveis_perfil = ["log_renda_domiciliar_pc", "escolaridade_adulto", "trabalhou_30d", "banheiros_por_pessoa", "comodos_por_pessoa",
                    "agua_rede", "internet_domiciliar", "internet_banda_larga_fixa", "possui_computador", "ivs_df", "qtd_estabelecimentos"]
def fmt_p(v): return "N/A" if pd.isna(v) else (f"{v:.4e}" if v < 1e-4 else f"{v:.4f}")
res_kw = []
for var in variaveis_perfil:
    if var not in df_A_mc.columns: continue
    grupos = [g[var].dropna().values for _, g in df_A_mc.groupby(CL)]
    grupos = [g for g in grupos if len(g) > 0]
    if len(grupos) > 1 and len(np.unique(np.concatenate(grupos))) > 1:
        H, pv = kruskal(*grupos)
        res_kw.append({"Variável": var, "Estatística H": round(H, 4), "p-valor": fmt_p(pv), "Significativo (5%)": "Sim" if pv < 0.05 else "Não"})
print(pd.DataFrame(res_kw).to_string(index=False))

# ==============================================================================
# ETAPA 6.1 — SENSIBILIDADE DA PARTIÇÃO À CONSTRUÇÃO DE `escolaridade_adulto`
# ==============================================================================
# A PDAD traz a escolaridade por categoria. Os níveis completos usam a duração oficial (9, 12, 16 anos, como no IDH/PNUD);
# apenas os níveis INCOMPLETOS são aproximações. Aqui varia-se esses valores e o corte de idade, refaz-se o PCA da D2
# (mesma regra de Monte Carlo) e compara-se a partição K_FINAL com a de referência pelo ARI.
print("\nEtapa 6.1 — Sensibilidade a `escolaridade_adulto` (valores dos níveis incompletos e corte de idade)...\n")
_base = dict(MAPA_ANOS_ESTUDO)
def _mapa(f, m, sup): return {**_base, "fundamental incompleto": f, "médio incompleto": m, "superior incompleto": sup}
CENARIOS_ESC = {
    "Referência (4,5 / 10,5 / 14; idade >= 25)": (_base, 25),
    "Baixo (2 / 10 / 13)": (_mapa(2, 10, 13), 25),
    "Alto (6 / 11 / 15)": (_mapa(6, 11, 15), 25),
    "Extremo baixo (1 / 9,5 / 12,5)": (_mapa(1, 9.5, 12.5), 25),
    "Extremo alto (8 / 11,5 / 15,5) [implausível]": (_mapa(8, 11.5, 15.5), 25),
    "Ordinal 0-6 (sem converter em anos)": ({k: i for i, k in enumerate(_base)}, 25),
    "Idade >= 20": (_base, 20),
    "Idade >= 18": (_base, 18),
}
_d2_vars = DIMENSOES_REF["D2_Humano_Trabalho"]
_cols_outras = [c for c in comp_cols if not c.startswith("D2_Humano_Trabalho")]
linhas_sens = []
particoes_cen = []
for nome_c, (mapa_c, idade_c) in CENARIOS_ESC.items():
    _m = df_mor.copy()
    _anos = _m[ESC].astype("string").str.strip().str.lower().map(mapa_c)
    _m["escolaridade_adulto"] = np.where(_m["idade_num"] >= idade_c, _anos, np.nan)
    _esc = media_ponderada_vetorizada(_m, CHAVES, ["escolaridade_adulto"], PESO_MOR)
    _A = df_matriz_A[CHAVES + [v for v in _d2_vars if v != "escolaridade_adulto"]].merge(_esc, on=CHAVES, how="left")
    _X = StandardScaler().fit_transform(_A[_d2_vars].apply(pd.to_numeric, errors="coerce").fillna(0))
    _k = max(1, int(np.sum(np.linalg.eigvalsh(np.corrcoef(_X, rowvar=False))[::-1] > monte_carlo_parallel_analysis(*_X.shape))))
    _Z = np.column_stack([df_componentes_mc[_cols_outras].to_numpy(), PCA(n_components=_k).fit_transform(_X)])
    _lab = KMeans(K_FINAL, n_init=50, random_state=SEED).fit(_Z).labels_
    particoes_cen.append(_lab)
    linhas_sens.append({"Cenário": nome_c, "Comp. totais": _Z.shape[1], "Tamanhos": sorted(np.bincount(_lab), reverse=True),
                        "Silhueta": round(silhouette_score(_Z, _lab), 3),
                        "ARI vs referência": round(adjusted_rand_score(kmeans_final.labels_, _lab), 3)})
df_sens_esc = pd.DataFrame(linhas_sens)
print(df_sens_esc.to_string(index=False))
print("\nRedação: ARI próximo de 1 = a partição não depende do valor atribuído aos níveis incompletos. "
      "O cenário 'extremo alto' aproxima o nível incompleto do completo e contradiz a definição das categorias.")

# ==============================================================================
# ETAPA 6.2 — MATRIZ DE CO-ASSOCIAÇÃO (CONSENSO): QUAIS SEGMENTOS SÃO "NÚCLEO" E QUAIS SÃO "FRONTEIRA"
# ==============================================================================
# Duas fontes de variação: (a) subamostragem de 80% dos segmentos e (b) os 8 cenários de construção de `escolaridade_adulto`.
# Para cada par de segmentos mede-se a fração de partições em que ficam no mesmo cluster. O escore de cada segmento é a
# co-associação média com os demais segmentos do SEU cluster de referência. Limiar de 0,75 é convenção (não é teste estatístico).
print("\nEtapa 6.2 — Co-associação (consenso) entre subamostras e cenários...\n")
LIM_NUCLEO = 0.75
N_ALERTA_SEG = 10   # segmentos com menos domicílios na amostra são marcados com * nos mapas
ref = kmeans_final.labels_
same_sub, both_sub = np.zeros((N, N)), np.zeros((N, N))
for b in range(N_SUBAMOSTRAS):
    idx = np.sort(RNG.choice(N, int(0.8 * N), replace=False))
    lb = KMeans(K_FINAL, n_init=10, random_state=1000 + b).fit(X_pca[idx]).labels_
    same_sub[np.ix_(idx, idx)] += (lb[:, None] == lb[None, :])
    both_sub[np.ix_(idx, idx)] += 1
co_sub = np.divide(same_sub, both_sub, out=np.zeros_like(same_sub), where=both_sub > 0)
co_cen = np.mean([(pl[:, None] == pl[None, :]).astype(float) for pl in particoes_cen], axis=0)

def _rotulo(r):
    sub = (r[SUB].replace("Características Urbanas", "Moradia e Lazer").replace("Agricultura ", "Agr. ")
           .replace("Assentamento/Agrovila", "Assentamento"))
    return f"{r[RA]} — {sub}"
rotulos = [_rotulo(r) for _, r in df_A_mc.iterrows()]
linhas_co = []
for i in range(N):
    mates = [j for j in range(N) if ref[j] == ref[i] and j != i]
    e_sub = float(np.mean(co_sub[i, mates])) if mates else np.nan
    e_cen = float(np.mean(co_cen[i, mates])) if mates else np.nan
    linhas_co.append({"Segmento": rotulos[i], "Cluster_ref": int(ref[i]), "n_amostra": int(df_A_mc.loc[i, "amostra_n"]),
                      "Escore_subamostras": round(e_sub, 2), "Escore_cenarios": round(e_cen, 2),
                      "Classe": "Núcleo" if (min(e_sub, e_cen) >= LIM_NUCLEO if mates else False) else "Fronteira"})
df_co = pd.DataFrame(linhas_co).sort_values(["Cluster_ref", "Escore_subamostras"], ascending=[True, False]).reset_index(drop=True)
with pd.option_context("display.max_rows", None, "display.width", 200):
    print(df_co.to_string(index=False))
print("\nResumo por cluster de referência:")
print(df_co.groupby("Cluster_ref").agg(n_segmentos=("Segmento", "size"), nucleo=("Classe", lambda x: int((x == "Núcleo").sum())),
                                       fronteira=("Classe", lambda x: int((x == "Fronteira").sum())),
                                       escore_sub_medio=("Escore_subamostras", "mean"), escore_cen_medio=("Escore_cenarios", "mean")).round(2).to_string())
print(f"\nSegmentos 'núcleo' (escore >= {LIM_NUCLEO} nas duas fontes): {(df_co['Classe'] == 'Núcleo').sum()} de {N}")
df_co.to_csv(FIG_DIR / "coassociacao_segmentos.csv", index=False)
_ord = np.argsort(ref, kind="stable")
fig, axes = plt.subplots(1, 2, figsize=(20, 9))
for ax, M_, tit in ((axes[0], co_sub, "Subamostras (80%)"), (axes[1], co_cen, "Cenários de escolaridade")):
    sns.heatmap(M_[np.ix_(_ord, _ord)], cmap="viridis", vmin=0, vmax=1, ax=ax, xticklabels=False,
                yticklabels=[rotulos[i] for i in _ord], cbar_kws={"label": "co-associação"})
    ax.set_title(f"Co-associação — {tit}", fontweight="bold"); ax.tick_params(axis="y", labelsize=6)
salvar("coassociacao_segmentos")

# ==============================================================================
# ETAPA 7 — LISA (BASELINE E SENSIBILIDADE)
# ==============================================================================
print("\nEtapa 7 — Análise espacial local (LISA)...\n")

def carregar_shapefile():
    if not PATH_SHP_ZIP.exists():
        print(f"Shapefile não encontrado em {PATH_SHP_ZIP}: cenários k-NN e mapas serão pulados.")
        return None
    pasta = "shape_df_extracao"
    if not os.path.exists(pasta):
        with zipfile.ZipFile(PATH_SHP_ZIP, "r") as z: z.extractall(pasta)
    shp = glob.glob(f"{pasta}/**/*.shp", recursive=True)
    if not shp: return None
    g = gpd_read(shp[0])
    col = next((c for c in g.columns if c.lower() in ["ra", "nome", "ra_nome", "nome_ra", "regiao", "nm_ra"]), g.columns[0])
    g["ra_norm"] = g[col].apply(limpar_nome_ra_base).replace({"nacleo bandeirante": "nucleo bandeirante"})
    return g

def gpd_read(path):
    import geopandas as gpd
    return gpd.read_file(path)

gdf_df = carregar_shapefile()

def media_ra_ponderada(df, cols):
    """Média por RA ponderada pelos domicílios expandidos de cada segmento (mesma regra dos mapas e do Moran/LISA)."""
    return df.groupby("ra_norm").apply(lambda g: pd.Series({c: np.average(g[c], weights=g[col_peso]) for c in cols}))
df_ra = media_ra_ponderada(df_A_mc, ["D4_Tecnologico_Conectividade_MC_PC1", "internet_domiciliar"]).join(df_A_mc.groupby("ra_norm")[col_peso].sum()).reset_index()
RAs = df_ra["ra_norm"].tolist()
assert df_ra[["D4_Tecnologico_Conectividade_MC_PC1", "internet_domiciliar"]].notna().all().all(), "NaN no vetor do LISA"  # [CORREÇÃO 6]

VIZINHOS_DF = {
    "planaltina": ["sobradinho", "sobradinho ii", "paranoa", "fercal"],
    "sobradinho": ["planaltina", "sobradinho ii", "fercal", "paranoa", "varjao", "itapoa", "lago norte"],
    "sobradinho ii": ["sobradinho", "fercal", "brazlandia"], "fercal": ["sobradinho", "sobradinho ii", "planaltina", "brazlandia"],
    "paranoa": ["planaltina", "sobradinho", "itapoa", "sao sebastiao", "jardim botanico", "lago norte"],
    "itapoa": ["paranoa", "sobradinho"], "sao sebastiao": ["paranoa", "jardim botanico", "santa maria"],
    "jardim botanico": ["sao sebastiao", "paranoa", "lago sul", "santa maria"],
    "brazlandia": ["sobradinho ii", "fercal", "ceilandia", "taguatinga"],
    "ceilandia": ["brazlandia", "sol nascente por do sol", "samambaia", "taguatinga"],
    "sol nascente por do sol": ["ceilandia", "samambaia", "taguatinga"],
    "samambaia": ["ceilandia", "sol nascente por do sol", "taguatinga", "riacho fundo ii", "recanto das emas"],
    "recanto das emas": ["samambaia", "riacho fundo ii", "gama", "santa maria"],
    "gama": ["recanto das emas", "riacho fundo ii", "park way", "santa maria"],
    "santa maria": ["gama", "park way", "jardim botanico", "sao sebastiao", "recanto das emas"],
    "riacho fundo ii": ["samambaia", "recanto das emas", "gama", "riacho fundo", "park way"],
    "riacho fundo": ["riacho fundo ii", "park way", "nucleo bandeirante", "guara", "taguatinga"],
    "park way": ["gama", "santa maria", "jardim botanico", "lago sul", "nucleo bandeirante", "riacho fundo", "riacho fundo ii", "guara"],
    "nucleo bandeirante": ["riacho fundo", "park way", "guara", "candangolandia"],
    "guara": ["riacho fundo", "park way", "nucleo bandeirante", "candangolandia", "sia", "taguatinga", "aguas claras", "cruzeiro", "plano piloto"],
    "taguatinga": ["ceilandia", "sol nascente por do sol", "samambaia", "riacho fundo", "guara", "aguas claras", "vicente pires", "brazlandia"],
    "aguas claras": ["taguatinga", "guara", "vicente pires"], "vicente pires": ["taguatinga", "aguas claras"],
    "lago sul": ["park way", "jardim botanico", "plano piloto", "cruzeiro", "lago norte", "candangolandia"],
    "lago norte": ["sobradinho", "paranoa", "plano piloto", "lago sul"],
    "plano piloto": ["guara", "cruzeiro", "lago sul", "lago norte", "sobradinho", "varjao", "sia"],
    "cruzeiro": ["guara", "plano piloto", "lago sul", "sia"], "varjao": ["sobradinho", "plano piloto"],
    "sia": ["guara", "plano piloto"], "candangolandia": ["nucleo bandeirante", "guara", "lago sul"],
}
adj = {ra: [v for v in VIZINHOS_DF.get(ra, []) if v in RAs] for ra in RAs}
for k in adj:
    if not adj[k]: adj[k] = [RAs[0]]
W_Queen_manual = libpysal.weights.W(adj, id_order=RAs); W_Queen_manual.transform = "R"   # versão digitada (só se não houver shapefile)
W_Queen = W_Queen_manual
FONTE_QUEEN = "lista manual de vizinhos (sem shapefile)"
if gdf_df is not None:
    # Contiguidade Queen calculada a partir dos polígonos das RAs (substitui a lista manual, que divergia do mapa em ~19 vínculos)
    _gq = gdf_df[gdf_df["ra_norm"].isin(RAs)].drop_duplicates("ra_norm").set_index("ra_norm").loc[RAs].reset_index()
    W_Queen = libpysal.weights.Queen.from_dataframe(_gq, idVariable="ra_norm", use_index=False)
    assert W_Queen.id_order == RAs and not W_Queen.islands, f"Queen com ilhas: {W_Queen.islands}"
    _dif = sum(set(W_Queen.neighbors[r]) != set(W_Queen_manual.neighbors[r]) for r in RAs)
    print(f"OK Queen calculada do shapefile ({len(RAs)} RAs, sem ilhas). RAs cujos vizinhos diferem da lista manual: {_dif}.")
    W_Queen.transform = "R"; FONTE_QUEEN = "shapefile das RAs"
    import json; json.dump({k: list(v) for k, v in W_Queen.neighbors.items()}, open("vizinhos_queen.json", "w"), ensure_ascii=False)   # usado pelo programa 18 (join count)

W_KNN = None
if gdf_df is not None:
    import geopandas as gpd
    cent = gdf_df.set_index("ra_norm")["geometry"].centroid.to_dict()
    gdf_ra = gpd.GeoDataFrame(df_ra.assign(geometry=df_ra["ra_norm"].map(cent)), geometry="geometry")
    W_KNN = KNN.from_dataframe(gdf_ra, k=4); W_KNN.transform = "R"

cenarios = [
    {"nome": "BASELINE DO TCC (Queen, PCA D4, FDR)", "w": W_Queen, "var": "D4_Tecnologico_Conectividade_MC_PC1", "fdr": True, "is_baseline": True},
    {"nome": "SENSIBILIDADE DE VIZINHANÇA (k-NN, PCA D4, FDR)", "w": W_KNN, "var": "D4_Tecnologico_Conectividade_MC_PC1", "fdr": True},
    {"nome": "EXPLORATÓRIO (Queen, PCA D4, p<0,05 sem FDR)", "w": W_Queen, "var": "D4_Tecnologico_Conectividade_MC_PC1", "fdr": False},
    {"nome": "INFRAESTRUTURA BRUTA (Queen, internet_domiciliar, FDR)", "w": W_Queen, "var": "internet_domiciliar", "fdr": True},
    {"nome": "FRICÇÃO RURURBANA MÚLTIPLA (k-NN, PCA D4, p<0,05 sem FDR)", "w": W_KNN, "var": "D4_Tecnologico_Conectividade_MC_PC1", "fdr": False},
]
MAPA_Q = {1: "Q1 (Alto-Alto): Polo de Alta Conectividade", 2: "Q2 (Baixo-Alto): Outlier Tecnológico",
          3: "Q3 (Baixo-Baixo): Polo de Exclusão Digital", 4: "Q4 (Alto-Baixo): Outlier Tecnológico"}
lisa_res = {}
for cen in cenarios:
    print("=" * 100, f"\n{cen['nome']}\n", "=" * 100, sep="")
    if cen["w"] is None:
        print("Pulado (requer shapefile para a matriz k-NN).\n"); continue
    y = df_ra[cen["var"]].to_numpy(dtype=float)
    lm = esda.moran.Moran_Local(y, cen["w"], permutations=999, seed=SEED)
    mask = multipletests(lm.p_sim, alpha=0.05, method="fdr_bh")[0] if cen["fdr"] else lm.p_sim < 0.05
    lisa_res[cen["nome"]] = (lm.q, lm.p_sim, mask)
    if cen.get("is_baseline"):
        df_ra["LISA_Q"], df_ra["LISA_Significativo"] = lm.q, mask
    if not mask.any():
        print("Nenhum polo significativo neste cenário.\n"); continue
    t = pd.DataFrame({"Região Administrativa (RA)": [r.title() for r in np.array(RAs)[mask]],
                      "Quadrante": [MAPA_Q[q] for q in lm.q[mask]], "p-Valor (Simulado)": lm.p_sim[mask]})
    print(t.to_string(index=False), "\n")

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
NUM_RA = numerar_ras(RAs)

# ---- Tabela LaTeX com a numeração das RAs (usada nos mapas) ----
_ras_ord = sorted(NUM_RA, key=NUM_RA.get); _n = len(_ras_ord); _lin = (_n + 2) // 3
_linhas = []
for _i in range(_lin):
    _cel = []
    for _c in range(3):
        _k = _c * _lin + _i
        _cel.append(f"{_k + 1} & {NOME_RA.get(_ras_ord[_k], _ras_ord[_k].title()).replace('&', chr(92) + '&')}" if _k < _n else " & ")
    _linhas.append(" & ".join(_cel) + r" \\")
(FIG_DIR / "Tabela_num_ras.tex").write_text("\\begin{table}[htbp]\n\\centering\n\\caption{Numeração das Regiões Administrativas nos mapas}\\label{tab:num_ras}\n"
    "\\begin{tabular}{r l r l r l}\n\\hline\n" + "\n".join(_linhas) + "\n\\hline\n\\end{tabular}\n\\par\\footnotesize Fonte: elaboração própria. Ordem alfabética.\n\\end{table}\n", encoding="utf-8")

# ---- Moran global (RA) e mapas LISA --------------------------------------------
print("=" * 100, "\nAUTOCORRELAÇÃO ESPACIAL GLOBAL (Moran's I, 999 permutações; unidade = RA, média dos segmentos)\n", "=" * 100, sep="")
df_ra_ext = media_ra_ponderada(df_A_mc, comp_cols + ["internet_domiciliar", "agua_rede", "rua_pavimentada", "solucao_sanitaria_estruturada",
                                                    "log_renda_domiciliar_pc", "escolaridade_adulto", "ivs_df"]).loc[RAs]
linhas_mi = []
for var in df_ra_ext.columns:
    for nome_w, w in (("Queen", W_Queen), ("k-NN (k=4)", W_KNN)):
        if w is None: continue
        np.random.seed(SEED); mi = esda.moran.Moran(df_ra_ext[var].to_numpy(float), w, permutations=999)
        linhas_mi.append({"Variável (média por RA)": var, "Vizinhança": nome_w, "I de Moran": round(mi.I, 3), "p (simulado)": round(mi.p_sim, 3)})
df_moran = pd.DataFrame(linhas_mi); print(df_moran.to_string(index=False))
df_moran.to_csv(FIG_DIR / "moran_global.csv", index=False)

if gdf_df is not None:
    COR_Q = {1: "#d7191c", 2: "#abd9e9", 3: "#2c7bb6", 4: "#fdae61"}
    NOME_Q = {1: "Alto-Alto", 2: "Baixo-Alto", 3: "Baixo-Baixo", 4: "Alto-Baixo"}
    import matplotlib.patches as mpatches
    def _rotular(ax, g):
        rotular_numeros(ax, g, NUM_RA, fs=10)
    chaves = [k for k in lisa_res if "BASELINE" in k or "EXPLORAT" in k or "FRICÇÃO" in k]
    TITULOS = {"BASELINE": "(a) Queen, com correção FDR", "EXPLORAT": "(b) Queen, sem correção",
               "FRICÇÃO": "(c) k-NN (k = 4), sem correção"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8))
    axes = axes.ravel()
    for ax, k in zip(axes, chaves):
        q, psim, mask = lisa_res[k]
        gdf_df.plot(ax=ax, color="white", edgecolor="#bbbbbb", linewidth=.4)
        base = gdf_df[gdf_df["ra_norm"].isin(RAs)].copy().set_index("ra_norm").loc[RAs].reset_index()
        base.plot(ax=ax, color="#e6e6e6", edgecolor="black", linewidth=.6)
        sig = base.assign(q=q, sig=mask)
        for qq in (1, 2, 3, 4):
            sub = sig[(sig["sig"]) & (sig["q"] == qq)]
            if len(sub): sub.plot(ax=ax, color=COR_Q[qq], edgecolor="black", linewidth=.8)
        rotular_numeros(ax, base, NUM_RA, fs=9, dmin=0.06)
        ax.set_axis_off()
        ax.set_title(next(v for kk, v in TITULOS.items() if kk in k), fontsize=13, fontweight="bold")
    for ax in axes[len(chaves):-1]: ax.set_axis_off()
    axes[-1].set_axis_off()
    axes[-1].legend(handles=[mpatches.Patch(color=COR_Q[i], label=NOME_Q[i]) for i in (1, 2, 3, 4)] + [mpatches.Patch(color="#e6e6e6", label="Não significativo")],
                    loc="center", fontsize=12, title="LISA do componente de conectividade (D4)", title_fontsize=12, frameon=False)
    salvar("mapa_lisa")

    # Heterogeneidade intra-RA: quantos clusters distintos convivem dentro de cada RA
    het = df_A_mc.groupby("ra_norm").agg(n_clusters=(CL, "nunique"), n_segmentos=("ra_norm", "size")).reset_index()
    gh = gdf_df.merge(het, on="ra_norm", how="left")
    gh["n_clusters"] = gh["n_clusters"].map(lambda v: f"{int(v)}" if pd.notna(v) else np.nan)
    fig, ax = plt.subplots(figsize=(11, 9))
    gdf_df.plot(ax=ax, color="white", edgecolor="#bbbbbb", linewidth=.4)
    gh.dropna(subset=["n_clusters"]).plot(column="n_clusters", ax=ax, cmap="OrRd", edgecolor="black", linewidth=.7, categorical=True, legend=True,
                                          legend_kwds={"title": "Nº de clusters\nna RA", "loc": "lower left"})
    _rotular(ax, gdf_df); ax.set_axis_off(); ax.set_title("Heterogeneidade intra-RA: clusters distintos entre os segmentos da RA", fontsize=13, fontweight="bold")
    salvar("mapa_heterogeneidade_ra")
    het.sort_values("n_clusters", ascending=False).to_csv(FIG_DIR / "heterogeneidade_ra.csv", index=False)
    print("\nRAs por nº de clusters distintos entre seus segmentos:"); print(het.sort_values(["n_clusters", "n_segmentos"], ascending=False).to_string(index=False))

# ==============================================================================
# GRÁFICOS DOS CLUSTERS
# ==============================================================================
sns.set_theme(style="whitegrid")
tabela_cruzada.plot(kind="bar", stacked=True, figsize=(12, 7), colormap="viridis", edgecolor="black")
plt.title(f"Composição dos grupos (K={K_FINAL}) por tipologia original", fontsize=16, fontweight="bold")
plt.xlabel("Grupo"); plt.ylabel("Domicílios expandidos"); plt.xticks(range(K_FINAL), [NOMES_GRUPOS[i].replace(" — ", "\n", 1) for i in range(K_FINAL)], rotation=0)
plt.legend(title="Tipologia original", bbox_to_anchor=(1.05, 1), loc="upper left"); salvar("clusters_composicao_tipologia")

plt.figure(figsize=(10, 6))
sns.heatmap(tabela_cruzada, annot=True, fmt=".0f", cmap="YlGnBu", linewidths=.5, cbar_kws={"label": "Domicílios"})
plt.title("Domicílios por Cluster e Tipologia"); salvar("clusters_heatmap_tipologia")

# ==============================================================================
# MAPA (se o shapefile estiver disponível)
# ==============================================================================
dom_cluster_ra = (df_A_mc.groupby(["ra_norm", CL])[col_peso].sum().reset_index()
                  .sort_values(col_peso, ascending=False).drop_duplicates("ra_norm"))
if gdf_df is not None:
    import matplotlib as mpl, matplotlib.patches as mpatches
    ras_pdad, ras_shp = set(df_A_mc["ra_norm"]), set(gdf_df["ra_norm"])
    print("RAs da PDAD sem geometria:", sorted(ras_pdad - ras_shp) or "nenhuma")
    gdf_mapa = gdf_df.merge(dom_cluster_ra[["ra_norm", CL]], on="ra_norm", how="left")
    fig, ax = plt.subplots(figsize=(11, 7.5))
    gdf_df.plot(ax=ax, color="lightgrey", edgecolor="white", linewidth=.5)
    gdf_mapa.dropna(subset=[CL]).plot(column=CL, ax=ax, cmap="viridis", edgecolor="black", linewidth=.8, legend=False)
    cl_u = sorted(dom_cluster_ra[CL].unique()); cm = mpl.colormaps["viridis"].resampled(max(len(cl_u), 2))
    ax.legend(handles=[mpatches.Patch(color=cm(i / max(len(cl_u) - 1, 1)), label=NOMES_GRUPOS[int(c)]) for i, c in enumerate(cl_u)],
              title="Grupo predominante", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=12, title_fontsize=12)
    rotular_numeros(ax, gdf_df, NUM_RA, fs=11, dmin=0.045)
    plt.title("Tipologia Territorial: Grupo dominante por RA", fontsize=15, fontweight="bold"); ax.set_axis_off(); salvar("mapa_cluster_dominante")

    # ---- Painéis por subpopulação: cada RA colorida pelo cluster do SEU segmento (unidade real da análise) ----
    CORES_CL = {0: "#440154", 1: "#3b528b", 2: "#21918c", 3: "#5ec962", 4: "#fde725"}
    _cmap = df_co.set_index("Segmento")["Classe"]
    df_A_mc["classe_coassoc"] = [_cmap[r] for r in rotulos]
    subs = [x for x in ["Agricultura Familiar", "Agricultura Empresarial", "Assentamento/Agrovila", "Características Urbanas"] if x in set(df_A_mc[SUB])]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6)); axes = axes.ravel()
    for ax, sub in zip(axes, subs):
        gdf_df.plot(ax=ax, color="white", edgecolor="#bbbbbb", linewidth=.4)
        sg = df_A_mc[df_A_mc[SUB] == sub]
        sem = gdf_df[gdf_df["ra_norm"].isin(set(RAs) - set(sg["ra_norm"]))]
        if len(sem): sem.plot(ax=ax, color="#eeeeee", edgecolor="#999999", linewidth=.5)
        for _, r in sg.iterrows():
            geom = gdf_df[gdf_df["ra_norm"] == r["ra_norm"]]
            if geom.empty: continue
            fr = r["classe_coassoc"] == "Fronteira"
            geom.plot(ax=ax, facecolor=CORES_CL[int(r[CL])], edgecolor="black" if not fr else "white", linewidth=.8, hatch="///" if fr else None)
            if r["amostra_n"] < N_ALERTA_SEG:
                c = geom.geometry.iloc[0].representative_point(); ax.annotate("*", (c.x, c.y), fontsize=16, ha="center", color="red", fontweight="bold", zorder=7)
        rotular_numeros(ax, gdf_df[gdf_df["ra_norm"].isin(sg["ra_norm"])], NUM_RA, fs=10, dmin=0.06); ax.set_axis_off()
        ax.set_title(f"{sub.replace('Características Urbanas', 'Moradia e Lazer')} ({len(sg)} segmentos)", fontsize=14, fontweight="bold")
    h = [mpatches.Patch(color=CORES_CL[i], label=NOMES_GRUPOS[i]) for i in range(K_FINAL)]
    h += [mpatches.Patch(facecolor="white", edgecolor="black", hatch="///", label="Segmento de fronteira (co-associação < 0,75)"),
          mpatches.Patch(color="#eeeeee", label="RA sem segmento desta subpopulação (n < 5)"), plt.Line2D([], [], marker="*", color="red", linestyle="", markersize=12, label=f"n < {N_ALERTA_SEG} domicílios")]
    fig.legend(handles=h, loc="lower center", ncol=3, fontsize=11)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.95, bottom=0.10, hspace=0.10, wspace=0.0)
    plt.savefig(FIG_DIR / "mapa_clusters_por_subpopulacao.png", dpi=200, bbox_inches="tight", pad_inches=0.05); plt.show(); plt.close()

    # ---- Janela ampliada: região central, onde os polígonos são pequenos ----
    CENTRAL = ["taguatinga", "aguas claras", "nucleo bandeirante", "riacho fundo", "riacho fundo ii", "park way", "scia estrutural",
               "vicente pires", "samambaia", "sol nascente por do sol", "guara", "candangolandia", "sia"]
    zb = gdf_df[gdf_df["ra_norm"].isin(CENTRAL)].total_bounds; mx, my = (zb[2] - zb[0]) * .06, (zb[3] - zb[1]) * .06
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    for ax, tipo in zip(axes, ("cluster", "lisa")):
        gdf_df.plot(ax=ax, color="white", edgecolor="#bbbbbb", linewidth=.4)
        if tipo == "cluster":
            for _, r in dom_cluster_ra.iterrows():
                gdf_df[gdf_df["ra_norm"] == r["ra_norm"]].plot(ax=ax, color=CORES_CL[int(r[CL])], edgecolor="black", linewidth=.8)
            ax.set_title("Grupo dominante (domicílios expandidos) — região central", fontweight="bold")
            ax.legend(handles=[mpatches.Patch(color=CORES_CL[i], label=NOMES_GRUPOS[i]) for i in range(K_FINAL)], loc="lower left", fontsize=9)
        else:
            k = [x for x in lisa_res if "FRICÇÃO" in x][0]; q, psim, mask = lisa_res[k]
            gdf_df[gdf_df["ra_norm"].isin(RAs)].plot(ax=ax, color="#e6e6e6", edgecolor="black", linewidth=.6)
            for ra, qq, m in zip(RAs, q, mask):
                if m: gdf_df[gdf_df["ra_norm"] == ra].plot(ax=ax, color=COR_Q[int(qq)], edgecolor="black", linewidth=.9)
            ax.set_title("LISA D4 (k-NN, p<0,05 sem FDR) — região central", fontweight="bold")
            ax.legend(handles=[mpatches.Patch(color=COR_Q[i], label=NOME_Q[i]) for i in (1, 2, 3, 4)] + [mpatches.Patch(color="#e6e6e6", label="Não significativo")], loc="lower left", fontsize=9)
        _jan = gdf_df[gdf_df["ra_norm"].isin(RAs) & gdf_df.geometry.representative_point().apply(
            lambda c: zb[0] - mx < c.x < zb[2] + mx and zb[1] - my < c.y < zb[3] + my)]
        rotular_numeros(ax, _jan, NUM_RA, fs=11, dmin=0.06)
        ax.set_xlim(zb[0] - mx, zb[2] + mx); ax.set_ylim(zb[1] - my, zb[3] + my); ax.set_axis_off()
    salvar("mapa_zoom_central")

# ==============================================================================
# ETAPA 8 — PERFIL DOS CLUSTERS E GAPS POR RA
# ==============================================================================
print("\nEtapa 8 — Perfil dos clusters e diagnóstico por RA...\n")
variaveis_diag = ["renda_domiciliar_pc", "escolaridade_adulto", "internet_banda_larga_fixa", "possui_computador",
                  "agua_rede", "solucao_sanitaria_estruturada", "rua_pavimentada"]
def media_pond(df, col, w):
    d = df[[col, w]].dropna()
    return np.nan if d.empty or d[w].sum() == 0 else np.average(d[col], weights=d[w])

perfil_cl = pd.DataFrame({c: {v: media_pond(g, v, col_peso) for v in variaveis_diag} for c, g in df_A_mc.groupby(CL)}).T
perfil_cl.index.name = CL
perfil_ra = pd.DataFrame({r: {v: media_pond(g, v, col_peso) for v in variaveis_diag} for r, g in df_A_mc.groupby("ra_norm")}).T
print("PERFIL MÉDIO DOS CLUSTERS\n", perfil_cl.round(4).to_string(), "\n")
print("PERFIL MÉDIO DAS RAs\n", perfil_ra.round(4).to_string(), "\n")

diag = perfil_ra.join(dom_cluster_ra.set_index("ra_norm")[CL], how="inner")
gaps = pd.DataFrame({r: {v: (row[v] / perfil_cl.loc[int(row[CL]), v] - 1) if pd.notna(row[v]) and perfil_cl.loc[int(row[CL]), v] > 0 else np.nan
                         for v in variaveis_diag} for r, row in diag.iterrows()}).T
print("DIAGNÓSTICO DE VULNERABILIDADES (gap < -5% vs. cluster dominante)")
for r in gaps.index:
    d = gaps.loc[r].dropna(); d = d[d < -0.05].sort_values()
    c = int(diag.loc[r, CL])
    if d.empty: print(f"{r.title()} (Cluster {c}): desempenho equilibrado ou superior à média do grupo.")
    else:
        print(f"{r.title()} (Cluster {c}): deficiências relativas em:")
        for v, g in d.items(): print(f"   - {v}: {g * 100:.1f}% abaixo do esperado para o Cluster {c}")

# Diferença ABSOLUTA em relação à média do cluster dominante (pp para proporções; R$ para renda; anos para escolaridade).
# A diferença relativa explodia quando a média do cluster era muito baixa (ex.: pavimentação no cluster 4 = 2%).
UNID = {"renda_domiciliar_pc": ("R$", 1), "escolaridade_adulto": ("anos", 1)}
dif = pd.DataFrame({r: {v: (row[v] - perfil_cl.loc[int(row[CL]), v]) * (1 if v in UNID else 100) if pd.notna(row[v]) else np.nan
                        for v in variaveis_diag} for r, row in diag.iterrows()}).T
z = dif / dif.std()                                    # cor = diferença em desvios-padrão entre RAs (comparável entre colunas)
ordem = diag.loc[dif.index, CL].sort_values(kind="stable").index
dif, z = dif.loc[ordem], z.loc[ordem]
rotulo_y = [f"{i.title()} (C{int(diag.loc[i, CL])})" for i in dif.index]
rotulo_x = [f"{v}\n({UNID[v][0]})" if v in UNID else f"{v}\n(pp)" for v in variaveis_diag]
ann = dif.map(lambda x: "" if pd.isna(x) else f"{x:+.0f}" if abs(x) >= 10 else f"{x:+.1f}")
plt.figure(figsize=(14, 10))
sns.heatmap(z, cmap="RdYlGn", center=0, vmin=-2.5, vmax=2.5, annot=ann, fmt="", linewidths=.5, xticklabels=rotulo_x, yticklabels=rotulo_y,
            cbar_kws={"label": "Diferença em desvios-padrão entre RAs (cor)"})
plt.title("Desempenho da RA vs. média do seu cluster dominante (valores: diferença absoluta)"); plt.xticks(rotation=45, ha="right")
salvar("desempenho_ra_vs_cluster")
dif.round(2).to_csv(FIG_DIR / "desempenho_ra_vs_cluster_diferencas.csv")
print("\nAnálise concluída.")
