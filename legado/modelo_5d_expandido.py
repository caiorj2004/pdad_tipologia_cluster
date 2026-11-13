# -*- coding: utf-8 -*-
"""
Análise de Cluster do PDAD Rural - Modelo Expandido (5 Dimensões)
Inclui as variáveis da dimensão D5_Social_Institucional.
Adaptado para ambiente local (Git/Python nativo).
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re
from pathlib import Path
import os
import glob
import unicodedata
import zipfile
import geopandas as gpd
from scipy.stats import kruskal
from sklearn.cluster import KMeans
from sklearn.metrics import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
    pairwise_distances
)
from sklearn.mixture import GaussianMixture
from statsmodels.stats.multitest import multipletests
import matplotlib.patches as mpatches
import matplotlib as mpl
import libpysal
import esda

warnings.filterwarnings('ignore')

# ==============================================================================
# CONFIGURAÇÃO DE CAMINHOS LOCAIS DO REPOSITÓRIO GIT
# ==============================================================================
BASE_DIR = Path(__file__).parent.resolve()
DATA_DIR = BASE_DIR / "dados"

PATH_DOM = DATA_DIR / "pdadr_domicilios.xlsx"
PATH_MOR = DATA_DIR / "pdadr_moradores.xlsx"
ARQUIVO_ZIP_SHP = DATA_DIR / "regioes_administrativas.zip"
PASTA_SHP = DATA_DIR / "shape_df_extracao"

# ==============================================================================
# FUNÇÕES SUPLEMENTARES DE PADRONIZAÇÃO GEOGRÁFICA DO DF
# ==============================================================================
DE_PARA_RAS_TEMP = {
    "sol nascente": "sol nascente por do sol", "por do sol": "sol nascente por do sol",
    "sol nascente por do sol": "sol nascente por do sol", "sol nascente e por do sol": "sol nascente por do sol",
    "nucleo bandeirante": "nucleo bandeirante", "bandeirante": "nucleo bandeirante",
    "n bandeirante": "nucleo bandeirante", "ra viii nucleo bandeirante": "nucleo bandeirante",
    "scia": "scia estrutural", "estrutural": "scia estrutural", "scia estrutural": "scia estrutural",
    "riacho fundo i": "riacho fundo", "riacho fundo": "riacho fundo", "riacho fundo ii": "riacho fundo ii",
    "varjao": "varjao", "fercal": "fercal", "itapoan": "itapoa", "itapoa": "itapoa",
    "sra": "sia", "sia": "sia", "brasilia": "plano piloto", "plano piloto": "plano piloto",
}

def limpar_nome_ra_base(texto):
    if pd.isna(texto): return ""
    s = unicodedata.normalize("NFKD", str(texto)).encode("ASCII", "ignore").decode("utf-8").lower()
    s = re.sub(r"^ra\s+[ivxlcdm0-9]+\s*[-_–—]?\s*", "", s)
    s = re.sub(r"^regiao\s+administrativa\s+(de\s+)?", "", s)
    s = re.sub(r"[/_\-\.]", " ", s)
    s = re.sub(r"\be\b", " ", s)
    s = " ".join(s.split())
    return DE_PARA_RAS_TEMP.get(s, s)

# ==============================================================================
# PIPELINE PRINCIPAL 5D
# ==============================================================================
def main():
    print("🚀 Executando Fase 1: Saneamento, Tratamento de Pesos e Agregação (5D)...\n")

    for p in [PATH_DOM, PATH_MOR, ARQUIVO_ZIP_SHP]:
        if not p.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {p}. Coloque na pasta 'dados/'.")

    df_dom = pd.read_excel(PATH_DOM)
    df_mor = pd.read_excel(PATH_MOR)

    df_dom.columns = df_dom.columns.astype(str).str.strip()
    df_mor.columns = df_mor.columns.astype(str).str.strip()

    ID, RA, SUB, ESTRATO = "domicilio_id", "a01ra", "a01subpopulacao", "a01estratoV2"
    PESO_DOM, PESO_MOR = "peso_dom", "peso_mor"
    NPESS, ESC = "npessoas_qtd_de_pessoas_no_domicilio", "escolaridade"
    CHAVES = [RA, SUB]

    for df in [df_dom, df_mor]:
        for col in [RA, SUB, ESTRATO]:
            df[col] = (df[col].astype("string").str.replace(r"\s+", " ", regex=True)
                       .str.strip().str.replace("Assentamento/ Agrovila", "Assentamento/Agrovila", regex=False))

    df_dom[PESO_DOM] = pd.to_numeric(df_dom[PESO_DOM], errors="coerce")
    df_mor[PESO_MOR] = pd.to_numeric(df_mor[PESO_MOR], errors="coerce")

    N_H = {
        "Agricultura Empresarial - Leste": 5087, "Agricultura Empresarial - Norte": 2968, "Agricultura Empresarial - Oeste": 2613, "Agricultura Empresarial - Sul": 4663,
        "Agricultura Familiar - Central Adjacente I": 595, "Agricultura Familiar - Central Adjacente II": 688, "Agricultura Familiar - Leste": 1483, "Agricultura Familiar - Norte": 6097, "Agricultura Familiar - Oeste": 6998, "Agricultura Familiar - Sul": 1216,
        "Assentamento/Agrovila - Central Adjacente I": 106, "Assentamento/Agrovila - Leste": 1213, "Assentamento/Agrovila - Norte": 318, "Assentamento/Agrovila - Oeste": 465, "Assentamento/Agrovila - Sul": 351,
        "Características Urbanas - Central Adjacente I": 218, "Características Urbanas - Central Adjacente II": 5452, "Características Urbanas - Leste": 9664, "Características Urbanas - Norte": 5556, "Características Urbanas - Oeste": 5780, "Características Urbanas - Sul": 5009,
    }

    obs = df_dom.groupby(ESTRATO).agg(n_h=(ID, "count"), w_orig=(PESO_DOM, "first")).reset_index()
    obs["N_h"] = obs[ESTRATO].map(N_H)
    obs["w_esperado"] = obs["N_h"] / obs["n_h"]
    obs["fator_escala"] = 10 ** np.round(np.log10(obs["w_orig"] / obs["w_esperado"]))
    obs["peso_dom_corrigido"] = obs["w_orig"] / obs["fator_escala"]
    df_dom = df_dom.merge(obs[[ESTRATO, "peso_dom_corrigido"]], on=ESTRATO, how="left")

    def parse_binaria(series):
        s = series.astype("string").str.strip().str.lower()
        res = pd.Series(np.nan, index=series.index, dtype=float)
        res[s.isin(["sim", "1"])] = 1.0
        res[s.isin(["não", "nao", "0"])] = 0.0
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
            xw = np.where(valid, x * w, 0.0)
            w_val = np.where(valid, w, 0.0)
            temp_df = pd.DataFrame({"xw": xw, "w_val": w_val})
            grouped = temp_df.groupby([df[c] for c in group_cols], dropna=False).sum()
            result_df[var] = np.where(grouped["w_val"] > 0, grouped["xw"] / grouped["w_val"], np.nan)
        return result_df

    df_dom["ra_norm"] = df_dom[RA].map(limpar_nome_ra_base)
    MAPA_IVS_DF = {"scia estrutural": 0.72, "sol nascente por do sol": 0.60, "fercal": 0.55, "varjao": 0.53, "itapoa": 0.53, "sudoeste octogonal": 0.09, "aguas claras": 0.10, "cruzeiro": 0.12, "sia": 0.13, "lago sul": 0.14}
    df_dom["ivs_df"] = df_dom["ra_norm"].map(MAPA_IVS_DF).fillna(0.34)

    volumetria_ra = df_dom.groupby("ra_norm")["peso_dom_corrigido"].sum()
    df_dom["qtd_estabelecimentos"] = df_dom["ra_norm"].map(volumetria_ra)
    df_dom["pct_estabelecimentos"] = df_dom["ra_norm"].map(volumetria_ra / volumetria_ra.sum())

    df_dom["renda_domiciliar_pc"] = np.maximum(pd.to_numeric(df_dom.get("renda_domiciliar_per_capita", df_dom.get("renda_domiciliar_calculada", 0) / df_dom[NPESS]), errors="coerce"), 0.0)
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
    df_dom["acesso_energetico_estruturado"] = combine_max(df_dom, ["B16_1_energia_eletrica_rede_geral_ceb", "B16_3_energia_eletrica_proprio_energia_solar", "B16_4_energia_eletrica_outras_fontes_renovaveis"], parse_binaria)
    df_dom["destino_com_coleta_direta"] = combine_max(df_dom, ["B17_1_tratamento_do_lixo_coleta_seletiva_direta", "B17_2_tratamento_do_lixo_coleta_convencional_direta_ou_nao_seletiva"], parse_binaria)
    df_dom["rua_pavimentada"] = parse_binaria(df_dom["B18_1_infraest_urbana_a_rua_de_acesso_principal_asfaltada_pavimentada"])
    df_dom["iluminacao_publica"] = parse_binaria(df_dom["B18_5_infraest_urbana_na_rua_tem_iluminacao"])
    df_dom["drenagem"] = parse_binaria(df_dom["B18_6_infraest_urbana_na_rua_tem_drenagem_de_agua_de_chuva"])
    df_dom["ponto_onibus"] = parse_binaria(df_dom["B20_8_nas_proximidades_existe_ponto_de_onibus"])
    df_dom["internet_domiciliar"] = parse_binaria(df_dom["C05_o_domicilio_ou_alguem_do_domicilio_possuia_acesso_a_internet_no_ultimo_mes"])
    df_dom["internet_banda_larga_fixa"] = parse_binaria(df_dom["C05_1_1_acesso_por_internet_banda_larga_fixa"])
    df_dom["internet_movel"] = parse_binaria(df_dom["C05_1_2_acesso_por_internet_de_rede_de_celular"])
    df_dom["possui_computador"] = combine_max(df_dom, ["C07_12_qtd_de_microcomputadores_desktop", "C07_13_qtd_de_notebooks_laptops"])

    # VARIÁVEIS INSTITUCIONAIS - D5
    col_inc = "B22_1_a_sua_propriedade_ja_recebeu_incentivos_a_producao"
    df_dom["incentivo_producao"] = parse_binaria(df_dom[col_inc]) if col_inc in df_dom.columns else np.nan

    MAPA_ESC = {"sem instrução": 1, "alfabetização de jovens e adultos": 2, "ensino fundamental incompleto": 3, "ensino fundamental completo": 4, "ensino médio incompleto": 5, "ensino médio completo": 6, "superior incompleto": 7, "superior completo": 8, "especialização de nível superior": 9, "mestrado": 10, "doutorado": 11}
    df_mor["idade_num"] = pd.to_numeric(df_mor.get("E02_idade__Idade_informada", df_mor.get("E02_idade")), errors="coerce")
    df_mor["escolaridade_num"] = pd.to_numeric(df_mor[ESC], errors="coerce").fillna(df_mor[ESC].astype("string").str.strip().str.lower().map(MAPA_ESC))
    df_mor["escolaridade_adulto"] = np.where(df_mor["idade_num"] >= 25, df_mor["escolaridade_num"], np.nan)
    df_mor["trabalhou_30d"] = parse_binaria(df_mor["I05__Trabalhou_nos_últimos_30_dias"])
    df_mor["tem_ctps"] = parse_binaria(df_mor["I17__Possui_CTPS_assinada"])
    df_mor["contribui_previdencia"] = parse_binaria(df_mor["I18__Contribui_para_alguma_Previdência_Social_Pública"])
    df_mor["internet_uso_3m"] = parse_binaria(df_mor["F04__Acessou_a_Internet_nos_últimos_3_meses"])
    df_mor["internet_computador"] = parse_binaria(df_mor["F05_1__Utilizou_microcomputador_para_acessar_a_internet"])

    # VARIÁVEIS INSTITUCIONAIS MORADOR - D5
    col_ben = "I24_2__Recebe_ou_já_recebeu_benefícios_de_programas_sociais"
    col_val = "I22_1__Valor__R___recebido_no_mês_passado_de_benefícios_sociais__Bolsa_família__BPC_LOAS__Bolsa_de_estudo__valor"
    col_aux = "F06_7__Pedido_de_auxílio_benefícios_governamentais"
    df_mor["recebe_beneficios_sociais"] = parse_binaria(df_mor[col_ben]) if col_ben in df_mor.columns else np.nan
    df_mor["valor_beneficios_sociais"] = pd.to_numeric(df_mor[col_val], errors="coerce") if col_val in df_mor.columns else np.nan
    df_mor["pedido_auxilio_governamental"] = parse_binaria(df_mor[col_aux]) if col_aux in df_mor.columns else np.nan

    df_dom["tipologia_original"] = df_dom[SUB].astype(str)

    VAR_MOR = ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia", "recebe_beneficios_sociais", "valor_beneficios_sociais", "pedido_auxilio_governamental"]
    VAR_DOM = ["renda_domiciliar_pc", "possui_automovel", "possui_motocicleta", "possui_bicicleta", "possui_geladeira", "possui_freezer", "possui_lavadora", "possui_secadora", "possui_tv", "possui_ar_condicionado", "possui_microondas", "comodos_por_pessoa", "banheiros_por_pessoa", "agua_rede", "agua_poco_artesiano", "agua_poco_cisterna", "captacao_chuva", "caixa_agua", "solucao_sanitaria_estruturada", "acesso_energetico_estruturado", "destino_com_coleta_direta", "rua_pavimentada", "iluminacao_publica", "drenagem", "ponto_onibus", "internet_domiciliar", "internet_banda_larga_fixa", "internet_movel", "possui_computador", "ivs_df", "qtd_estabelecimentos", "pct_estabelecimentos", "incentivo_producao"]

    matriz_mor_A = media_ponderada_vetorizada(df_mor, CHAVES, VAR_MOR, PESO_MOR)
    matriz_dom_A = media_ponderada_vetorizada(df_dom, CHAVES, VAR_DOM, "peso_dom_corrigido")
    matriz_dom_A["log_renda_domiciliar_pc"] = np.where(matriz_dom_A["renda_domiciliar_pc"] >= 0, np.log1p(matriz_dom_A["renda_domiciliar_pc"]), np.nan)

    counts_A = df_dom.groupby(CHAVES).agg(amostra_n=(ID, "count"), domicilios_expandidos=("peso_dom_corrigido", "sum")).reset_index()
    df_matriz_A = counts_A.merge(matriz_dom_A, on=CHAVES, how="left").merge(matriz_mor_A, on=CHAVES, how="left")
    df_matriz_A = df_matriz_A[df_matriz_A["amostra_n"] >= 5].reset_index(drop=True)

    sub_dominante = df_dom.groupby(CHAVES)[SUB].agg(lambda x: x.mode()[0] if not x.mode().empty else "Desconhecido").reset_index(name="tipologia_original")
    df_matriz_A = df_matriz_A.merge(sub_dominante, on=CHAVES, how="left")

    print(f"\n✅ FASE 1 CONCLUÍDA: Matriz A de Referência gerada com {len(df_matriz_A)} segmentos territoriais válidos.")

    print("\n🚀 Executando Etapa 4 — PCA por Dimensão (5D) (Monte Carlo Corrigido)...\n")
    DIMENSOES_5D = {
        "D1_Financeiro_Ativos": ["log_renda_domiciliar_pc", "possui_automovel", "possui_motocicleta", "possui_geladeira", "possui_freezer", "possui_lavadora", "possui_secadora", "possui_tv", "possui_ar_condicionado", "possui_microondas"],
        "D2_Humano_Trabalho": ["escolaridade_adulto", "trabalhou_30d", "tem_ctps", "contribui_previdencia", "ivs_df"],
        "D3_Fisico_Habitacional_Territorial": ["comodos_por_pessoa", "banheiros_por_pessoa", "agua_rede", "agua_poco_artesiano", "agua_poco_cisterna", "captacao_chuva", "caixa_agua", "solucao_sanitaria_estruturada", "acesso_energetico_estruturado", "destino_com_coleta_direta", "rua_pavimentada", "iluminacao_publica", "drenagem", "ponto_onibus", "qtd_estabelecimentos", "pct_estabelecimentos"],
        "D4_Tecnologico_Conectividade": ["internet_domiciliar", "internet_banda_larga_fixa", "internet_movel", "possui_computador", "internet_uso_3m", "internet_computador"],
        "D5_Social_Institucional": ["recebe_beneficios_sociais", "valor_beneficios_sociais", "pedido_auxilio_governamental", "incentivo_producao"]
    }

    def monte_carlo_parallel_analysis_corrigido(X_std, n_iterations=1000, random_state=42):
        np.random.seed(random_state)
        n_samples, n_features = X_std.shape
        simulated_eigenvalues = np.zeros((n_iterations, n_features))
        for i in range(n_iterations):
            X_sim = np.random.normal(size=(n_samples, n_features))
            simulated_eigenvalues[i, :] = np.linalg.eigvalsh(np.corrcoef(X_sim, rowvar=False))[::-1]
        return np.percentile(simulated_eigenvalues, 95, axis=0)

    df_componentes_mc = pd.DataFrame(index=df_matriz_A.index)
    df_A_mc = df_matriz_A.copy()
    df_A_mc["ra_norm"] = df_A_mc[RA].map(limpar_nome_ra_base)
    resumo_pca_dual = []

    for dim, vars_list in DIMENSOES_5D.items():
        valid_cols = [v for v in vars_list if v in df_matriz_A.columns and df_matriz_A[v].notna().sum() >= 3 and df_matriz_A[v].dropna().nunique() > 1]
        if not valid_cols: continue

        X_A_std = StandardScaler().fit_transform(df_matriz_A[valid_cols].apply(pd.to_numeric, errors="coerce").fillna(0))
        eigen_obs = np.linalg.eigvalsh(np.corrcoef(X_A_std, rowvar=False))[::-1]
        eigen_sim_p95 = monte_carlo_parallel_analysis_corrigido(X_A_std)

        # Retenção forçada para preservar o sinal governamental
        k_mc = int(np.sum(eigen_obs > eigen_sim_p95))
        if k_mc == 0: k_mc = 1 
        
        k_kai = max(1, int(np.sum(eigen_obs >= 1.0)))

        scores_mc = PCA(n_components=k_mc).fit_transform(X_A_std)
        for j in range(k_mc):
            col_name = f"{dim}_MC_PC{j+1}"
            df_A_mc[col_name] = scores_mc[:, j]
            df_componentes_mc[col_name] = scores_mc[:, j]

        pca_full = PCA().fit(X_A_std)
        for j in range(len(eigen_obs)):
            resumo_pca_dual.append({"dimensao": dim, "componente": j + 1, "eigenvalue_obs": round(eigen_obs[j], 4), "eigenvalue_mc_p95": round(eigen_sim_p95[j], 4), "retido_mc": j < k_mc, "retido_kaiser": j < k_kai, "variancia_explicada_pct": round(pca_full.explained_variance_ratio_[j] * 100, 2), "variancia_acumulada_pct": round(np.cumsum(pca_full.explained_variance_ratio_ * 100)[j], 2)})

    print("\n" + "=" * 110)
    print("📊 COMPARAÇÃO DUAL DE RETENÇÃO (MONTE CARLO CORRIGIDO VS. KAISER)")
    print("=" * 110)
    print(pd.DataFrame(resumo_pca_dual).to_string(index=False))

    print("\n🚀 Executando Etapa 5 — Benchmarking Aberto de K (2 a 10) com Espaço Purificado...\n")
    X_pca = df_componentes_mc.apply(pd.to_numeric).to_numpy()
    N, p = X_pca.shape

    def compute_gap_statistic(X, k, n_refs=10):
        disp_real = np.sum(np.min(pairwise_distances(X, KMeans(k, n_init=10, random_state=42).fit(X).cluster_centers_)**2, axis=1))
        ref_disps = np.zeros(n_refs)
        for i in range(n_refs):
            random_X = np.random.uniform(np.min(X, axis=0), np.max(X, axis=0), size=X.shape)
            km_sim = KMeans(k, n_init=10, random_state=42).fit(random_X)
            ref_disps[i] = np.sum(np.min(pairwise_distances(random_X, km_sim.cluster_centers_)**2, axis=1))
        return np.mean(np.log(ref_disps)) - np.log(disp_real)

    resultados_k = []
    for k in range(2, min(11, N)):
        km_sim = KMeans(n_clusters=k, n_init=20, random_state=42).fit(X_pca)
        try: bic_val = GaussianMixture(n_components=k, covariance_type='diag', random_state=42).fit(X_pca).bic(X_pca)
        except: bic_val = np.nan
        resultados_k.append({"K": k, "Silhouette": round(silhouette_score(X_pca, km_sim.labels_), 4), "Gap_Statistic": round(compute_gap_statistic(X_pca, k), 4), "BIC": round(bic_val, 2) if not np.isnan(bic_val) else np.nan, "Calinski_Harabasz": round(calinski_harabasz_score(X_pca, km_sim.labels_), 2), "Davies_Bouldin": round(davies_bouldin_score(X_pca, km_sim.labels_), 4)})

    print("\n" + "=" * 110)
    print("📌 TABELA COMPLETA DE BENCHMARKING DE K (2 A 10)")
    print("=" * 110)
    print(pd.DataFrame(resultados_k).to_string(index=False))

    print("\n🧪 Executando Validação de Ciências Sociais: Teste de Referência Nula (Permutação)...\n")
    def executar_teste_referencia_nula(X_real, k_escolhido, n_permutacoes=1000, random_state=42):
        np.random.seed(random_state)
        n_samples, n_features = X_real.shape
        silhueta_real = silhouette_score(X_real, KMeans(n_clusters=k_escolhido, n_init=20, random_state=random_state).fit(X_real).labels_)
        silhuetas_nulas = []
        for i in range(n_permutacoes):
            X_nulo = np.copy(X_real)
            for col in range(n_features): np.random.shuffle(X_nulo[:, col])
            silhuetas_nulas.append(silhouette_score(X_nulo, KMeans(n_clusters=k_escolhido, n_init=10, random_state=i).fit(X_nulo).labels_))
        silhuetas_nulas = np.array(silhuetas_nulas)
        return silhueta_real, np.mean(silhuetas_nulas), np.sum(silhuetas_nulas >= silhueta_real) / n_permutacoes

    for k_test in [6, 8]:
        sil_real, med_nula, p_val = executar_teste_referencia_nula(X_pca, k_test, n_permutacoes=1000)
        print(f"📊 RESULTADOS PARA K = {k_test} (p-valor: {p_val:.4f} | Silhueta Observada: {sil_real:.4f})")

    print("\n🚀 Executando Etapa 6 — K-Means Final (K=6) e Cruzamento Populacional...\n")
    kmeans_final = KMeans(n_clusters=6, n_init=50, random_state=42).fit(X_pca)
    df_A_mc["Cluster_K6"] = kmeans_final.labels_
    col_peso_dom = "domicilios_expandidos"

    tabela_cruzada_tipologia = df_A_mc.pivot_table(index="Cluster_K6", columns="tipologia_original", values=col_peso_dom, aggfunc="sum", fill_value=0)
    print("📌 Domicílios Expandidos por Cluster (K=6) e Tipologia Original:")
    print(tabela_cruzada_tipologia.to_string())

    print("\n🚀 Executando Etapa 7 — Análise Espacial Local (LISA)...\n")
    df_ra_espacial = df_A_mc.groupby('ra_norm').agg({"D4_Tecnologico_Conectividade_MC_PC1": "mean", col_peso_dom: "sum"}).reset_index()
    VIZINHOS_DF = {"planaltina": ["sobradinho", "sobradinho ii", "paranoa", "fercal"], "sobradinho": ["planaltina", "sobradinho ii", "fercal", "paranoa", "varjao", "itapoa", "lago norte"], "sobradinho ii": ["sobradinho", "fercal", "brazlandia"], "fercal": ["sobradinho", "sobradinho ii", "planaltina", "brazlandia"], "paranoa": ["planaltina", "sobradinho", "itapoa", "sao sebastiao", "jardim botanico", "lago norte"], "itapoa": ["paranoa", "sobradinho"], "sao sebastiao": ["paranoa", "jardim botanico", "santa maria"], "jardim botanico": ["sao sebastiao", "paranoa", "lago sul", "santa maria"], "brazlandia": ["sobradinho ii", "fercal", "ceilandia", "taguatinga"], "ceilandia": ["brazlandia", "sol nascente por do sol", "samambaia", "taguatinga"], "sol nascente por do sol": ["ceilandia", "samambaia", "taguatinga"], "samambaia": ["ceilandia", "sol nascente por do sol", "taguatinga", "riacho fundo ii", "recanto das emas"], "recanto das emas": ["samambaia", "riacho fundo ii", "gama", "santa maria"], "gama": ["recanto das emas", "riacho fundo ii", "park way", "santa maria"], "santa maria": ["gama", "park way", "jardim botanico", "sao sebastiao", "recanto das emas"], "riacho fundo ii": ["samambaia", "recanto das emas", "gama", "riacho fundo", "park way"], "riacho fundo": ["riacho fundo ii", "park way", "nucleo bandeirante", "guara", "taguatinga"], "park way": ["gama", "santa maria", "jardim botanico", "lago sul", "nucleo bandeirante", "riacho fundo", "riacho fundo ii", "guara"], "nucleo bandeirante": ["riacho fundo", "park way", "guara", "candangolandia"], "guara": ["riacho fundo", "park way", "nucleo bandeirante", "candangolandia", "sia", "taguatinga", "aguas claras", "cruzeiro", "plano piloto"], "taguatinga": ["ceilandia", "sol nascente por do sol", "samambaia", "riacho fundo", "guara", "aguas claras", "vicente pires", "brazlandia"], "aguas claras": ["taguatinga", "guara", "vicente pires"], "vicente pires": ["taguatinga", "aguas claras"], "lago sul": ["park way", "jardim botanico", "plano piloto", "cruzeiro", "lago norte", "candangolandia"], "lago norte": ["sobradinho", "paranoa", "plano piloto", "lago sul"], "plano piloto": ["guara", "cruzeiro", "lago sul", "lago norte", "sobradinho", "varjao", "sia"], "cruzeiro": ["guara", "plano piloto", "lago sul", "sra"], "varjao": ["sobradinho", "plano piloto"], "sia": ["guara", "plano piloto"], "candangolandia": ["nucleo bandeirante", "guara", "lago sul"]}
    
    RAs_presentes = df_ra_espacial['ra_norm'].tolist()
    adj_dict = {ra: [v for v in VIZINHOS_DF.get(ra, []) if v in RAs_presentes] or [RAs_presentes[0]] for ra in RAs_presentes}
    w = libpysal.weights.W(adj_dict)
    w.transform = 'R'
    lm = esda.moran.Moran_Local(df_ra_espacial["D4_Tecnologico_Conectividade_MC_PC1"].fillna(0).to_numpy(), w, permutations=999, seed=42)
    rejeita, _, _, _ = multipletests(lm.p_sim, alpha=0.05, method='fdr_bh')
    df_ra_espacial["LISA_Significativo"] = rejeita
    df_ra_espacial["LISA_Q"] = lm.q
    
    print(f"  - Hotspots (Q3 + Sig): {np.sum((lm.q == 3) & rejeita)}")
    print(f"  - Coldspots (Q1 + Sig): {np.sum((lm.q == 1) & rejeita)}")

    print("\n🚀 Detalhamento Espacial (LISA cruzado com Tipologia)...")
    df_lisa_sig = df_ra_espacial[df_ra_espacial["LISA_Significativo"] == True][["ra_norm", "LISA_Q"]]
    if not df_lisa_sig.empty:
        mapa_quadrantes = {1: "Q1 (Alto-Alto)", 2: "Q2 (Baixo-Alto)", 3: "Q3 (Baixo-Baixo)", 4: "Q4 (Alto-Baixo)"}
        df_lisa_detalhe = df_A_mc.merge(df_lisa_sig, on="ra_norm", how="inner")
        df_lisa_detalhe["Tipo_LISA"] = df_lisa_detalhe["LISA_Q"].map(mapa_quadrantes)
        print(df_lisa_detalhe.pivot_table(index="Tipo_LISA", columns="tipologia_original", values=col_peso_dom, aggfunc="sum", fill_value=0).to_string())

    print("\n🚀 Gerando Gráficos e Mapas (Feche a janela de cada gráfico para continuar)...\n")
    sns.set_theme(style="whitegrid")
    tabela_cruzada_tipologia.plot(kind='bar', stacked=True, figsize=(12, 7), colormap='viridis', edgecolor='black')
    plt.title('Composição dos Clusters (5D) por Tipologia Original', fontsize=16, fontweight='bold', pad=15)
    plt.tight_layout(); plt.show()

    if not PASTA_SHP.exists():
        with zipfile.ZipFile(ARQUIVO_ZIP_SHP, 'r') as zip_ref: zip_ref.extractall(PASTA_SHP)
    
    gdf_df = gpd.read_file(list(PASTA_SHP.glob("**/*.shp"))[0])
    col_nome_ra = next((col for col in gdf_df.columns if col.lower() in ['ra', 'nome', 'ra_nome', 'nome_ra', 'regiao', 'nm_ra']), gdf_df.columns[0])
    gdf_df['ra_norm'] = gdf_df[col_nome_ra].apply(limpar_nome_ra_base).replace({'nacleo bandeirante': 'nucleo bandeirante'})
    
    dom_cluster_ra = df_A_mc.groupby(['ra_norm', 'Cluster_K6'])[col_peso_dom].sum().reset_index().sort_values(col_peso_dom, ascending=False).drop_duplicates('ra_norm')
    gdf_mapa = gdf_df.merge(dom_cluster_ra[['ra_norm', 'Cluster_K6']], on='ra_norm', how='left')

    fig, ax = plt.subplots(1, 1, figsize=(14, 10))
    gdf_df.plot(ax=ax, color='lightgrey', edgecolor='white', linewidth=0.5)
    gdf_mapa.dropna(subset=['Cluster_K6']).plot(column='Cluster_K6', ax=ax, cmap='viridis', edgecolor='black', linewidth=0.8, legend=False)
    
    clusters_unicos = sorted(dom_cluster_ra['Cluster_K6'].dropna().unique())
    cmap_viridis = mpl.colormaps['viridis'].resampled(max(len(clusters_unicos), 2))
    ax.legend(handles=[mpatches.Patch(color=cmap_viridis(idx / (len(clusters_unicos) - 1)) if len(clusters_unicos) > 1 else cmap_viridis(0), label=f'Cluster {int(c)}') for idx, c in enumerate(clusters_unicos)], title='Cluster Predominante', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Tipologia Territorial 5D: Cluster Dominante por RA', fontsize=18, fontweight='bold', pad=20)
    ax.set_axis_off(); plt.tight_layout(); plt.show()

if __name__ == "__main__":
    main()
