# -*- coding: utf-8 -*-
"""PDAD Rural 2022 — Saneamento dos dados (camada silver -> camada gold).

Versão consolidada e revisada do notebook "Análise de Qualidade dos Dados e Saneamento".
Entrada : Dicionario_aplicado_PDADR_2022.xlsx (abas 'PDADR_2022_Domicilios' e 'PDADR 2022 Moradores').
Saída   : pdadr_domicilios.xlsx e pdadr_moradores.xlsx (colunas originais + pares *_classe / *_valor
          + renda_individual_total, renda_domiciliar_calculada, renda_domiciliar_per_capita).

O que mudou em relação ao notebook anterior (resumo; detalhes na seção 0 e nos comentários):
  1. Blocos duplicados removidos: havia duas definições de PATTERN_FAT_FINGER, codigos_sujos,
     aplicar_dupla_dimensao e duas execuções da dupla dimensão e da renda individual; a segunda
     sobrescrevia a primeira. Agora há UMA função e UMA tabela de regras por domínio, que reproduz o
     comportamento efetivo anterior (verificado contra os xlsx publicados).
  2. Código 66666 ("Não tem outros rendimentos") em I22_1/I22_2, digitado como 666, 6666, 666666 ou
     misturado com 9: agora vale R$ 0,00. Antes virava erro e anulava a renda do morador inteiro.
  3. Quase-códigos (77774 em I20; 36666 e 66663 em I22): viram NaN (erro de digitação de código).
  4. Correções manuais por domicílio (D2089, D2150) removidas: já eram cobertas pela trava geral de
     R$ 100 mil por fonte e a regra de D2150 anulava I21 de TODOS os moradores do domicílio.
  5. A "Validade 100%" (val_limpo_perc = 100.00, digitada à mão) e a "Completude" (tabelas com número de
     colunas diferente) foram substituídas por contagens por classe e por verificações que interrompem
     a execução (assert) se restar código administrativo, valor >= R$ 100 mil ou quase-código.
  6. Exportação corrigida (o bloco `with` estava comentado e o corpo não: IndentationError).
  7. Nomes de RA em E23_1/I08 unificados com a grafia de a01ra (antes: 'Recanto Das Emas' vs 'Recanto das Emas').
  8. Adicionadas as colunas de transparência n_moradores_renda_ausente e renda_dom_parcial.
"""

# %% [0] Bibliotecas e configuração ------------------------------------------------------------------
import os
import re

import numpy as np
import pandas as pd

pd.set_option("display.max_columns", None)

try:                                     # no Colab/Jupyter, display() mostra tabelas formatadas
    display                              # noqa: B018
except NameError:                        # fora do notebook (terminal), cai para print
    display = print

# Colab: monta o Google Drive (necessário para ler/gravar em /content/drive/MyDrive)
if "google.colab" in __import__("sys").modules and not os.path.isdir("/content/drive/MyDrive"):
    from google.colab import drive
    drive.mount("/content/drive")

PATH_SILVER = os.environ.get("PDADR_SILVER", "/content/drive/MyDrive/Dicionario_aplicado_PDADR_2022.xlsx")
SHEET_DOM = "PDADR_2022_Domicilios"
SHEET_MOR = "PDADR 2022 Moradores"
OUT_DIR = os.environ.get("PDADR_OUT_DIR", "/content/drive/MyDrive")
EXPORTAR = os.environ.get("PDADR_EXPORTAR", "1") == "1"


# %% [1] Carga e verificações iniciais ----------------------------------------------------------------
df_domicilios = pd.read_excel(PATH_SILVER, sheet_name=SHEET_DOM)
df_moradores = pd.read_excel(PATH_SILVER, sheet_name=SHEET_MOR)
print(f"Domicílios: {df_domicilios.shape} | Moradores: {df_moradores.shape}")

print(f"Linhas duplicadas — Domicílios: {df_domicilios.duplicated().sum()} | Moradores: {df_moradores.duplicated().sum()}")
assert df_domicilios["domicilio_id"].is_unique, "domicilio_id duplicado na tabela de domicílios"


def top_nulos(df, nome, n=10):
    nul = df.isnull().sum()
    t = pd.DataFrame({"Faltantes": nul, "Percentual (%)": nul / len(df) * 100})
    t = t[t["Faltantes"] > 0].sort_values("Percentual (%)", ascending=False).head(n)
    print(f"\nTop {n} colunas com mais nulos — {nome}:")
    display(t.round(2))


top_nulos(df_domicilios, "Domicílios")
top_nulos(df_moradores, "Moradores")


# %% [2] Padronização sintática -----------------------------------------------------------------------
df_domicilios.columns = [c.replace(".", "_") for c in df_domicilios.columns]
df_moradores.columns = [c.replace(".", "_") for c in df_moradores.columns]

df_domicilios = df_domicilios.rename(columns={
    "C010_1": "C10_1_local_de_compra_no_ultimo_mes_de_produtos_de_alimentacao_higiene_limpeza",
    "C010_2": "C10_2_local_de_compra_nos_ultimos_12_meses_de_eletrodomesticos",
    "C010_3": "C10_3_local_de_compra_nos_ultimos_12_meses_de_materiais_de_construcao_manutencao",
    "C010_4": "C10_4_local_de_compra_nos_ultimos_12_meses_de_servicos_exceto_limpeza_domestica",
})

# Regiões Administrativas (grafia única, usada também em E23_1 e I08)
MAP_RA = {
    1: "Plano Piloto", 2: "Gama", 3: "Taguatinga", 4: "Brazlândia", 5: "Sobradinho",
    6: "Planaltina", 7: "Paranoá", 8: "Núcleo Bandeirante", 9: "Ceilândia",
    10: "Guará", 11: "Cruzeiro", 12: "Samambaia", 13: "Santa Maria",
    14: "São Sebastião", 15: "Recanto das Emas", 16: "Lago Sul", 17: "Riacho Fundo",
    18: "Lago Norte", 19: "Candangolândia", 20: "Águas Claras", 21: "Riacho Fundo II",
    22: "Sudoeste e Octogonal", 23: "Varjão", 24: "Park Way", 25: "SCIA",
    26: "Sobradinho II", 27: "Jardim Botânico", 28: "Itapoã", 29: "SIA",
    30: "Vicente Pires", 31: "Fercal", 32: "Sol Nascente/Pôr do Sol", 33: "Arniqueira",
}
MAP_RA_FIX = {**MAP_RA, **{str(k): v for k, v in MAP_RA.items()}, **{f"{float(k):.1f}": v for k, v in MAP_RA.items()}}
VARIANTES_RA = {
    "Sol Nascente / Pôr do Sol": "Sol Nascente/Pôr do Sol",
    "Sol Nascente / Por do Sol": "Sol Nascente/Pôr do Sol",
    "Sol Nascente/Por do Sol": "Sol Nascente/Pôr do Sol",
    "Recanto Das Emas": "Recanto das Emas",
}
for _df in (df_domicilios, df_moradores):
    if "a01ra" in _df.columns:
        s = _df["a01ra"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
        s = s.str.replace(r"\s*/\s*", "/", regex=True)
        _df["a01ra"] = s.replace(MAP_RA_FIX).replace(VARIANTES_RA)

# E23_1 (Anexo 9) e I08 (Anexo 12) — apenas moradores
MAP_E23_1 = {0: "Na Região Administrativa", 88888: "Não sabe", 99999: "Não se aplica", **MAP_RA}
MAP_I08 = {
    **MAP_RA,
    0: "Na Região Administrativa", 34: "Vários locais",
    41: "Águas Lindas de Goiás", 42: "Alexânia", 43: "Cidade Ocidental", 44: "Cristalina",
    45: "Cocalzinho de Goiás", 46: "Formosa", 47: "Luziânia", 48: "Novo Gama", 49: "Padre Bernardo",
    50: "Planaltina", 51: "Santo Antônio do Descoberto", 52: "Valparaíso de Goiás", 53: "Anápolis",
    54: "Goiânia", 55: "Outros municípios de Goiás", 56: "No domicílio", 57: "Outros locais",
    88888: "Não sabe", 99999: "Não se aplica",
}
df_moradores["E23_1"] = df_moradores["E23_1"].replace(MAP_E23_1)
df_moradores["I08"] = df_moradores["I08"].replace(MAP_I08)

# Posse de animais (C11): 'Não tem' / 'Não se aplica' -> classe; quantidade não numérica -> 0
def tratar_posse_animais(df, col_base, nome):
    df[f"C11_{nome}_posse"] = "Tem"
    df.loc[df[col_base] == "Não tem", f"C11_{nome}_posse"] = "Não tem"
    df.loc[df[col_base] == "Não se aplica", f"C11_{nome}_posse"] = "Não se aplica"
    df[f"C11_{nome}_valor"] = pd.to_numeric(df[col_base], errors="coerce").fillna(0).astype(int)
    return df


for col, nome in [("C11_1_qtd_de_gatos", "gato"), ("C11_2_qtd_de_cachorros", "cachorro"), ("C11_3_qtd_de_aves", "ave"),
                  ("C11_4_qtd_de_peixes", "peixe"), ("C11_5_qtd_de_outros_tipos_de_animais", "outros")]:
    df_domicilios = tratar_posse_animais(df_domicilios, col, nome)
print("Padronização sintática concluída.")


# %% [3] Regras de códigos administrativos e função única de "Dupla Dimensão" -------------------------
# Dicionário IPEDF: I04_1_1, I04_2_1, I20, I21, I22_1, I22_2 usam 77777 = Recusa, 88888 = Não sabe,
# 99999 = Não se aplica; I22_1 e I22_2 têm ainda 66666 = "Não tem outros rendimentos".
# Nos microdados os códigos aparecem com dígitos a mais/menos (8888, 666666, 7777777...): daí a regex.
FAT_FINGER = re.compile(r"^[6789]{4,}$")          # repetições de 6/7/8/9 com 4+ dígitos -> erro/recusa
NAO_TEM_I22 = re.compile(r"^[69]{3,}$")           # 66666 e variantes (666, 6666, 666666, 66699...) em I22 -> 0

CODIGOS_SUJOS = {
    "Recusa", "Não sabe", "Não Sabe",
    "77777", "88888", "8888", "88885", "88880",
    "666666", "777777", "888888", "999999", "7777777", "8888888", "9999999", "7777778", "8888889",
}
# Quase-códigos: um dígito trocado de 77777 / 66666 (confirmados na inspeção dos microdados).
QUASE_CODIGOS = {"77774", "36666", "66663"}

SEM_RENDA = {
    "Não se aplica", "99999", "99999.0",
    "Tem trabalho principal, mas não recebeu no mês anterior",
    "Tem outros trabalhos, mas não recebeu no mês anterior",
    "Não recebeu/Não tem outros rendimentos", "Não tem outros rendimentos", "Não tem", "Não recebeu",
}
SEM_DESPESA = {"Não se aplica", "Não paga", "Não tem", "99999", "99999.0"}
SEM_BENS = {"Não se aplica", "Não tem", "0"}

CLASSE_ZERO, CLASSE_ERRO, CLASSE_OK = "Isento / Ausente", "Código Administrativo / Erro", "Valor Informado"


def classificar(val, zero, nulos, zero_regex=None):
    """Devolve (classe, valor). Prioridade: zero estrutural > código/erro > número > texto."""
    if pd.isna(val):
        return np.nan, np.nan
    s = str(val).strip()
    limpo = s.split(".")[0]
    if s in zero or limpo in zero or (zero_regex is not None and zero_regex.match(limpo)):
        return CLASSE_ZERO, 0.0
    if s in nulos or limpo in nulos or FAT_FINGER.match(limpo):
        return CLASSE_ERRO, np.nan
    try:
        return CLASSE_OK, float(s)
    except ValueError:
        return s, np.nan            # texto livre (ex.: 'Não se aplica' fora do zero estrutural): fica na classe


def dupla_dimensao(df, coluna, zero=(), nulos=CODIGOS_SUJOS, zero_regex=None):
    if coluna not in df.columns:
        print(f"Aviso: coluna '{coluna}' não encontrada; ignorada.")
        return df
    res = [classificar(v, set(zero), set(nulos), zero_regex) for v in df[coluna]]
    df[f"{coluna}_classe"] = [r[0] for r in res]
    df[f"{coluna}_valor"] = [r[1] for r in res]
    return df


# %% [4] Moradores ------------------------------------------------------------------------------------
COLUNAS_RENDA = [
    "I04_1_1__Valor__R___da_aposentadoria_recebida_no_mês_passado",
    "I04_2_1__Valor__R___da_pensão_recebida_no_mês_passado",
    "I20__Valor__R___recebido_no_mês_passado_pelo_Trabalho_Principal",
    "I21__Valor__R___recebido_no_mês_passado_pelos_outros_trabalhos",
    "I22_1__Valor__R___recebido_no_mês_passado_de_benefícios_sociais__Bolsa_família__BPC_LOAS__Bolsa_de_estudo_",
    "I22_2__Valor__R___recebido_no_mês_passado_de_outros_rendimentos_Brutos__Doações__Aluguéis__Aplicações_Financeiras_",
]
for col in COLUNAS_RENDA:
    eh_i22 = col.startswith(("I22_1", "I22_2"))
    df_moradores = dupla_dimensao(
        df_moradores, col, zero=SEM_RENDA, nulos=CODIGOS_SUJOS | QUASE_CODIGOS,
        zero_regex=NAO_TEM_I22 if eh_i22 else None,
    )

for col in ["I07__Qtd__de_horas_semanais_que_trabalha_no_trabalho_principal",
            "I23__Qtd__de_horas_dedica_semanalmente_aos_afazeres_domésticos"]:
    df_moradores = dupla_dimensao(df_moradores, col)                        # 'Não se aplica' fica NaN (condicional)

for col in ["F01__Qtd__de_aparelho_s__celular_es__para_uso_pessoal",
            "F02_1__Qtd__de_linhas_pré_pagas", "F02_2__Qtd__de_linhas_pós_pagas"]:
    df_moradores = dupla_dimensao(df_moradores, col, zero={"Não tem"})

for col in ["E19__Tempo_de_moradia_no_DF", "E20__Tempo_de_moradia_em_área_rural"]:
    df_moradores = dupla_dimensao(df_moradores, col, zero={"Menos de 1 ano"})

for col in ["J02_1__Qtd__de_filhos_homens_nascidos_vivos", "J02_2__Qtd__de_filhas_mulheres_nascidas_vivas",
            "J03_1__Qtd__de_filhos_homens_moram_no_domicílio", "J03_2__Filhas_mulheres_moram_no_domicílio",
            "J04_2__Idade_estimada_de_nascimento_do_último_a__filho_a__nascido_a__vivo_a_"]:
    df_moradores = dupla_dimensao(df_moradores, col)                        # 'Não se aplica' (homens/crianças) fica NaN

# E05_1: vazamento de chave do sistema (ex.: '15306-3') no lugar de 'Sim'
COL_E05_1 = "E05_1__A_responsabilidade_pelo_domicílio_é_compartilhada_com_mais_algum_morador_do_domicílio"
if COL_E05_1 in df_moradores.columns:
    validos = {"Não", "Não se aplica", "Não sabe", "Sim"}
    df_moradores[COL_E05_1] = df_moradores[COL_E05_1].apply(
        lambda x: "Sim" if pd.notna(x) and str(x).strip() not in validos else x)


# %% [5] Domicílios -----------------------------------------------------------------------------------
for col in ["B05_valor_mensal_pago_da_prestacao_no_ultimo_mes", "B06_valor_mensal_do_aluguel_pago_no_ultimo_mes"]:
    df_domicilios = dupla_dimensao(df_domicilios, col, zero=SEM_DESPESA)

for col in ["C01_1_quantidade_de_automoveis", "C02_1_quantidade_de_motocicletas", "C03_1_quantidade_de_bicicletas",
            "C04_1_quantidade_de_carrocas", "C04_2_quantidade_de_cavalos", "C08_1_qtd_de_mensalistas",
            "C09_1_qtd_de_diaristas"] + [c for c in df_domicilios.columns if c.startswith("C07_")
                                         and not c.endswith(("_classe", "_valor"))]:
    df_domicilios = dupla_dimensao(df_domicilios, col, zero=SEM_BENS)
print("Dupla Dimensão aplicada em moradores e domicílios.")


# %% [6] Correções de escala e anomalias extremas na renda ---------------------------------------------
COLS_VALOR = [f"{c}_valor" for c in COLUNAS_RENDA]
COLS_CLASSE = [f"{c}_classe" for c in COLUNAS_RENDA]

# 6.1 Erro de escala do salário mínimo de 2022 (R$ 1.212,00) digitado sem separador decimal
n_escala = 0
for col in COLS_VALOR:
    for errado, certo in ((121200.0, 1212.0), (242400.0, 2424.0)):
        mk = df_moradores[col] == errado
        n_escala += int(mk.sum())
        df_moradores.loc[mk, col] = certo
print(f"Valores de escala corrigidos (121200->1212; 242400->2424): {n_escala}")

# 6.2 Trava geral: fonte individual >= R$ 100.000 é tratada como erro de digitação (NaN).
#     (Cobre os casos que antes eram corrigidos "à mão" por id: D2089, R$ 366.666 em I22_2; D2150, R$ 1.000.000 em I21.)
n_trava = 0
for col, cls in zip(COLS_VALOR, COLS_CLASSE):
    mk = df_moradores[col] >= 100000.0
    n_trava += int(mk.sum())
    df_moradores.loc[mk, col] = np.nan
    df_moradores.loc[mk, cls] = CLASSE_ERRO
print(f"Fontes individuais >= R$ 100 mil convertidas em NaN: {n_trava}")


# %% [7] Renda individual, domiciliar e per capita ------------------------------------------------------
# Regra: se QUALQUER fonte do morador é recusa/não sabe/erro, a renda individual é NaN (perda declarada);
# caso contrário, é a soma das fontes com valor >= 0 (0,0 se nenhuma fonte tiver valor).
CLASSES_PERDA = ["Recusa", "Não sabe", "Não Sabe", CLASSE_ERRO]
val = df_moradores[COLS_VALOR]
cls = df_moradores[COLS_CLASSE]
algum_perdido = cls.isin(CLASSES_PERDA).any(axis=1)
soma = val.where(val.notna() & (val >= 0)).sum(axis=1, min_count=1).fillna(0.0)
df_moradores["renda_individual_total"] = soma.where(~algum_perdido)

# Domicílio = soma dos moradores com renda conhecida (min_count=1: NaN só se TODOS forem NaN).
# ATENÇÃO: é um piso. Se um morador tem renda NaN e outro não, a soma ignora o primeiro.
g = df_moradores.groupby("domicilio_id")["renda_individual_total"]
agreg = pd.DataFrame({
    "renda_domiciliar_calculada": g.sum(min_count=1),
    "n_moradores_renda_ausente": g.apply(lambda s: int(s.isna().sum())),
}).reset_index()
agreg["renda_dom_parcial"] = (agreg["n_moradores_renda_ausente"] > 0) & agreg["renda_domiciliar_calculada"].notna()

for c in ["renda_domiciliar_calculada", "renda_domiciliar_per_capita", "n_moradores_renda_ausente", "renda_dom_parcial"]:
    if c in df_domicilios.columns:
        df_domicilios = df_domicilios.drop(columns=[c])
df_domicilios = df_domicilios.merge(agreg, on="domicilio_id", how="left")

n_linhas = df_moradores.groupby("domicilio_id").size()
assert (n_linhas.reindex(df_domicilios["domicilio_id"]).values == df_domicilios["npessoas_qtd_de_pessoas_no_domicilio"].values).all(), \
    "npessoas difere do número de moradores: o denominador da renda per capita estaria errado"
df_domicilios["renda_domiciliar_per_capita"] = (
    df_domicilios["renda_domiciliar_calculada"] / df_domicilios["npessoas_qtd_de_pessoas_no_domicilio"])
print("Renda individual, domiciliar e per capita calculadas.")


# %% [8] Saneamento final de texto e nomes de colunas -----------------------------------------------------
for col in ["E04_1_1", "E04_2_1"]:
    if col in df_moradores.columns:
        df_moradores[col] = df_moradores[col].astype(str).replace(
            {"99999": "Não se aplica", "99999.0": "Não se aplica", "nan": "Não se aplica", "None": "Não se aplica"})
if "E23_1" in df_moradores.columns:
    df_moradores["E23_1"] = df_moradores["E23_1"].astype(str).replace(
        {"99999": "Não se aplica", "99999.0": "Não se aplica", "88888": "Não sabe", "88888.0": "Não sabe", "nan": "Não se aplica"})
for df in (df_moradores, df_domicilios):
    if "insegurança_alimentar_x" in df.columns:
        df.drop(columns=["insegurança_alimentar_x"], inplace=True)

df_moradores.rename(columns={
    "E02_dt_nascimento": "E02_dt_nascimento__Data_de_Nascimento_dd_mm_aaaa", "E02_idade": "E02_idade__Idade_informada",
    "E03_sexo": "E03_sexo__Sexo_informado", "E04_1_1": "E04_1_1__Descricao_outro_genero",
    "E04_2_1": "E04_2_1__Descricao_outra_orientacao_sexual", "E16pais": "E16_pais__Pais_de_nascimento",
    "E16uf": "E16_uf__UF_de_nascimento", "E16municipio": "E16_municipio__Municipio_de_nascimento",
    "E18pais": "E18_pais__Pais_ultima_moradia", "E18uf": "E18_uf__UF_ultima_moradia",
    "E18municipio": "E18_municipio__Municipio_ultima_moradia", "E21": "E21__RA_onde_morava_em_julho_2021",
    "E21_pais": "E21_pais__Pais_onde_morava_em_julho_2021", "E21_estado": "E21_estado__UF_onde_morava_em_julho_2021",
    "E21_municipio": "E21_municipio__Municipio_onde_morava_em_julho_2021",
    "E23_1": "E23_1__Local_intencao_novo_domicilio", "G02": "G02__Local_atendimento_saude",
    "H05": "H05__Local_unidade_ensino", "I08": "I08__Local_trabalho_principal",
    "J04_1": "J04_1__Mes_ano_nascimento_ultimo_filho",
}, inplace=True)

for col in ["C05_1_1_acesso_por_internet_banda_larga_fixa", "C05_1_2_acesso_por_internet_de_rede_de_celular"]:
    if col in df_domicilios.columns:
        df_domicilios[col] = df_domicilios[col].astype(str).replace(
            {"99999": "Não se aplica", "99999.0": "Não se aplica", "nan": "Não se aplica", "None": "Não se aplica"})
for col in ["B14_8_a_agua_possui_cheiro_gosto_ou_sabor_diferente", "B14_9_possui_outorga_para_uso_da_agua"]:
    if col in df_domicilios.columns:
        df_domicilios[col] = df_domicilios[col].astype(str).replace(
            {"99999": "Não se aplica", "99999.0": "Não se aplica", "88888": "Não sabe", "88888.0": "Não sabe",
             "77777": "Recusa", "77777.0": "Recusa", "nan": "Não se aplica"})

for col in ["a01ra", "a01subpopulacao", "a01upt"]:
    for df in (df_moradores, df_domicilios):
        if col in df.columns:
            df[col] = df[col].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
print("Saneamento final concluído.")


# %% [9] Auditoria (contagens reais e verificações que interrompem a execução) -----------------------------
print("\n=== Células por classe nas seis colunas de renda (moradores) ===")
resumo = pd.DataFrame({c.split("__")[0]: df_moradores[f"{c}_classe"].value_counts(dropna=False) for c in COLUNAS_RENDA})
display(resumo.fillna(0).astype(int))

# 9.1 Nenhum código administrativo, quase-código, negativo ou valor >= R$ 100 mil em _valor
CODIGOS_NUM = {66666, 77777, 88888, 99999, 77774, 36666, 66663, 88885, 88880,
               666666, 777777, 888888, 999999, 7777777, 8888888, 9999999}
for c in COLS_VALOR:
    v = df_moradores[c].dropna()
    assert not v.isin(CODIGOS_NUM).any(), f"Código administrativo remanescente em {c}"
    assert (v >= 0).all(), f"Valor negativo em {c}"
    assert (v < 100000).all(), f"Valor >= R$ 100 mil remanescente em {c}"


# 9.2 Varredura de quase-códigos (distância de edição 1 de um código) entre os valores aceitos
def dist1(a, b):
    if len(a) == len(b):
        return sum(x != y for x, y in zip(a, b)) <= 1
    if abs(len(a) - len(b)) == 1:
        s, l = sorted((a, b), key=len)
        return any(l[:i] + l[i + 1:] == s for i in range(len(l)))
    return False


for c in COLUNAS_RENDA:
    ok = df_moradores[df_moradores[f"{c}_classe"] == CLASSE_OK][c].astype(str).str.strip().str.split(".").str[0]
    sus = ok[ok.str.match(r"^\d{4,}$") & ok.apply(lambda x: any(dist1(x, k) and x != k for k in
                                                              ("66666", "77777", "88888", "99999")))]
    assert sus.empty, f"Possível quase-código aceito como valor em {c}: {sus.unique().tolist()}"

# 9.3 Coerência da renda individual
assert (df_moradores["renda_individual_total"].isna() == algum_perdido).all()
assert (df_moradores["renda_individual_total"].dropna() >= 0).all()

# 9.4 Ausência de renda e domicílios parciais
n_dom = len(df_domicilios)
n_nan = df_domicilios["renda_domiciliar_per_capita"].isna().sum()
n_par = int(df_domicilios["renda_dom_parcial"].sum())
print(f"\nMoradores com renda individual NaN: {df_moradores['renda_individual_total'].isna().sum()} "
      f"({df_moradores['renda_individual_total'].isna().mean():.1%})")
print(f"Domicílios sem renda (todos NaN): {n_nan} ({n_nan / n_dom:.1%} de {n_dom})")
print(f"Domicílios com renda PARCIAL (algum morador NaN; renda é um piso): {n_par} ({n_par / n_dom:.1%})")
assert n_nan / n_dom < 0.50, "Mais da metade dos domicílios sem renda"

pc = df_domicilios["renda_domiciliar_per_capita"].dropna()
print(f"Renda per capita (R$): média {pc.mean():,.2f} | mediana {pc.median():,.2f} | máx {pc.max():,.2f} | "
      f"variância do log1p {np.log1p(pc).var():.4f}")

# 9.5 Valores extremos aceitos (>= R$ 20.000 em uma fonte) — para revisão manual, não são alterados
ext = []
for c in COLUNAS_RENDA:
    v = df_moradores[f"{c}_valor"]
    for i in v[v >= 20000].index:
        ext.append((df_moradores.at[i, "domicilio_id"], c.split("__")[0], v[i]))
print(f"\nValores aceitos >= R$ 20.000 por fonte ({len(ext)}); maiores: "
      f"{sorted(ext, key=lambda t: -t[2])[:8]}")


# %% [9b] Pareto: códigos administrativos/erros por variável de renda (contagem real, após o saneamento) ------
try:
    import matplotlib.pyplot as plt
    par = pd.Series({c.split("__")[0]: int((df_moradores[f"{c}_classe"] == CLASSE_ERRO).sum()) for c in COLUNAS_RENDA})
    par = par.sort_values(ascending=False)
    acum = par.cumsum() / par.sum() * 100
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.bar(par.index, par.values, color="C0")
    ax1.set_ylabel("Células tratadas como erro/recusa (NaN)")
    ax2 = ax1.twinx()
    ax2.plot(par.index, acum.values, color="C1", marker="D")
    ax2.axhline(80, color="red", ls="--", alpha=.5)
    ax2.set_ylabel("Percentual acumulado (%)")
    plt.title("Pareto: erros e recusas nas variáveis de renda (Série I)")
    plt.show()
except Exception as e:                  # o gráfico é opcional
    print("Pareto não gerado:", e)


# %% [10] Exportação -----------------------------------------------------------------------------------
if EXPORTAR:
    os.makedirs(OUT_DIR, exist_ok=True)
    try:
        import xlsxwriter  # noqa: F401
        motor = "xlsxwriter"
    except ImportError:
        motor = "openpyxl"
    with pd.ExcelWriter(os.path.join(OUT_DIR, "pdadr_moradores.xlsx"), engine=motor) as w:
        df_moradores.to_excel(w, index=False)
    with pd.ExcelWriter(os.path.join(OUT_DIR, "pdadr_domicilios.xlsx"), engine=motor) as w:
        df_domicilios.to_excel(w, index=False)
    print(f"Arquivos exportados em {OUT_DIR}")
