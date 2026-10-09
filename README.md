# 🚜 Agrotec-DF: Tipologia e Clusterização Territorial

Este projeto integra o esforço de inteligência territorial do **Agrotec-DF**, focado em mapear e analisar o ecossistema de agricultura familiar no Distrito Federal utilizando os microdados da **PDAD Rural 2022** (Pesquisa Distrital por Amostra de Domicílios).
O objetivo principal desta modelagem é aplicar técnicas robustas de aprendizado de máquina não supervisionado e estatística espacial para identificar agrupamentos (clusters) de vulnerabilidade e infraestrutura nas Regiões Administrativas (RAs) do DF.

---

## 🏗️ Estrutura do Repositório

```text
pdad_tipologia_cluster/
│
├── dados/                              # Diretório de dados brutos e shapes
│   ├── pdadr_domicilios.xlsx           # Microdados domiciliares (PDAD Rural)
│   ├── pdadr_moradores.xlsx            # Microdados de moradores (PDAD Rural)
│   └── regioes_administrativas.zip     # Shapefiles das fronteiras das RAs do DF
│
├── modelo_4d_referencia.py             # Pipeline base com 4 dimensões de capital
├── modelo_5d_expandido.py              # Pipeline expandido com dimensão Social
├── requirements.txt                    # Dependências e versões do ecossistema
└── README.md                           # Documentação principal

```

---

## 🧠 Modelos Desenvolvidos

A pesquisa propõe duas abordagens modulares para a construção da tipologia territorial, baseadas na teoria dos Capitais:

1. **Modelo 4D (Referência):**
Concentra-se na extração de sinais de quatro capitais primários:
* **Financeiro/Ativos:** Renda domiciliar per capita e posse de bens.
* **Humano/Trabalho:** Escolaridade adulta, formalização do trabalho e IVS (Índice de Vulnerabilidade Social).
* **Físico/Habitacional:** Infraestrutura de água, esgoto, energia, pavimentação e adensamento.
* **Tecnológico/Conectividade:** Acesso à internet fixa/móvel e posse de computadores.


2. **Modelo 5D (Expandido):**
Integra uma quinta dimensão para avaliar o impacto de políticas públicas:
* **Social/Institucional:** Recebimento de benefícios sociais, valores monetários de auxílios e incentivos à produção. *(Nota metodológica: Esta dimensão possui retenção forçada no pipeline de redução de dimensionalidade).*

---

## 🔬 Metodologia Estatística e Espacial

Para garantir a ausência de viés estocástico e a validade estrutural dos agrupamentos, o pipeline implementa as seguintes etapas automatizadas:

* **Ponderação Populacional:** Expansão amostral dos domicílios respeitando o universo $N_h$ dos 21 estratos geográficos do IPEDF. Supressão de segmentos com $n < 5$ para evitar ruído.
* **Extração de Sinais Analíticos (PCA + Monte Carlo):** Em vez de reter componentes principais por variância arbitrária, o algoritmo utiliza **Análise Paralela de Monte Carlo** (1.000 iterações iterativas), retendo apenas os eixos (autovalores) que superam o percentil 95% do ruído aleatório.
* **Clusterização K-Means Dinâmica:** Agrupamento das Regiões Administrativas otimizado para $K=6$, baseado em métricas de validação intrínseca (Silhueta, Davies-Bouldin, Calinski-Harabasz).
* **Autocorrelação Espacial Local (LISA):** Validação topológica utilizando matrizes de contiguidade *Queen* (fronteiras reais via Geopandas/PySAL), com correção rigorosa de múltiplos testes **FDR (False Discovery Rate)** de Benjamini-Hochberg ($\alpha = 0.05$) para controle de falsos positivos espaciais.

---

## 🚀 Como Executar Localmente

### 1. Clonar o Repositório

Abra o terminal e baixe o código fonte:

```bash
git clone https://github.com/SEU_USUARIO/pdad_tipologia_cluster.git
cd pdad_tipologia_cluster

```

### 2. Configurar o Ambiente Virtual

Recomenda-se o uso de um ambiente virtual para isolar as dependências geoespaciais, que costumam ser sensíveis a conflitos de versão:

```bash
python -m venv venv

# Ativar no Windows:
venv\Scripts\activate
# Ativar no Linux/Mac:
source venv/bin/activate

```

### 3. Instalar Dependências

```bash
pip install -r requirements.txt

```

### 4. Executar os Modelos

Para rodar o Modelo de Referência e gerar as estatísticas e o mapa coroplético:

```bash
python modelo_4d_referencia.py

```

Para rodar a versão com a dimensão institucional:

```bash
python modelo_5d_expandido.py

```

---

## 🛠️ Stack Tecnológico

* **Manipulação e Engenharia de Dados:** `pandas`, `numpy`
* **Geoprocessamento e Estatística Espacial:** `geopandas`, `libpysal`, `esda`
* **Machine Learning:** `scikit-learn` (PCA, KMeans, GaussianMixture)
* **Inferência e Validação Estatística:** `scipy`, `statsmodels` (FDR, Kruskal-Wallis)
* **Visualização:** `matplotlib`, `seaborn`
