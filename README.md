# Tipologia dos territórios rurais do DF — PDAD Rural 2022 (reprodutibilidade e testes de estresse)

Repositório de apoio ao TCC *Tipologia dos Territórios Rurais do Distrito Federal: Uma Análise de Clustering Socioespacial com a PDAD Rural 2022* (IESB, Ciência de Dados e IA).
Ele contém o código, os dados e as saídas necessárias para reproduzir **cada modelo** discutido no TCC e para **testá-los sob estresse com a mesma bateria**, conforme a avaliação crítica de 12/11/2026.

> **Nenhum modelo foi eleito.** O repositório existe para comparar o modelo avaliado (`M0_antigo`) com o modelo revisado segundo a avaliação (`M1_revisado`) e com variantes intermediárias, e para mostrar quais conclusões se sustentam em todas elas.

## Estrutura

```
dados/gold/            microdados corrigidos (saída do programa 01) + IVS-DF (IPEDF)
dados/                 versão anterior das planilhas (sem as correções dos códigos 66666 etc.) e shapefile das RAs
legado/                scripts originais dos modelos 4D e 5D (Colab); mantidos para consulta
pipeline/              programas numerados e executores (veja abaixo)
resultados/<variante>/ saídas de cada variante (texto e CSV; figuras só de M0_antigo_K5 e M1_revisado_K3)
resultados/comparacao_variantes.{csv,md}   métricas de todas as variantes lado a lado
docs/                  rastreabilidade da avaliação, denominadores e notas
```

## Como reproduzir

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt          # versões originais do projeto (Colab, 2025)
# Testado em Python 3.13 com: pandas 3.0, numpy 2.5, scikit-learn 1.9, scipy 1.18, statsmodels 0.15, geopandas 1.2 (ver requirements-testado.txt)

bash pipeline/rodar_variante.sh M1_revisado_K3   # uma variante (cerca de 8 min)
bash pipeline/rodar_todas.sh 2                   # todas as variantes de pipeline/variantes.tsv, 2 em paralelo
python pipeline/22_comparar_variantes.py         # consolida em resultados/comparacao_variantes.md
```

Todas as sementes são fixas (`SEED = 42`; réplicas do bootstrap usam a semente da réplica). Os resultados do K-Means dependem do número de inicializações; o executor registra os parâmetros em `config.txt`.

### Variantes (`pipeline/variantes.tsv`)

| Grupo | Variantes | O que muda |
|---|---|---|
| Modelo avaliado | `M0_antigo_K5`, `M0_antigo_K3` | matriz original: 34 variáveis, IVS-DF e escala da RA dentro, água em 3 proporções, geladeira e secadora dentro, renda = ln da média |
| Modelo revisado | `M1_revisado_K3`, `M1_revisado_K5` | as 5 mudanças da avaliação juntas (flags `ctx,agua,esgoto,renda,dur`) |
| Uma mudança por vez | `A_so_ctx`, `A_so_agua`, `A_so_esgoto`, `A_so_renda`, `A_so_dur` | M0 + uma só mudança |
| Tira uma mudança | `L_sem_ctx`, `L_sem_agua`, `L_sem_esgoto`, `L_sem_renda`, `L_sem_dur` | M1 − uma mudança |
| Sensibilidade | `S_winsor95`, `S_winsor100`, `S_clr` | winsorização em P95 ou ausente; clr (3 colunas) no lugar de ilr |
| Mais réplicas | `*_b500n50` | bootstrap com 500 réplicas e 50 inicializações |

Flags do modelo (`V2`): `ctx` (IVS-DF e escala da RA saem da matriz), `agua` e `esgoto` (escada de serviço rede / individual adequado / precário, em ilr), `renda` (média de ln(1+renda) com winsorização), `dur` (geladeira e secadora saem; D1 vira "bens duráveis").

## Programas (`pipeline/`)

| Programa | Função |
|---|---|
| `01_saneamento_pdadr_2022.py` | saneamento da planilha-dicionário do IPEDF → `dados/gold/*.xlsx`. **Requer `Dicionario_aplicado_PDADR_2022.xlsx` (não incluído; peça ao autor ou ao IPEDF).** As saídas já estão em `dados/gold/`. |
| `02_analise_exploratoria_ponderada.py` | análise exploratória ponderada |
| `03_analise_cluster_4d.py` | **modelo**: agregação por segmento, PCA por dimensão, K-Means, nulos, co-associação, Moran/LISA. Parâmetros por variáveis de ambiente (`V2`, `K_FINAL`, `WINSOR_P`, `COMP`) |
| `04`–`13` | comparação de métodos, Ward, estabilidade cruzada, variáveis de apoio, Kruskal-Wallis, tabelas, mapas, robustez em tandem, perfis por RA |
| `16_robustez_amostral_e_pesos.py` | bootstrap de domicílios (fluxo inteiro refeito), ICs dos perfis, estabilidade contra nulo de covariância, encolhimento, "não classificável", K-Means ponderado |
| `17_ablacao_matriz.py` | ablação das mudanças da matriz |
| `18_variancia_ami_espacial.py` | decomposição de variância RA × subpopulação, AMI, join count |
| `19_validacao_externa.py` | EBIA, CadÚnico, benefícios, escritura, outorga |
| `20_figuras_perfil_validacao.py` | figuras de perfis e validação |
| `21_stress_estrutural.py` | testes de estresse S1–S10 (sementes, tamanho × extremidade, Wilson, aparado, arquétipos, deixa-um-fora, algoritmos, corte de n, puxão da RA, eixos) |
| `24_nulos_bootstrap_metades_territorial.py` | nulo do ARI do bootstrap, replicabilidade por metades de domicílios e nulo de atribuição domiciliar dentro do estrato (AMI de RA/subpopulação) |
| `25_multinivel_moran_ajustado_validacao_ordenada.py` | MixedLM (subpopulação fixa, RA aleatória), Moran I ajustado por composição e validação externa ordenada |
| `26_decomposicao_domiciliar_mista.py` | decomposição no nível do domicílio como a avaliação especifica: logit misto (binários) e MixedLM (renda), efeito fixo de subpopulação + RA + RA×subpopulação aleatórios |
| `22_comparar_variantes.py` | consolidação das variantes |

## Dados

* `dados/gold/pdadr_domicilios.xlsx` e `pdadr_moradores.xlsx`: microdados da PDAD Rural 2022 (IPEDF) com as correções de códigos administrativos (66666, 77777, 88888, 99999 e variantes) documentadas no programa 01. São as planilhas usadas nos resultados do TCC.
* `dados/pdadr_*.xlsx` (raiz): **versão anterior**, sem essas correções; mantida só porque os scripts de `legado/` as leem.
* `dados/gold/Base_dados-IVS.xlsx`: IVS-DF (IPEDF); usado **apenas para descrição** no modelo revisado.
* `dados/regioes_administrativas.zip`: limites das RAs.

Fonte dos dados: IPEDF, PDAD Rural 2022 (https://www.ipe.df.gov.br/).

## Limitações conhecidas

* O programa 01 não roda sem a planilha-dicionário (silver), que não está no repositório.
* A mediana ponderada de renda por segmento não foi implementada (a avaliação aceita a média de ln com winsorização).
* O K-Means aparado e a análise de arquétipos são implementações próprias, não os pacotes `tclust` e `archetypes` do R.
* Na validação deixa-uma-RA-fora, o PCA é mantido fixo (só o K-Means é refeito).
* Os resultados dependem da versão das bibliotecas; as diferenças esperadas são pequenas (ordem de rótulos), mas não foram medidas.


Auditoria contra a avaliação: `docs/auditoria_pipeline_vs_avaliacao.md`. A coluna opcional `extras` de `pipeline/variantes.tsv` aceita `BOOT_NIVEL=estrato;BOOT_RETENCAO=regra`, `AGUA_ESCADA=...`, `ESGOTO_ESCADA=...`.
