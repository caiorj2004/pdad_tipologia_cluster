# Rastreabilidade: avaliação da orientadora (12/11/2026) → código → saída

Cada ponto da avaliação aponta para o script que o executa e para o arquivo (por variante, em `resultados/<variante>/`) que o prova.
"Status" diz apenas se o teste existe e roda; **não diz que o resultado foi favorável**. A comparação entre variantes está em `resultados/comparacao_variantes.md` e nenhum modelo foi eleito.

| # | Ponto da avaliação | Script | Saída | Status |
|---|---|---|---|---|
| 1 | IVS-DF e escala da RA fora da matriz | `03` (flag `ctx`) | `config.txt`, `saida_03.txt` | implementado; testado como `A_so_ctx` e `L_sem_ctx` |
| 2 | Água e esgoto por escada de serviço (PNSR/JMP), em ilr ou clr | `03` (flags `agua`, `esgoto`; `COMP=ilr|clr`) | `S_clr`, `A_so_agua`, `A_so_esgoto` | implementado |
| 3 | Renda por média de ln(1+x) com winsorização no P99 | `03` (flag `renda`; `WINSOR_P`) | `S_winsor95`, `S_winsor100` | implementado (mediana ponderada **não** implementada: a avaliação aceita "mediana **ou** média de log") |
| 4 | Geladeira e secadora fora; D1 = bens duráveis | `03` (flag `dur`) | `A_so_dur`, `L_sem_dur` | implementado |
| 5 | Bootstrap de domicílios refazendo todo o fluxo, com perfis e ICs | `16` (a) | `saida_16.txt`, `perfis_grupos_ic.csv` | implementado; `NINIT_BOOT` separa ótimo local de incerteza amostral |
| 6 | Encolhimento e/ou "não classificável" | `16` (d), (e) | `robustez_amostral_segmentos.csv` | implementado |
| 7 | ARI de estabilidade contra o nulo de covariância | `16` (c) | `estabilidade_vs_nulo_cov.csv` | implementado |
| 8 | Destino do grupo 2 decidido pelo bootstrap | `16`, `21` S4 | `saida_16.txt`, `stress_S4_*.csv` | resultado por variante; decisão final **em aberto** |
| 9 | Decomposição de variância RA × subpopulação (multinível) | `18` | `variancia_ra_sub.csv` | implementado (ANOVA e MixedLM) |
| 10 | AMI dos rótulos com RA e subpopulação | `18` | `ami_rotulos.csv` | implementado |
| 11 | Média por RA ponderada; join count dos grupos | `03`, `18` | `moran_global.csv`, `joincount.csv` | implementado |
| 12 | Validação externa (EBIA, benefícios, posse do imóvel) | `19` | `validacao_externa*.csv` | implementado; posse do imóvel = escritura (B07) e documentação (B04) |
| 13 | Validação com técnicos da Emater-DF | – | – | **pesquisa futura** (permitido pela avaliação) |
| 14 | Protocolo de busca bibliográfica | (texto do TCC, Cap. 3) | – | no texto; bases pagas (Scopus/WoS) **não** consultadas |
| 15 | Interpretação como hipóteses; nomes substantivos | (texto do TCC) | – | no texto |
| E1 | Spearman tamanho amostral × distância ao centro (C1) | `21` S2 | `stress_S2_tamanho_vs_ponta.csv` | implementado |
| E2 | ICs de Wilson dos segmentos com n < 12 | `21` S3 | `stress_S3_wilson_pequenos.csv` | implementado |
| E3 | K-Means ponderado | `16` (b), `21` S7 | `robustez_pesos_*.csv`, `stress_S7_algoritmos.csv` | implementado |
| E4 | K-Means aparado (trimmed), α = 6% e 10% | `21` S4 | `stress_S4_kmeans_aparado.csv` | implementado (implementação própria, sem o pacote `tclust` do R) |
| E5 | Análise de arquétipos | `21` S5 | `stress_S5_arquetipos*.csv` | implementado (implementação própria por NNLS alternado) |
| E6 | K alternativo nomeado (K = 3 ou 4) | `variantes.tsv`, `03` | `M1_revisado_K3`, `M1_revisado_K5` | K = 2 a 8 comparados em todas as variantes |
| E7 | Denominador da CTPS explícito | `03` (linha com `tem_ctps`) | `docs/denominadores.md` | documentado |
| E8 | Reprodutibilidade: links dos códigos e dos dados | este repositório | `README.md` | **Falta o DOI/arquivamento (Zenodo) e o link do artefato** |
| E9 | Interação RA × subpopulação (efeito mecânico da RA) | `21` S9 | `stress_resumo.csv` | implementado |

## Rodadas de estresse (todas aplicadas a cada variante)

1. **Validade interna**: silhueta, CH, DB, nulo embaralhado e nulo de covariância, para K = 2 a 8 (`03`).
2. **Amostragem**: subamostras de 80% dos segmentos; bootstrap de domicílios com o fluxo inteiro refeito; estabilidade contra nulo de covariância (`03`, `16`).
3. **Algoritmo**: sementes únicas, Ward, ligação média, GMM, K-Means ponderado, aparado, arquétipos (`21` S1, S4, S5, S7).
4. **Estrutura do dado**: tamanho × extremidade, corte de n, deixa-uma-RA-fora, deixa-uma-subpopulação-fora (`21` S2, S6, S8).
5. **Confusão territorial**: AMI, variância, puxão da RA, join count, Moran (`18`, `21` S9).
6. **Construção da matriz**: ablação de cada ajuste (A_so_*, L_sem_*), winsorização e clr (`variantes.tsv`).
7. **Validade externa**: EBIA, CadÚnico, benefícios, escritura, outorga, IVS-DF (`19`).
