# Rastreabilidade: avaliação da orientadora (12/11/2026) → código → saída

Cada ponto da avaliação aponta para o script que o executa e para o arquivo (por variante, em `resultados/<variante>/`) que o prova.
"Status" diz apenas se o teste existe e roda; **não diz que o resultado foi favorável**. A comparação entre variantes está em `resultados/comparacao_variantes.md`; a leitura detalhada, em `docs/rodadas_stress_test.md` (tabelas em `docs/tabelas_rodadas.md`). O texto do TCC adota M2 como principal e M0, M1 e M2b como comparação e sensibilidade (ver `docs/criticas_e_modelo_hibrido.md`).

| # | Ponto da avaliação | Script | Saída | Status |
|---|---|---|---|---|
| 1 | IVS-DF e escala da RA fora da matriz | `03` (flag `ctx`) | `config.txt`, `saida_03.txt` | implementado; testado como `A_so_ctx` e `L_sem_ctx` |
| 2 | Água e esgoto por escada de serviço (PNSR/JMP), em ilr ou clr | `03` (flags `agua`, `esgoto`; `COMP=ilr|clr`) | `S_clr`, `A_so_agua`, `A_so_esgoto` | implementado |
| 3 | Renda por média de ln(1+x) com winsorização no P99 | `03` (flag `renda`; `WINSOR_P`) | `S_winsor95`, `S_winsor100` | implementado (mediana ponderada **não** implementada: a avaliação aceita "mediana **ou** média de log") |
| 4 | Geladeira e secadora fora; D1 = bens duráveis | `03` (flag `dur`) | `A_so_dur`, `L_sem_dur` | implementado |
| 5 | Bootstrap de domicílios refazendo todo o fluxo, com perfis e ICs | `16` (a) | `saida_16.txt`, `perfis_grupos_ic.csv` | implementado; `NINIT_BOOT` separa ótimo local de incerteza amostral |
| 6 | Encolhimento e/ou "não classificável" | `16` (d), (e) | `robustez_amostral_segmentos.csv` | implementado |
| 7 | ARI de estabilidade contra o nulo de covariância | `16` (c) | `estabilidade_vs_nulo_cov.csv` | implementado |
| 8 | Destino do grupo 2 decidido pelo bootstrap | `16`, `21` S4 | `saida_16.txt`, `stress_S4_*.csv` | resultado por variante; no texto do TCC, o K final é 3 e a regra é: segmentos com menos de 10 domicílios são não classificáveis e os de coligação < 0,60, de fronteira |
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
| E6 | K alternativo nomeado (K = 3 ou 4) | `variantes.tsv`, `03` | `M0_antigo_K3/K4/K5`, `M1_revisado_K3/K4/K5` | K = 2 a 8 comparados em todas as variantes; K=3 e K=4 nomeados para M0 e M1 |
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

## Rodadas adicionadas depois da auditoria (programas 24, 25, 26, 27)

| Ponto | Programa | Saída | Status |
|---|---|---|---|
| Nulo do bootstrap (ARI observado contra dados sem grupos, mesma covariância) | `24` (1) | `nulos_24_resumo.csv` | implementado; resolução do p = 0,01 (100 matrizes nulas) |
| Replicabilidade por metades dos domicílios | `24` (2) | `nulos_24_resumo.csv` | implementado |
| Nulo de atribuição domiciliar dentro do estrato (AMI com RA e subpopulação) | `24` (3) | `nulos_24_resumo.csv` | implementado |
| Multinível por segmento (CCI da RA), Moran I ajustado pela composição, validação ordenada | `25` | `multinivel_25.csv`, `moran_ajustado_25.csv`, `validacao_ordenada_25.csv` | implementado; sem correção para comparações múltiplas |
| Decomposição de variância no nível do domicílio (logit misto por VB; ln da renda em MixedLM) | `26` | `decomposicao_domiciliar_26*.csv` | implementado em Python; sem pesos amostrais; VB subestima a incerteza |
| Encolhimento por n efetivo de Kish (variante d2) | `16` | `robustez_amostral_segmentos.csv` | implementado |
| Bootstrap por estrato e retenção de componentes redecidida por réplica | `16` (`BOOT_NIVEL`, `BOOT_RETENCAO`) | `*_bootEstratoRegra` | implementado para M0 K=5 e M1 K=3 |
| Escadas alternativas de água e esgoto | `03` (`AGUA_ESCADA`, `ESGOTO_ESCADA`) | variantes `E_*` | implementado |
| K = 4 para M0 e M1 | `variantes.tsv` | `M0_antigo_K4`, `M1_revisado_K4` | implementado |
| Tabelas T1 a T8 e ARI entre partições | `27` | `docs/tabelas_rodadas.md`, `resultados/ari_entre_particoes.csv` | implementado |
| Validação externa sob incerteza da partição | `28` | `validacao_externa_incerteza_28*.csv` | implementado; desfechos fixos, só a partição varia; sem correção para múltiplos testes |
| Jackknife por domicílio; Emater-DF | – | – | **não feito** |

| Validação externa condicional à RA (permutação dentro da RA) e eta² da RA | `29` | `validacao_condicional_RA_29*.csv`, `eta2_RA_escores_29.csv` | implementado |
| Híbrido M2 (`ctx`, `renda`, `dur`; água e esgoto como em M0), K = 3, 4, 5 | `variantes.tsv` | `M2_semEscadas_K3/K4/K5`; análise em `docs/criticas_e_modelo_hibrido.md` | implementado; definido depois de ver `L_sem_agua` |

| Água como ilr das três fontes originais, sem escada (M2b), K = 3 | `03` (`AGUA_RAWILR=1`), `variantes.tsv` | `M2b_aguaILR_K3` | implementado e executado com a bateria completa em K = 3 (`reproduzir_tudo.py --variantes M2b_aguaILR_K3`, mais `09`, `12` e `29` à parte); K = 4 e 5 não executados |
| Programas 09 e 12 sensíveis à configuração da rodada (`K_FINAL`, `V2`, `AGUA_RAWILR`) | `09`, `12` | `kw_todas.csv`, `robustez_tandem.csv` | corrigido: o 12 usava as variáveis do M1 no M2 e os números da sensibilidade à PCA foram refeitos |
