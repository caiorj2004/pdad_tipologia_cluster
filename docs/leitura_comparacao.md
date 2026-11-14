# Leitura da comparação de variantes (sem eleger vencedor)

Base: `resultados/comparacao_variantes.{csv,md}`, 27 variantes completas (etapas 03, 16, 18, 19, 21, 23, 24, 25, 26), geradas em 2026-11-13. Esta leitura descreve o que os números mostram; **não escolhe modelo**.

## O que é igual em todas as variantes
- Multinível por segmento, Moran I (bruto ≈ 0,003; ajustado pela composição ≈ 0) e decomposição domiciliar (RA ≈ 5,7%, subpopulação ≈ 0,2%, segmento ≈ 12%, residual ≈ 82%) dependem dos dados, não da partição. Não discriminam modelos; descrevem o território.
- Winsorização a 1,00 vs 0,99 vs 0,95 não muda a partição de M1 K=3 (tamanhos 24/17/8 nos três casos). A escolha de winsor não decide nada.
- Aumentar o bootstrap de 200 para 500 réplicas (50 inicializações) mexe pouco no ARI médio (M0 K=5: 0,458 → 0,472; M1 K=3: 0,381 → 0,381). O Monte Carlo está estável.
- Bootstrap por estrato com regra de retenção redecidada por réplica dá resultado próximo ao padrão (M1 K=3: 0,381 vs 0,373).

## Contra os nulos
- Nenhum K final bate o nulo de covariância com p < 0,05 em M1. Em M0 K=5, p = 0,110. Mesmo sem esse corte, A_so_dur, A_so_renda e A_so_agua têm alguns K isolados abaixo de 0,05.
- Contra o nulo do bootstrap, o ARI observado não passa do p95 do nulo em nenhuma variante principal (M1 K=3: 0,381 vs p95 0,687; M0 K=5: 0,458 vs 0,551). A estabilidade dos grupos não se distingue do que se obteria com ruído de mesma covariância.
- Metades (split-half): ARI observado médio de 0,23 a 0,32, abaixo do p95 do nulo em todas.
- Estabilidade de subamostra a 80%: de 0,48 (M1 K=4) a 0,71. K=2 é sempre mais estável (0,86 a 0,98) e tem a maior silhueta (≈ 0,31 a 0,33).

## Entre os modelos centrais
| | M0 K=5 | M0 K=3 | M1 K=3 | M1 K=4 | M1 K=5 |
|---|---|---|---|---|---|
| silhueta | 0,236 | 0,221 | 0,253 | 0,224 | 0,212 |
| ARI bootstrap médio | 0,458 | 0,558 | 0,381 | 0,396 | 0,368 |
| ARI subamostra 80% | 0,664 | 0,658 | 0,658 | 0,480 | 0,566 |
| menor grupo (% dom.) | 0,5 | 22,4 | 22,5 | 18,0 | 8,1 |
| ARI deixando uma RA fora (refazendo o PCA) | 0,805* | 0,903* | 0,736 | 0,457 | 0,664 |

*Valores de S6b na coluna `S6b_ARI_medio_sem_uma_RA`.

- M0 K=5 tem um grupo com 0,5% dos domicílios expandidos (3 segmentos). Isso é fragilidade de interpretação mesmo com métricas aceitáveis.
- M1 K=3 tem a maior silhueta do conjunto de K final, mas o menor ARI de bootstrap entre os K=3. M0 K=3 tem bootstrap mais alto e silhueta menor. Há um trade-off, não um dominante.
- M1 K=4 e K=5 perdem em subamostra e em deixar-uma-RA-fora.

## Estrutura espacial e composição
- AMI entre partição e RA é baixa (0,05 a 0,16); com subpopulação, em geral similar. Em M1 K=3, o AMI com RA não passa do nulo de atribuição domiciliar (p = 0,68), enquanto o de subpopulação fica perto do limiar (p = 0,06). Em M0 K=5 é o inverso (RA p = 0,05; subpop p = 0,29). Não há leitura "territorial" ou "setorial" única.
- Correlação de Spearman entre tamanho do segmento e distância ao centro: negativa em todas as variantes (−0,33 a −0,63; p ≤ 0,02). Segmentos pequenos tendem a ficar mais longe dos centros. Isso persiste em todo o conjunto e precisa ser dito no texto.
- Validação externa ordenada: M1 K=3 sem associação significativa; M0 K=3/K=5 e A_so_dur têm 4–5 desfechos com p < 0,05. Com 49 segmentos e vários desfechos, isso não corrige para comparações múltiplas.

## Sensibilidade à construção
- Retirar uma dimensão (L_sem_*) ou usá-la sozinha (A_so_*) altera tamanhos e estabilidade, mas nenhuma variante isolada domina em todas as métricas.
- Variantes da escada de água e esgoto (E_*): partições K=3 parecidas (20/17/12 a 22/16/11); AMI com subpopulação varia bastante (0,015 a 0,231), o que mostra que essa métrica é sensível à definição dessas variáveis.
- CLR vs ILR: S_clr tem silhueta 0,257 e partição 23/18/8, ligeiramente diferente de ILR.

## Limites
- 49 segmentos dão pouco poder; ausência de significância não é prova de ausência de estrutura.
- Nenhum K=3 bate o nulo de covariância; a tipologia deve ser apresentada como descrição analítica, não como descoberta de grupos naturais.
- Pendências fora do pipeline: validação externa sob incerteza da partição, jackknife por domicílio, comparação com a Emater.
- Intervalos do modelo misto por variacional Bayes subestimam a incerteza; o modelo domiciliar não usa pesos.
