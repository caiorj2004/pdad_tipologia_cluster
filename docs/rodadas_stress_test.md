# Rodadas extras de stress test: análise por modelo, K e parâmetro

Base: 27 variantes completas (etapas 03, 16, 18, 19, 21, 23, 24, 25, 26), geradas em 2026-11-13, sem falhas. Todos os números abaixo vêm de `docs/tabelas_rodadas.md` (tabelas T1 a T8, geradas por `pipeline/27_tabelas_rodadas.py` a partir de `resultados/comparacao_variantes.csv`). Este documento **descreve e pesa evidências; não elege modelo**. O texto do TCC não foi alterado.

## 0. Como ler

**O que cada variante é.** M0 é o modelo antigo (34 variáveis, nenhum ajuste da avaliação). M1 é o revisado (30 variáveis, cinco ajustes ligados). Os ajustes são os cinco flags de `variantes.tsv`:

| flag | o que muda na matriz |
|---|---|
| `ctx` | IVS-DF e número de estabelecimentos saem da matriz (viram contexto externo) |
| `agua` | água por escada de serviço (adequação + tipo), em ilr ou clr |
| `esgoto` | esgoto por escada de serviço, idem |
| `renda` | renda = média de ln(1+x) winsorizada, no lugar de ln da renda domiciliar per capita |
| `dur` | geladeira e secadora saem; D1 passa a "bens duráveis" |

- `A_so_X` = M0 **mais apenas** o ajuste X (K = 3).
- `L_sem_X` = M1 **menos apenas** o ajuste X (K = 3).
- `S_*` = M1 K=3 com winsor 0,95 / 1,00 ou composição clr. `E_*` = M1 K=3 com uma definição alternativa de escada (cisterna precária, nascente individual, "gambiarra" na rede, esgoto rudimentar individual).
- `*_b500n50`, `*_bootEstratoRegra` = mesmo modelo central, só muda o bootstrap (500 réplicas × 50 inicializações; reamostragem por estrato com retenção de componentes redecidida por réplica).

> **Correção.** A versão anterior de `leitura_comparacao.md` dizia que `A_so_*` e `L_sem_*` "retiram uma dimensão ou a usam sozinha". Isso estava errado: são ablações de **ajustes**, não de dimensões. O arquivo foi reescrito.

**Convenções das métricas** (ver tabelas):
- *p vs nulo de covariância*: fração de matrizes de 49 segmentos com a mesma covariância e **sem grupos** que atingem silhueta ≥ à observada. p alto significa que a silhueta observada é compatível com ausência de grupos.
- *Bootstrap*: ARI entre a partição de referência e a partição refeita em cada réplica de domicílios (fluxo inteiro refeito). *Nulo do bootstrap* (programa 24): mesmo procedimento sobre dados sintéticos sem grupos; "p vs nulo" = fração do nulo com ARI ≥ ao observado, em 100 matrizes nulas (resolução 0,01).
- *Fronteira*: segmento que fica no próprio grupo em menos de 60% das réplicas. *Não classificável*: n < 10 (6 segmentos, igual em todas as variantes).
- *ARI entre partições (T8)*: não depende da rotulação.

## 1. Achados que valem para quase todas as variantes

1. **Sem sinal de agrupamento contra o nulo de covariância no K final.** Das 27 linhas, nenhuma tem p_cov < 0,05 no K escolhido (T2). O menor p é 0,110 (M0 K=5) entre os modelos centrais. A estrutura que existe é de **contínuo com concentrações**, não de grupos separados.
2. **Estabilidade abaixo do nulo.** Em nenhuma variante o ARI do bootstrap observado supera o p95 do nulo (0 de 27, T3), nem o de metades (0 de 27). Para M1 K=3, o nulo atinge ARI ≥ ao observado em 96% das matrizes (p = 0,960); para M0 K=5, 39,6%. Leitura cautelosa: os grupos não são mais reprodutíveis do que seriam em dados sem grupos, **mas** o nulo gaussiano tem limites (seção 7).
3. **K=2 é sempre o mais estável.** ARI de bootstrap K=2 de 0,58 a 0,73 e silhueta ≈ 0,31 a 0,33, contra 0,37 a 0,57 em K=3. A passagem K=2→3 custa estabilidade em todas as variantes; mas p_cov em K=2 também é alto (≥ 0,158). K=2 não ganha respaldo estatístico adicional; só é mais fácil de reproduzir.
4. **Tamanho do segmento × distância ao centro.** Spearman negativo e p ≤ 0,02 em todas as 27 linhas (−0,33 a −0,63): segmentos com poucos domicílios amostrados ficam mais longe do centro, portanto mais "extremos". Isso é uma assinatura de ruído amostral na matriz e contamina a leitura dos grupos periféricos.
5. **Ótimos locais dominam sementes únicas.** Com uma só semente, só 0 a 10% das execuções atingem o ótimo (T4: % sementes no ótimo), e o ARI médio entre semente única e referência é 0,45 a 0,62. Os resultados centrais usam muitas inicializações; ainda assim, com 10 inicializações, só 7% a 73% das vezes se chega ao ótimo, conforme a variante.
6. **Decomposição de variância (dados, não partição).** No nível do domicílio: RA ≈ 5,7%, subpopulação ≈ 0,2%, segmento (RA × subpopulação) ≈ 12,0%, resíduo ≈ 82,1% (média de cinco indicadores: quatro binários e ln da renda; por indicador, o resíduo vai de 54% em água da rede a 97% em renda, e o componente de segmento chega a 46% em água da rede). Moran I bruto 0,003 e ajustado ≈ 0. Esses números são **idênticos em todas as variantes**, porque não dependem da partição. Descrevem o território; não discriminam modelos.
7. **Parâmetros que não mexem em nada.** Winsor 0,99, 1,00 e 0,95 geram a mesma partição de M1 K=3 (ARI = 1 entre as três; T8). Bootstrap de 200 vs 500 réplicas muda o ARI médio em ≤ 0,014 (M0 K=5: 0,458 → 0,472; M1 K=3: 0,381 → 0,381); reamostragem por estrato com retenção redecidida muda ≤ 0,008 no ARI (M0 K=5: 0,458 → 0,450; M1 K=3: 0,381 → 0,373) e muda a fronteira (M1 K=3: 21 → 17 segmentos).

## 2. Modelo antigo (M0, 34 variáveis)

M0 usa 6 componentes. Nas três linhas (K=3, 4, 5) a matriz e a ordenação dos segmentos são as mesmas; só muda o K. IVS-DF e estabelecimentos **estão na matriz**, portanto a associação com o IVS-DF em T7 (p ≈ 0) é parcialmente circular em M0 e nas variantes `A_so_*` sem `ctx`.

### M0 K=3
- Partição: [19, 17, 13] segmentos; menor grupo = 22,4% dos domicílios expandidos.
- Validade interna: silhueta 0,221; p_emb 0,002; **p_cov 0,635** (o pior de M0).
- Estabilidade: ARI subamostra 80% = 0,658; bootstrap 0,558 (p5 0,341); p vs nulo 0,347; metades 0,322 vs nulo p95 0,553; encolhimento 0,892 (2 segmentos mudam).
- Estrutura: sem uma RA (refaz PCA) 0,903, mínimo 0,378 (Ceilândia); 4% das RAs derrubam o ARI abaixo de 0,5. Fronteira: 7 segmentos (16,1% dos domicílios).
- Território: AMI com RA 0,112 (p perm. 0,007; p atrib. domiciliar 0,069); pares na mesma RA caem no mesmo grupo em 52,6% contra 32,6% esperado (p = 0,018).
- Validação: ordenada, 4 desfechos com p < 0,05, |rho| médio 0,320.
- **Sustenta:** a maior estabilidade de bootstrap e de deixar-uma-RA-fora entre os modelos centrais; grupos de tamanho balanceado.
- **Fragiliza:** silhueta compatível com ausência de grupos (p_cov 0,635); associação com RA acima do esperado pelo desenho (parte "território"); IVS-DF dentro da matriz.
- Equivalência: a partição é idêntica à de `A_so_renda` (ARI = 1): o ajuste de renda sozinho não altera a partição de M0 K=3.

### M0 K=4
- Partição [19, 18, 9, 3]; menor grupo = **0,53%** dos domicílios (3 segmentos).
- Silhueta 0,247 (a maior entre as linhas centrais de M0); p_cov 0,128; ARI80 0,645; bootstrap 0,503.
- **Encolhimento derruba o ARI para 0,404 (17 segmentos mudam)**; aparado 6% = 0,449. É a variante mais sensível ao encolhimento entre as de M0.
- AMI com subpopulação 0,007 (p atrib. 0,901): grupos ortogonais à subpopulação.
- **Sustenta:** silhueta mais alta de M0 e p_cov marginal (0,128).
- **Fragiliza:** grupo minúsculo, instável sob encolhimento.

### M0 K=5 (modelo avaliado)
- Partição [14, 14, 10, 8, 3]; um grupo com **0,53%** dos domicílios (3 segmentos).
- Silhueta 0,236; p_emb 0,004; p_cov 0,110; ARI80 0,664.
- Bootstrap 0,458 (p5 0,259; com 500 réplicas, 0,472); obs. vs nulo p95 0,551, p = 0,396. Metades 0,283 vs nulo p95 0,458 (praticamente igual).
- Estrutura: sem uma RA (refaz PCA) 0,805, **mínimo 0,507 e 0% das RAs com ARI < 0,5** (o melhor valor das 27 linhas nesta métrica, compartilhado com as réplicas b500n50 e bootEstratoRegra, que têm a mesma partição); sem uma subpopulação 0,760. Encolhimento 0,760 (5 mudam). Fronteira: 14 segmentos (26,9%).
- Território: AMI com RA 0,159 (p perm. 0,001; p atrib. 0,050); pares na mesma RA no mesmo grupo 44,7% contra 22,1% esperado (p = 0,006); join count 15 contra 9,3 esperado (p = 0,022): é o modelo com maior estrutura espacial/territorial entre os centrais.
- Validação: 4 desfechos ordenados com p < 0,05; Kruskal-Wallis para EBIA (insegurança qualquer p = 0,031; moderada/grave p = 0,036) e benefício social (p = 0,027).
- Em K=7 o p_cov é 0,016 (único K de M0 abaixo de 0,05; silhueta 0,251).
- **Sustenta:** melhor robustez a retirar uma RA; maior ancoragem territorial e de validação externa.
- **Fragiliza:** grupo de 0,5% dos domicílios (3 segmentos), que deixa a leitura desse grupo dependente de poucos segmentos; 14 segmentos na fronteira; p_cov 0,110; IVS-DF na matriz; o nulo de bootstrap mostra estabilidade compatível com sem-grupos.
- Observação: parte do sinal territorial pode ser mecânico: variáveis de contexto da RA (IVS-DF, estabelecimentos) estão na matriz e são constantes dentro da RA.

## 3. Modelo revisado (M1, 30 variáveis)

M1 usa 5 componentes (PC1 de D1, D2, D4 e dois PCs de D3). IVS-DF e estabelecimentos **fora** da matriz. A retenção de componentes varia entre réplicas (S6b: "nº de PCs varia" = sim em M1), o que não ocorre em M0.

### M1 K=3
- Partição [24, 17, 8]; menor grupo = 22,5% dos domicílios (8 segmentos).
- Silhueta **0,253** (a maior entre K finais centrais); p_emb 0,006; p_cov 0,305; ARI80 0,658.
- Bootstrap **0,381** (p5 0,197); nulo p95 0,687, p = 0,960; metades 0,263 vs nulo p95 0,687.
- Fronteira: **21 segmentos, 47,9% dos domicílios** (52,7% com 500 réplicas; 44,4% com bootstrap por estrato).
- Estrutura: sem uma RA (refaz PCA) 0,736, mínimo 0,314 (Sobradinho; RA crítica varia entre Gama e Sobradinho conforme o procedimento), 26% das RAs com ARI < 0,5. **Sem uma subpopulação 0,446, mínimo 0,159**: a partição depende fortemente de uma subpopulação. Corte n ≥ 20 → ARI 0,213.
- Território: AMI com RA 0,071 (p atrib. 0,683; p perm. 0,040); com subpopulação 0,072 (p atrib. 0,059; p perm. 0,052). Join count 20 vs 15,4 (p = 0,083). O sinal territorial é fraco e, contra o nulo de atribuição domiciliar, o AMI com a RA não se distingue do esperado pelo desenho amostral.
- Validação externa: **0 desfechos ordenados significativos**, |rho| médio 0,053; Kruskal-Wallis: EBIA moderada/grave p = 0,013, CadÚnico p = 0,001, IVS-DF p < 0,001 (aqui o IVS-DF é externo à matriz); benefício social p = 0,164; escritura p = 0,065.
- **Sustenta:** maior silhueta; IVS-DF externo e associado; grupos de tamanho razoável; ARI de subamostra 80% de 0,658, no nível de M0.
- **Fragiliza:** pior ARI de bootstrap entre os K=3 centrais (0,381 vs 0,558 em M0 K=3); fronteira em quase metade dos domicílios; dependência de subpopulação; validação ordenada sem sinal; p_cov 0,305.
- Equivalências: partições idênticas às de `S_winsor95`, `S_winsor100`, `*_b500n50` e `*_bootEstratoRegra` (ARI = 1).
- Parentes próximos: `L_sem_esgoto` (ARI 0,943), `L_sem_dur` (0,872), `E_esgoto_rudimentar_individual` (0,792), `S_clr` (0,792), `E_agua_cisterna_precaria` (0,771).

### M1 K=4
- Partição [14, 13, 11, 11]; menor grupo 18,0%.
- Silhueta 0,224; p_emb 0,016; p_cov 0,535; **ARI80 0,480 (o menor de todas as 27 linhas)**.
- Bootstrap 0,396 (p5 0,236); obs vs nulo p95 0,603, p = 0,802. Encolhimento 0,471 (17 mudam). Ward 0,485; aparado 6% 0,306.
- **Sem uma RA 0,457 (refaz PCA), mínimo 0,211, 70% das RAs derrubam o ARI abaixo de 0,5** (a pior das 27). Sem uma subpop 0,584.
- Território: AMI com subpopulação 0,152 (p atrib. 0,020, p perm. 0,002); pares na mesma subpopulação 32,5% contra 23,6% (p = 0,001): em K=4 os grupos se alinham à subpopulação (setor), não à RA.
- Ordenada: 2 desfechos significativos.
- **Sustenta:** grupos bem balanceados; associação com subpopulação.
- **Fragiliza:** pior estabilidade de subamostra e de retirada de RA; sensível a encolhimento.

### M1 K=5
- Partição [14, 12, 10, 9, 4]; menor grupo 8,1%.
- Silhueta 0,212; **p_emb 0,060** (único K final central que não separa do nulo embaralhado); p_cov 0,565; ARI80 0,566; bootstrap 0,368; encolhimento 0,428 (18 mudam).
- Sem uma RA 0,664, mínimo 0,451 (Sol Nascente/Pôr do Sol); GMM 0,580 e Ward 0,755 (concordância com outros algoritmos acima dos outros K de M1).
- Território: AMI com subpopulação 0,143 (p perm. 0,003); com RA 0,119 (p perm. 0,003; p atrib. 0,267); join count 11 vs 8,3 (p = 0,168).
- Arquétipos: variância explicada 0,859 (a maior entre as linhas de M1).
- **Sustenta:** concordância com GMM/Ward; separa subpopulação.
- **Fragiliza:** silhueta menor; p_emb não significativo; bootstrap mais baixo.

### M0 vs M1 (T8)
- ARI entre M0 K=5 e M1 K=3: **0,402**; M0 K=3 vs M1 K=3: 0,376. As partições dos dois modelos concordam pouco: a revisão **muda** a tipologia, não só a ajusta.
- Entre M1 K=3 e M1 K=4: 0,434; K=5: 0,442. Os K de M1 não são "refinamentos" da mesma partição, já que o ARI entre eles é moderado.

## 4. Ablações de ajustes (K = 3)

### A_so_X (M0 + um ajuste)
| variante | partição | sil. | p_cov | ARI80 | boot | sem RA | ordenada p<0,05 | observação |
|---|---|---|---|---|---|---|---|---|
| A_so_ctx | [20,17,12] | 0,231 | 0,571 | 0,667 | 0,534 | 0,872 | 5 | ARI vs M1 K=3 = 0,332 |
| A_so_agua | [19,16,14] | 0,197 | 0,882 | 0,657 | 0,539 | 0,733 | 3 | menor silhueta do conjunto; **ARI após encolhimento 1,000** |
| A_so_esgoto | [20,15,14] | 0,246 | 0,309 | 0,625 | 0,516 | 0,813 | 5 | idêntica a A_so_dur |
| A_so_renda | [19,17,13] | 0,218 | 0,657 | 0,665 | 0,552 | 0,874 | 4 | idêntica a M0 K=3 |
| A_so_dur | [20,15,14] | 0,263 | 0,138 | 0,706 | 0,565 | 0,808 | 5 | maior silhueta e ARI80 das A_so; 3 K com p_cov < 0,05 |

- Nenhum ajuste isolado aproxima M0 de M1: o ARI com a partição de M1 K=3 fica entre 0,33 e 0,52. O efeito de M0 → M1 é **combinado**.
- `dur` e `esgoto` produzem a mesma partição (ARI = 1) quando aplicados sozinhos a M0.
- `A_so_dur` tem a melhor silhueta (0,263) e p_cov (0,138) do bloco e **3 valores de K com p_cov < 0,05** (K=7 e outros). Seu p_cov no K final (0,138) fica no patamar de M0 K=4 (0,128) e M0 K=5 (0,110), sem cruzar 0,05.
- `A_so_agua` e `A_so_ctx` não passam do nulo de covariância (0,882 e 0,571).

### L_sem_X (M1 menos um ajuste)
| variante | partição | sil. | p_cov | boot | fronteira (% dom.) | sem RA | ordenada p<0,05 | ARI vs M1 K=3 |
|---|---|---|---|---|---|---|---|---|
| L_sem_ctx | [19,16,14] | 0,199 | 0,880 | 0,515 | 1,8 | 0,805 | 3 | 0,515 |
| L_sem_agua | [19,17,13] | 0,230 | 0,635 | 0,519 | 13,1 | 0,756 | 3 | 0,242 |
| L_sem_esgoto | [23,17,9] | 0,256 | 0,262 | 0,410 | 36,1 | 0,751 | 0 | 0,943 |
| L_sem_renda | [20,16,13] | 0,219 | 0,792 | 0,428 | 20,7 | 0,670 | 0 | 0,684 |
| L_sem_dur | [22,18,9] | 0,244 | 0,337 | 0,380 | 40,8 | 0,867 | 0 | 0,872 |

- Retirar `ctx` ou `agua` de M1 muda muito a partição (ARI vs M1 K=3 de 0,515 e 0,242) e **eleva** a estabilidade de bootstrap (0,515 e 0,519 contra 0,381) e reduz a fronteira (1,8% e 13,1% dos domicílios contra 47,9%). Ou seja: a instabilidade de M1 K=3 está ligada à presença desses dois ajustes.
- Retirar `esgoto`, `renda` ou `dur` pouco altera a partição (ARI 0,684 a 0,943) e mantém a instabilidade (bootstrap 0,38 a 0,43).
- `L_sem_ctx` = `A_so_agua` (mesma partição, ARI = 1) e `L_sem_renda` = `E_agua_nascente_individual`.
- Em L_sem_ctx, o IVS-DF volta para a matriz: a associação com IVS (p ≈ 0) deixa de ser externa.
- A estabilidade maior das variantes sem `ctx`/`agua` **não é evidência de partição melhor**: p_cov é alto nelas (0,880 e 0,635), isto é, a silhueta é compatível com ausência de grupos.

## 5. Parâmetros de construção (M1 K=3)

| variante | mudança | partição | ARI vs M1 K=3 | boot | p_cov | observação |
|---|---|---|---|---|---|---|
| S_winsor95 | winsor 0,95 | [24,17,8] | 1,000 | 0,384 | 0,307 | sem efeito |
| S_winsor100 | sem winsor | [24,17,8] | 1,000 | 0,385 | 0,305 | sem efeito |
| S_clr | clr em vez de ilr | [23,18,8] | 0,792 | 0,452 | 0,260 | silhueta 0,257; IVS-DF p = 0,003; **associação ordenada: 0 desfechos** |
| E_agua_cisterna_precaria | cisterna precária na escada | [22,16,11] | 0,771 | 0,434 | 0,685 | AMI subpop 0,155 |
| E_agua_nascente_individual | nascente como individual | [20,16,13] | 0,684 | 0,446 | 0,709 | AMI subpop 0,231 (p atrib. 0,010) |
| E_agua_gambiarra_na_rede | gambiarra como rede | [20,17,12] | 0,332 | 0,506 | 0,697 | **ARI após encolhimento 0,316 (18 mudam)**; AMI subpop 0,015; 4 desfechos ordenados |
| E_esgoto_rudimentar_individual | rudimentar como individual | [20,17,12] | 0,792 | 0,398 | 0,836 | menor AMI com RA (0,050) |

- **Winsorização não importa** para M1 K=3. **A escolha ilr vs clr importa moderadamente** (ARI 0,792 entre as partições).
- A **classificação de uma fonte de água em "rede", "individual" ou "precária" muda a partição** (ARI 0,33 a 0,77 vs M1 K=3). O que mais muda é tratar "gambiarra" como rede (ARI 0,332). Isso é sensibilidade à codificação, não instabilidade amostral, e merece decisão metodológica explícita no texto.
- Nenhum desses K=3 supera o nulo de covariância (p_cov de 0,260 a 0,836; o menor é S_clr).

## 6. Bootstrap e retenção de componentes

| comparação | resultado |
|---|---|
| M0 K=5: 200×30 → 500×50 | ARI médio 0,458 → 0,472; nulo p95 0,551 → 0,557; fronteira 14 → 11 seg. |
| M1 K=3: 200×30 → 500×50 | ARI médio 0,381 → 0,381; fronteira 21 → 25 seg. (47,9% → 52,7% dos domicílios) |
| M0 K=5: domicílios → estrato + regra | ARI 0,458 → 0,450; fronteira 14 → 16 seg. |
| M1 K=3: domicílios → estrato + regra | ARI 0,381 → 0,373; fronteira 21 → 17 seg. |

O número de réplicas **converge**: o ARI médio e o p5 mudam pouco. A lista de segmentos de fronteira não é estável entre procedimentos (diferença de 4 segmentos em M1 K=3), então ela deve ser tratada como indicador aproximado.

## 7. Limites do pipeline (afetam todas as comparações)

1. **49 segmentos, 21 estratos.** Pouco poder estatístico; p alto não prova ausência de estrutura.
2. **Nulo gaussiano.** Gera matrizes com a covariância observada e ruído amostral extraído do bootstrap; não reproduz assimetria, caudas e valores nos limites das variáveis (proporções). O ARI observado e o do nulo vêm de rotinas diferentes (observado: realinhamento húngaro com NINIT_BOOT inicializações; nulo: PCA com número fixo de componentes e 10 inicializações). A comparação é indicativa. Com N_NULO = 100 a resolução do p é 0,01.
3. **Comparações múltiplas.** Há 7 desfechos de validação externa mais a ordenada, por variante, sem correção. Os p de Kruskal-Wallis tratam segmentos como independentes, mas o IVS-DF é constante dentro da RA (23 RAs): o p ≈ 0 do IVS-DF superestima a evidência.
4. **Circularidade parcial.** Em M0 e nas variantes sem `ctx`, IVS-DF e estabelecimentos estão na matriz. Em todas, renda e posse de bens estão na matriz e a associação com CadÚnico e benefícios é, em parte, mecânica.
5. **Modelo misto por variacional Bayes** subestima a incerteza; o modelo domiciliar não usa pesos amostrais.
6. **Deixar-uma-RA-fora** refaz o PCA (S6b) mas não refaz o redimensionamento do K; o número de componentes varia em M1.
7. **Não realizados:** validação externa sob incerteza da partição; jackknife por domicílio; comparação com técnicos da Emater-DF; mediana ponderada de renda.

## 7b. Validação externa sob incerteza da partição (programa 28)

O programa 19 testa os sete desfechos externos contra a partição de referência. O programa 28 mantém o desfecho de cada segmento fixo e refaz o Kruskal-Wallis com o rótulo de cada réplica do bootstrap (200 réplicas, já realinhadas). Pergunta: a associação sobrevive quando a partição é refeita com a incerteza amostral?

Percentual das réplicas com p < 0,05 (mediana do p entre parênteses):

| desfecho | M0 K=3 | M0 K=4 | M0 K=5 | M1 K=3 | M1 K=4 | M1 K=5 |
|---|---|---|---|---|---|---|
| EBIA insegurança qualquer | 32,5 (0,095) | 18,5 (0,136) | 14,5 (0,151) | 8,0 (0,227) | 6,5 (0,265) | 5,0 (0,331) |
| EBIA moderada/grave | 75,0 (0,022) | 62,5 (0,032) | 51,5 (0,047) | 53,0 (0,043) | 42,5 (0,066) | 30,0 (0,082) |
| benefício social | 84,0 (0,015) | 72,5 (0,026) | 63,5 (0,037) | 38,5 (0,081) | 30,5 (0,091) | 25,0 (0,092) |
| CadÚnico | 100 (<0,001) | 98,5 (0,001) | 100 (0,001) | 93,5 (0,004) | 94,5 (0,003) | 97,0 (0,003) |
| escritura registrada | 22,5 (0,106) | 28,5 (0,090) | 20,5 (0,106) | 26,0 (0,104) | 15,0 (0,129) | 9,0 (0,167) |
| documentação formal | 12,0 (0,204) | 7,5 (0,228) | 3,5 (0,268) | 4,5 (0,305) | 3,0 (0,319) | 3,0 (0,418) |
| outorga de água | 21,0 (0,102) | 15,5 (0,123) | 12,5 (0,136) | 19,5 (0,120) | 9,5 (0,172) | 10,0 (0,197) |

Leitura:
- **CadÚnico** é robusto em todos os modelos (93,5% a 100% das réplicas). Mas renda está na matriz de todos eles, e o CadÚnico é focalizado por renda: a associação é em parte mecânica.
- **EBIA moderada/grave e benefício social** sobrevivem na maior parte das réplicas em M0 (51,5% a 84%) e menos em M1 (25% a 53%).
- **EBIA qualquer, escritura, documentação formal e outorga** não sobrevivem em nenhum modelo (≤ 32%).
- Dentro de M1, a robustez cai de K=3 para K=5 em EBIA moderada/grave (53,0% → 30,0%) e em benefício social (38,5% → 25,0%).
- O epsilon² mediano médio (tamanho de efeito) é 0,11 em M0 e 0,07 em M1: efeitos pequenos nos dois.
- Entre as 27 linhas, o número de desfechos com p mediano < 0,05 vai de 1 a 3 (T9 em `tabelas_rodadas.md`). Esta rodada favorece M0 sobre M1 em validação externa, mas o contraste é modesto e sem correção para múltiplos testes.

Limites: o desfecho do segmento não é reamostrado; só a partição varia (a incerteza do desfecho está no IC do programa 19). 49 segmentos; sem correção para comparações múltiplas.

## 8. Quadro de trade-offs (sem ranking)

| critério | favorece | desfavorece |
|---|---|---|
| silhueta no K final | M1 K=3 (0,253), A_so_dur (0,263), S_clr (0,257) | M1 K=5 (0,212), A_so_agua (0,197), L_sem_ctx (0,199) |
| p vs nulo de covariância no K final | M0 K=5 (0,110), A_so_dur (0,138), M0 K=4 (0,128) | A_so_agua (0,882), L_sem_ctx (0,880), E_esgoto_rud. (0,836) |
| bootstrap (ARI médio) | M0 K=3 (0,558), A_so_dur (0,565) | M1 K=3, K=4, K=5 (0,38 a 0,40), L_sem_dur (0,380) |
| fração de domicílios na fronteira | A_so_agua (0,8%), L_sem_ctx (1,8%) | M1 K=3 (47,9%; 52,7% com 500 réplicas) |
| menor grupo (% dom.) | M0 K=3 e M1 K=3 (22%) | M0 K=4/K=5 (0,53%) |
| sensibilidade a encolhimento | A_so_agua (1,000), M0 K=3 (0,892) | M0 K=4 (0,404), E_agua_gambiarra (0,316), L_sem_agua (0,300) |
| retirar uma RA | M0 K=5 (0% RAs com ARI < 0,5), M0 K=3 (4%) | M1 K=4 (70%), L_sem_esgoto (35%) |
| sinal territorial (RA) | M0 K=5, K=4, A_so_esgoto (p perm. ≤ 0,003) | M1 K=3 (p atrib. 0,683), S_clr (0,752) |
| validação externa sob incerteza da partição (7b) | M0 K=3/K=5 (3 desfechos com p mediano < 0,05), A_so_* (3) | M1 K=4/K=5, S_clr (1) |
| validação externa ordenada | A_so_ctx, A_so_esgoto, A_so_dur (5 desfechos) | M1 K=3, S_*, L_sem_{esgoto,renda,dur}, E_* exceto gambiarra (0) |
| IVS-DF como validação **externa** | M1 e derivadas (p ≤ 0,003) | M0 e A_so_* sem ctx (circular) |

Nenhuma variante domina as demais em todos os critérios. M0 K=5 e M0 K=3 vão melhor em estabilidade e validação, mas têm o IVS-DF dentro da matriz e, no caso de K=5, um grupo de 0,5% dos domicílios. M1 K=3 vai melhor em silhueta e em desenho (contexto fora da matriz), mas é o menos estável entre os K=3. O trade-off entre **desenho metodologicamente correto (M1)** e **estabilidade/validação (M0)** não se resolve com os dados desta rodada.

## 9. O que se pode e o que não se pode afirmar

**Pode-se afirmar, com os números desta rodada:**
- Nenhum K final bate o nulo de covariância; a tipologia é uma partição descritiva de um contínuo.
- A partição muda substancialmente entre M0 e M1 (ARI 0,40) e entre codificações alternativas de água e esgoto (ARI 0,33 a 0,79).
- Winsor, número de réplicas e nível do bootstrap não alteram as conclusões.
- Segmentos com menos domicílios tendem a ficar nas pontas (efeito de tamanho).
- A maior parte da variância dos indicadores está no nível do domicílio (≈ 82% residual).

**Não se pode afirmar:**
- Que algum dos modelos é "o melhor" (trade-offs acima).
- Que existam "grupos naturais".
- Que o IVS-DF validou os grupos em M0 (circular).
- Que a validação externa é favorável ou desfavorável a M1 em geral: ela depende do modelo e do desfecho, sem correção para múltiplos testes.

## 10. Rodadas ainda abertas
- Jackknife por domicílio e por estrato.
- Comparação com a Emater-DF (pesquisa futura).
- Escolha de K por critério pré-registrado, antes de olhar os resultados externos.
