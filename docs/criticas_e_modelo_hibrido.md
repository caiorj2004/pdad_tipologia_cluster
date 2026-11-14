# As críticas pelo mérito próprio e o modelo híbrido M2

Números: `docs/tabelas_rodadas.md` (T1 a T10) e `resultados/comparacao_variantes.csv`. Este documento foi escrito antes da escolha. **Depois dele, o autor adotou M2 K=3 como modelo principal do texto do TCC**, com M0 e M1 como comparação e o M2b como terceira sensibilidade (seção 4). A análise abaixo não foi alterada.

## 1. Cada crítica avaliada por si, e não por ter vindo da avaliação

Para cada ajuste de M1 há três perguntas: o argumento se sustenta conceitualmente? Os dados desta pesquisa o confirmam? O ajuste muda a partição? Para a última pergunta, o ARI entre a partição com o ajuste e sem ele, a partir de M0 e a partir de M1 (T8):

| ajuste | M0 + ajuste vs M0 K=3 | M1 − ajuste vs M1 K=3 |
|---|---|---|
| `ctx` (IVS-DF e escala da RA fora) | 0,879 | 0,515 |
| `agua` (escada) | 0,843 | **0,242** |
| `esgoto` (escada) | 0,740 | 0,943 |
| `renda` (ln da média winsorizada) | 1,000 | 0,684 |
| `dur` (geladeira e secadora fora) | 0,740 | 0,872 |

A partir de M0, nenhum ajuste isolado muda muito a partição (0,74 a 1,00), mas os cinco juntos produzem ARI de 0,376 (M0 K=3 vs M1 K=3). O efeito é **não aditivo**: o que mais pesa em M1 é a água por escada (sem ela, ARI de 0,242 com o próprio M1 K=3), depois o contexto fora da matriz.

### ctx: IVS-DF e escala da RA fora da matriz
- **Argumento:** as duas variáveis são constantes dentro da RA; tornam segmentos da mesma RA parecidos por construção. O IVS-DF vem da PDAD urbana. Conceitualmente sólido, e independente do resultado: uma variável constante dentro do grupo não descreve o domicílio.
- **Os dados:** o eta² da RA nos escores, ajustado pelos graus de liberdade (T10), é 0,50 em M0 e 0,415 em M1. A diferença existe (0,085) mas é modesta: a maior parte da "RA" nos escores permanece sem o contexto, porque domicílios da mesma RA são parecidos em outras variáveis (infraestrutura, conectividade).
- **O que muda:** partição moderadamente (0,515 a partir de M1).
- **Conclusão:** sólida e pouco consequente nos resultados. Custo de manter: nenhum. Mantém IVS-DF como validação externa (p < 0,001).

### agua e esgoto: escada de serviço (PNSR/JMP)
Há três críticas misturadas:
1. **Fechamento composicional** (partes de um todo geram correlações espúrias no PCA). Matematicamente correta (Aitchison, 1986). Com três partes e PCA por dimensão, a distorção é limitada, mas existe.
2. **Direção normativa** ("poço = carência"). Correta como crítica à **interpretação** (nomes de grupos, recomendações). No agrupamento, a distância é simétrica e não carrega o sinal.
3. **Prescrever a escada como codificação.** Esta é a parte discutível: a escada embute no agrupamento um julgamento sobre o que é "adequado". Os dados mostram o custo: codificações defensáveis diferentes mudam a partição (ARI de 0,33 a 0,79 contra M1 K=3; T8 e seção 5 de `rodadas_stress_test.md`), principalmente "gambiarra na rede" (0,332).
- **Os dados:** a água é onde M1 é mais sensível (0,242) e onde nasce a instabilidade de M1: sem a escada de água, M1 K=3 passa de 47,9% para 13,1% dos domicílios na fronteira e de ARI de bootstrap 0,381 para 0,519.
- **Conclusão:** a crítica 1 e a 2 se sustentam; a 3 é uma escolha de medida com forte efeito, não uma correção neutra. O checklist da avaliação aceita "escada **ou** log-razão", o que deixa espaço para uma log-razão das fontes originais, sem escada (ver seção 4).

### renda: ln(1+x) por domicílio, winsorizada
- **Argumento:** ln da média é sensível a valores extremos; média de ln com winsorização protege. Sólido.
- **Os dados:** sem efeito em M0 (ARI 1,000), efeito moderado em M1 (0,684). Winsor 0,95, 0,99 e 1,00 dão a mesma partição.
- **Conclusão:** sólida, barata, manter.

### dur: geladeira e secadora fora
- **Argumento:** quase sem variância (geladeira ≈ 99%), só ruído. Sólido; D1 é "bens duráveis", não "capital financeiro".
- **Os dados:** efeito pequeno em M1 (0,872).
- **Conclusão:** sólida, pouco consequente, manter.

### K e o grupo de 0,5%
- A avaliação diz que K=5 é K=4 mais um grupo de ruído. Os dados confirmam: em M0, K=4 e K=5 têm um grupo de 0,53% dos domicílios (3 segmentos), e o encolhimento derruba o ARI em M0 K=4 (0,404). Em K=3 (M0, M1, M2) o menor grupo tem 22% dos domicílios.
- **Conclusão:** a crítica se sustenta nos dados.

### Resumo do mérito
| ajuste | fundamento conceitual | confirmado nos dados | consequência na partição |
|---|---|---|---|
| ctx | forte | em parte (eta² RA −0,085) | moderada |
| água/esgoto: fechamento e leitura | forte | n/a | n/a |
| água/esgoto: escada como codificação | discutível (julgamento normativo na medida) | sensível à codificação | **grande** |
| renda | forte | sem efeito em M0 | pequena a moderada |
| dur | forte | sim | pequena |
| K=5 = K=4 + ruído | forte | sim | n/a |

## 2. Antes de olhar o M0: a vantagem de validação dele vem do contexto da RA?

Se M0 valida melhor porque o IVS-DF e a escala da RA estão na matriz, a vantagem seria mecânica. O programa 29 testa isso: o Kruskal-Wallis dos desfechos entre grupos é refeito com permutação dos rótulos **apenas entre segmentos da mesma RA** (5000 sorteios).

| modelo | desfechos com p < 0,05 dentro da RA | EBIA mod./grave | benefício social | CadÚnico |
|---|---|---|---|---|
| M0 K=3 | 2 | 0,005 | 0,086 | 0,011 |
| M0 K=4 | 3 | 0,006 | 0,064 | 0,008 |
| M0 K=5 | 1 | 0,126 | 0,066 | 0,001 |
| M1 K=3 | 0 | 0,072 | 0,417 | 0,067 |
| M1 K=4 | 1 | 0,244 | 0,186 | 0,006 |
| M1 K=5 | 2 | 0,039 | 0,370 | 0,004 |

A vantagem de M0 K=3 e K=4 **sobrevive** à condição de comparar dentro da mesma RA; a de M1 K=3 desaparece (nenhum desfecho). Logo o contexto da RA **não explica** a diferença de validação. Quem carrega a vantagem é a partição, com a escada de água, e não o contexto. (O CadÚnico é mecânico em parte, pois a renda está em todas as matrizes.)

## 3. Modelo híbrido M2

**Definição (por mérito, não por desempenho).** M2 aplica os ajustes cujo fundamento é inequívoco e a consequência não é arbitrária: `ctx`, `renda`, `dur`. Mantém água e esgoto como em M0 (fontes originais e variável de solução sanitária), em vez da escada, que é uma escolha de medida com grande efeito. Variantes: `M2_semEscadas_K3`, `_K4`, `_K5`, com a mesma bateria e o mesmo bootstrap (200 × 30) dos modelos centrais.

**Aviso sobre o processo.** Defini M2 **depois** de ver que a escada de água é o que mais separa M1 de M0 (`L_sem_agua`). É um caminho de decisão condicionado aos dados. A defesa é que a regra de inclusão é conceitual (quais ajustes são inequívocos), não um ajuste por métrica; mesmo assim, o leitor deve saber.

### Resultados (T1 a T10)
| | M0 K=3 | M0 K=5 | M1 K=3 | **M2 K=3** | M2 K=4 | M2 K=5 |
|---|---|---|---|---|---|---|
| partição (segmentos) | 19/17/13 | 14/14/10/8/3 | 24/17/8 | 20/17/12 | 16/13/12/8 | 15/12/10/7/5 |
| menor grupo (% dom.) | 22,4 | 0,5 | 22,5 | 22,2 | 19,7 | 2,9 |
| silhueta | 0,221 | 0,236 | 0,253 | 0,236 | 0,219 | 0,222 |
| p (nulo de covariância) | 0,635 | 0,110 | 0,305 | 0,535 | 0,623 | 0,437 |
| ARI subamostra 80% | 0,658 | 0,664 | 0,658 | **0,696** | 0,605 | 0,583 |
| ARI bootstrap | 0,558 | 0,458 | 0,381 | 0,535 | 0,407 | 0,390 |
| bootstrap vs nulo (p) | 0,347 | 0,396 | 0,960 | 0,337 | 0,713 | 0,733 |
| % dom. na fronteira | 16,1 | 26,9 | 47,9 | **3,2** | 37,4 | 54,5 |
| ARI após encolhimento | 0,892 | 0,760 | 0,794 | 0,772 | 0,669 | 0,364 |
| sem uma RA (refaz PCA) | 0,903 | 0,805 | 0,736 | 0,902 | 0,817 | 0,754 |
| RAs com ARI < 0,5 | 4% | 0% | 26% | 4% | 4% | 13% |
| sem uma subpopulação | 0,760 | 0,760 | 0,446 | 0,803 | 0,512 | 0,611 |
| Spearman n × distância | −0,43 | −0,42 | −0,55 | −0,52 | −0,57 | −0,40 |
| AMI com RA (p atrib.) | 0,112 (0,069) | 0,159 (0,050) | 0,071 (0,683) | 0,090 (0,139) | 0,114 (0,168) | 0,159 (0,010) |
| desfechos p mediano < 0,05, sob incerteza da partição (7 testados) | 3 | 3 | 2 | 2 | 2 | 1 |
| desfechos p < 0,05 dentro da RA | 2 | 1 | 0 | 2 | 1 | 1 |
| validação ordenada (p < 0,05) | 4 | 4 | 0 | 5 | 2 | 5 |
| IVS-DF externo à matriz | não | não | sim | sim | sim | sim |

### O que M2 K=3 mostra
- **Não é um "meio-termo" de M0 e M1; é quase M0 K=3.** ARI com M0 K=3 = 0,879; com M1 K=3 = 0,332. Isso decorre de a escada de água ser o que separa M1 de M0.
- **Estabilidade de M0 mantida** (bootstrap 0,535, fronteira 3,2%, subamostra 0,696, deixar uma RA fora 0,902) **com** o contexto fora da matriz, renda e bens duráveis ajustados.
- **Validação externa de M0 mantida**: 2 desfechos sobrevivem dentro da RA (EBIA mod./grave p = 0,033; CadÚnico p = 0,015), e 5 desfechos ordenados significativos; IVS-DF agora externo (p < 0,001).
- **O que não melhora:** a silhueta é compatível com ausência de grupos (p = 0,535, pior que M0 K=5 e M1 K=3); o ARI de bootstrap continua dentro do que o nulo gera (p = 0,337); o efeito tamanho × distância persiste (−0,52); o sinal territorial fica fraco (AMI com RA p = 0,139).
- **K=4 e K=5 em M2** repetem o padrão de M0 e M1: perdem estabilidade (fronteira de 37% e 55%), e K=5 tem um grupo de 2,9% dos domicílios.

### O que M2 não resolve
- Não atende o item 2 do checklist (escada de água e esgoto **ou** log-razão) **se** a escada for lida como obrigatória. A justificativa teria de ser: o fechamento composicional é a crítica correta e a escada é uma escolha de medida com forte efeito; mostrar a sensibilidade (M1, `E_*`) em vez de adotá-la.
- Dá o mesmo resultado nulo de grupos que M0 e M1 contra o nulo de covariância.
- Foi definido depois dos dados (aviso acima).

## 4. Alternativa M2b (não executada)
Satisfaz o checklist ao pé da letra (log-razão) sem a escada normativa: clr/ilr das três fontes originais (rede, poço artesiano, poço/cisterna), com tratamento de zeros, mais `ctx`, `renda`, `dur`. Exige uma codificação nova em `03` (as variantes `COMP=clr` atuais já agrupam em níveis da escada). É a forma mais defensável de responder à crítica 1 sem aceitar a crítica 3. Não rodei porque envolve uma decisão sua (tratamento de zeros e definição do balanço).

## 5. Como decidir sem circularidade
1. Declarar o critério antes de escolher: (a) desenho metodológico, (b) estabilidade, (c) validação externa.
2. M1 é o único que cumpre (a) à risca, mas é o pior em (b) e (c).
3. M0 é o melhor em (b) e (c), mas viola (a) no IVS-DF e na escala da RA.
4. M2 K=3 cumpre (a) em três dos quatro ajustes, mantém (b) e (c) de M0, e expõe a escada como sensibilidade.
5. Qualquer escolha deve reportar M1 e M0 como análises de sensibilidade, e declarar que a estrutura é um contínuo sem separação contra o nulo de covariância.
