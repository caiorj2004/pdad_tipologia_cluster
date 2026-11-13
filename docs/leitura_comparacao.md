# Leitura da comparação de variantes (13/11/2026) — sem eleger modelo

Fonte: `resultados/comparacao_variantes.csv` (19 variantes). Todos os números abaixo vêm dessa tabela.
Bootstrap: 200 réplicas e 30 inicializações para M0_antigo_K5, M0_antigo_K3, M1_revisado_K3; 100 réplicas e 20 inicializações para as demais; as duas variantes `_b500n50` usam 500 e 50.

## 1. O que reproduz os diagnósticos da orientadora (modelo antigo, K=5)
| Diagnóstico dela | Reprodução aqui |
|---|---|
| p de covariância em K=5 = 0,110 | 0,110 |
| AMI dos rótulos com a RA 0,16 e com a subpopulação 0,06 | 0,159 e 0,061 |
| grupo 2 com 3 segmentos e 24 domicílios (0,5%) | menor grupo com 3 segmentos e 0,53% dos domicílios |
| Spearman tamanho × distância ao centro = −0,45 (4 variáveis, 49 segmentos) | −0,42 (matriz completa) |

## 2. O que vale em todas as variantes (não depende de qual modelo vença)
1. **Nenhuma variante com K=3 supera o nulo de covariância** (p de 0,14 a 0,88). Só aparecem K significativos em partes das variantes do tipo "uma mudança por vez" e no modelo antigo (K=7, p=0,016); o modelo revisado não tem nenhum (menor p = 0,13, em K=8).
2. **Segmentos pequenos são mais extremos em todas as variantes**: Spearman entre n e distância ao centro de −0,42 a −0,56, e Mann-Whitney (n < 20 mais distantes) com p < 0,001 em todas. Não foi resolvido por nenhuma mudança de matriz testada.
3. **Validação externa**: CadÚnico difere entre grupos (Kruskal-Wallis entre segmentos, p < 0,05) e o IVS-DF (fora da matriz nas variantes `ctx`) ordena os grupos (p < 0,01) em todas as variantes. A insegurança alimentar moderada ou grave difere (p < 0,05) em todas, exceto em `L_sem_agua` (p = 0,060). A escritura registrada raramente difere (p < 0,05 em 3 de 19).
4. **A winsorização da renda não altera nada**: `S_winsor95` e `S_winsor100` produzem a mesma partição de `M1_revisado_K3` (bootstrap 0,384 e 0,385 contra 0,381).
5. **Só K=2 é estável no bootstrap** (ARI médio de 0,58 a 0,73 em K=2 contra 0,37 a 0,57 para K=3 a 5).

## 2b. Onde os resultados divergem entre as variantes
| Métrica | Modelo antigo (K=5) | Revisado (K=3) | Observação |
|---|---|---|---|
| ARI médio do bootstrap | 0,46 (0,47 com 500 réplicas) | 0,38 | variantes intermediárias de 0,38 a 0,57 |
| Puxão da RA (pares da mesma RA no mesmo grupo, p) | 0,006 | 0,142 | recolocar IVS-DF e escala da RA no revisado (`L_sem_ctx`) leva a 0,005 e sobe o ARI do bootstrap para 0,515 |
| Partição sem uma subpopulação (ARI médio, PCA refeito) | 0,76 | 0,45 | no revisado, tirar "Características Urbanas" leva o ARI a 0,16 |
| Partição sem uma RA (ARI médio, PCA refeito) | 0,81 | 0,74 | RA crítica varia (Fercal, Sobradinho, Gama, Ceilândia) |
| Ótimo do K-Means com 10 inicializações | 34% das sementes | 18% | de 7% a 73% entre variantes |
| AMI com a RA / com a subpopulação | 0,159 / 0,061 | 0,071 / 0,072 | "mais associado à RA" vale para o antigo, não para o revisado |

## 3. Pontos que continuam em aberto (não decididos)
* **Qual modelo é melhor.** Em estabilidade e separação do nulo o antigo tem números melhores; em ausência de variáveis constantes na RA e de composições colineares o revisado é mais defensável. As métricas de estabilidade melhoram quando o IVS-DF e a escala da RA voltam à matriz (`L_sem_ctx`), mas junto com o puxão da RA (item 2b).
* **O revisado depende fortemente da subpopulação Características Urbanas** (S6b): isso pode ser um traço do território ou uma limitação da matriz.
* **Grupo pequeno do antigo (3 segmentos)** some no revisado, mas os segmentos pequenos continuam extremos (item 2.2).
* Testes ainda não feitos: refazer toda a bateria com o dicionário aplicado reprocessado; variante com mediana ponderada de renda (descartada por decisão do autor); variante deixando fora também os 6 segmentos não classificáveis.
