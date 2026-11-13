# Auditoria do pipeline contra a avaliação (13/11/2026)

Revisão do repositório comparando (i) o que a avaliação pede, incluindo o que vai além das 15 caixas do checklist, e (ii) o que o pipeline fazia. Nenhum modelo é eleito aqui.

## Desvios encontrados e o que foi feito

| # | Pedido/expectativa da avaliação | Situação encontrada | Correção |
|---|---|---|---|
| 1 | Bootstrap **dentro dos estratos** amostrais, B ≥ 500, 50 inicializações, alinhamento húngaro | O bootstrap reamostrava domicílios **dentro de cada segmento** (tamanho do segmento fixo); 200 réplicas × 30 inicializações nas variantes principais | `BOOT_NIVEL=estrato` no programa 16 (tamanho do segmento varia, redesenho se algum segmento ficar vazio). Variantes `*_bootEstratoRegra`; as variantes `*_b500n50` já cobrem B=500 e 50 inicializações |
| 2 | Estabilidade deve incluir a **incerteza na escolha do nº de componentes** | O PCA mantinha fixo o nº de componentes da referência em cada réplica | `BOOT_RETENCAO=regra`: refaz a análise paralela de Monte Carlo em cada réplica |
| 3 | Encolhimento com **n efetivo de Kish**, alvo = **média da subpopulação**, τ² por momentos | O programa 16 encolhia para a média global com n bruto | Seção (d2) do programa 16 (a seção (d) antiga foi mantida para comparação) |
| 4 | Dizer se a estabilidade medida é maior do que a esperada **sem estrutura** | Existia nulo para silhueta e para a subamostragem, mas **não** para o ARI do bootstrap | Programa 24, parte 1: nulo paramétrico (mesma covariância, mesmo ruído amostral por segmento) |
| 5 | Replicabilidade fora da amostra | Só subamostras de segmentos e bootstrap, ambos usando os mesmos domicílios | Programa 24, parte 2: **divisão em metades** de domicílios dentro de cada segmento, com nulo |
| 6 | "Homogeneidade territorial" não pode ser só efeito de composição/desenho | AMI era comparada a permutação de rótulos de **segmento**, que ignora o estrato | Programa 24, parte 3: redistribuição de domicílios entre os segmentos **do mesmo estrato** |
| 7 | Decomposição **multinível** (MixedLM), não só η² | Havia decomposição por ANOVA | Programa 25, parte 1: subpopulação fixa, RA aleatória, p de permutação dentro de subpopulação |
| 8 | Autocorrelação espacial sem confundir com mistura de subpopulações | Join count/Moran sem ajuste de composição | Programa 25, parte 2: Moran I bruto e ajustado pela média da subpopulação |
| 9 | Validação externa coerente com a ordem dos grupos | Só Kruskal-Wallis (não usa ordem) e Wilson | Programa 25, parte 3: Spearman entre a posição do grupo (ordenado por renda) e cada indicador, com permutação |
| 10 | Escada de água/esgoto (rede / individual adequado / precário) com critérios discutíveis | Uma única definição das classes | Variantes `E_*`: cisterna como precária, nascente/córrego como individual, "gambiarra" como rede, fossa rudimentar como individual (`AGUA_ESCADA`, `ESGOTO_ESCADA`) |
| 11 | K não deve ser decidido por conveniência | Havia K=3 e K=5 por modelo | Variantes `M0_antigo_K4` e `M1_revisado_K4` |
| 12 | Nome da dimensão D1 coerente com o conteúdo | O rótulo antigo (`D1_Financeiro_Ativos`) permanecia nas saídas | Renomeado para `D1_Bens_Duraveis` em todo o pipeline (as saídas foram regeradas) |
| 13 | Reprodutibilidade | Os scripts `rodar_variante.sh/rodar_todas.sh` estavam desatualizados e ignoravam o programa 23 | Passaram a chamar `reproduzir_tudo.py` |

## O que continua fora (e por quê)

- Mediana ponderada (a avaliação diz "winsorização **ou** mediana"; escolhida a winsorização).
- V de Cramér corrigido (substituído por AMI com permutação).
- Regressão logística em nível domiciliar com erro robusto por conglomerado e controles: **não implementada**; a validação externa continua no nível do segmento e é, portanto, descritiva.
- Validação externa sob incerteza da partição (rótulos de cada réplica do bootstrap) e comparação com Emater: pendentes.
- Teste de dependência em relação a um único domicílio/segmento influente (jackknife por domicílio): pendente.
- O ruído do nulo do bootstrap vem das réplicas do próprio bootstrap (preserva a correlação dentro do segmento), mas assume que a covariância "verdadeira" é a observada descontado o ruído; é uma aproximação.

## Limites a ter em mente ao ler as rodadas

1. Com 49 segmentos, qualquer métrica de partição tem pouco poder; o p de permutação não distingue "sem estrutura" de "estrutura fraca".
2. Variantes que mudam o conjunto de variáveis (M0 x M1) mudam também o espaço onde os nulos são calculados; compare as métricas **relativas ao nulo** (p, ARI observado − ARI nulo), não os valores brutos.
3. Todas as rodadas usam a mesma semente (42); a sensibilidade à semente é coberta apenas pela estatística de ótimo global do K-Means (programa 21).
