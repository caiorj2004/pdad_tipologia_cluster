# Leitura da comparação de variantes (sem eleger vencedor)

Resumo de uma página. A análise detalhada por modelo, K e parâmetro está em `docs/rodadas_stress_test.md`; todas as tabelas, em `docs/tabelas_rodadas.md` (geradas por `pipeline/27_tabelas_rodadas.py`). Base: 27 variantes completas em `resultados/comparacao_variantes.{csv,md}`.

**Nomes.** M0 = modelo antigo; M1 = revisado (cinco ajustes: `ctx`, `agua`, `esgoto`, `renda`, `dur`). `A_so_X` = M0 mais só o ajuste X; `L_sem_X` = M1 menos só o ajuste X (são ablações de **ajustes**, não de dimensões). `S_*`, `E_*` e `*_b500n50`/`*_bootEstratoRegra` variam parâmetros do M1 ou M0 centrais.

## O que vale para quase todas as variantes
- Nenhum K final bate o nulo de covariância com p < 0,05 (menor p entre os centrais: M0 K=5, 0,110).
- O ARI de bootstrap e o de metades não superam o p95 do nulo em nenhuma das 27 variantes.
- K=2 é o mais estável em todas, e K=3 custa estabilidade.
- Spearman entre tamanho do segmento e distância ao centro é negativo em todas (−0,33 a −0,63; p ≤ 0,02).
- Decomposição domiciliar (≈ 82% residual, 12% segmento, 5,7% RA, 0,2% subpopulação) e Moran ajustado (≈ 0) não dependem da partição.
- Winsor (0,95, 0,99, 1,00) e número de réplicas do bootstrap (200 vs 500) não mudam as conclusões.

## O que muda entre modelos
- M0 → M1 muda a partição (ARI entre M0 K=5 e M1 K=3 = 0,402). Nenhum ajuste isolado reproduz o efeito combinado (ARI entre 0,33 e 0,52).
- M0 tem melhor bootstrap, estabilidade ao retirar uma RA e validação ordenada, mas tem IVS-DF dentro da matriz e, em K=4 e K=5, um grupo de 0,5% dos domicílios.
- M1 K=3 tem a maior silhueta entre os K finais centrais (0,253) e o contexto fora da matriz, mas o menor bootstrap entre os K=3 (0,381) e 47,9% dos domicílios na fronteira.
- Retirar `ctx` ou `agua` de M1 eleva muito a estabilidade, mas o p de covariância sobe (0,880 e 0,635).
- A codificação de água e esgoto (ilr vs clr; "gambiarra" como rede ou não; cisterna, nascente) muda a partição (ARI de 0,33 a 0,79 contra M1 K=3).

Nenhum modelo domina em todos os critérios. A decisão de modelo continua em aberto.
