#!/usr/bin/env bash
# Uso: pipeline/rodar_variante.sh <nome_da_variante>     (lê a linha correspondente em pipeline/variantes.tsv)
# Roda, dentro de resultados/<nome>/: 03 (modelo) -> 16 (bootstrap, nulos, encolhimento) -> 18 (variância/AMI/espacial) -> 19 (validação externa) -> 21 (testes de estresse)
set -euo pipefail
NOME="$1"; RAIZ="$(cd "$(dirname "$0")/.." && pwd)"
LINHA="$(grep -v '^#' "$RAIZ/pipeline/variantes.tsv" | awk -F'\t' -v n="$NOME" '$1==n')"
[ -n "$LINHA" ] || { echo "variante '$NOME' não encontrada"; exit 1; }
IFS=$'\t' read -r _ V2 K WIN COMP NB NI FIG <<< "$LINHA"
[ "$V2" = "nenhuma" ] && V2=""
OUT="$RAIZ/resultados/$NOME"; rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
export TCC_BASE="$RAIZ/dados/gold" TCC_SHP="$RAIZ/dados/regioes_administrativas.zip" TCC_03="$RAIZ/pipeline/03_analise_cluster_4d.py"
export V2="$V2" K_FINAL="$K" WINSOR_P="$WIN" COMP="$COMP" N_BOOT="$NB" NINIT_BOOT="$NI" MPLBACKEND=Agg PYTHONHASHSEED=0
printf 'variante=%s V2=%s K=%s WINSOR_P=%s COMP=%s N_BOOT=%s NINIT_BOOT=%s\n' "$NOME" "$V2" "$K" "$WIN" "$COMP" "$NB" "$NI" > config.txt
python "$TCC_03"                                  > saida_03.txt 2>&1
python "$RAIZ/pipeline/16_robustez_amostral_e_pesos.py" > saida_16.txt 2>&1
python "$RAIZ/pipeline/18_variancia_ami_espacial.py"    > saida_18.txt 2>&1
python "$RAIZ/pipeline/19_validacao_externa.py"         > saida_19.txt 2>&1
python "$RAIZ/pipeline/21_stress_estrutural.py"         > saida_21.txt 2>&1
rm -rf shape_df_extracao
if [ "$FIG" != "1" ]; then find . -name '*.png' -delete; fi
echo "OK $NOME"
