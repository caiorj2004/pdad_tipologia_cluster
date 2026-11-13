#!/usr/bin/env bash
# Roda todas as variantes de pipeline/variantes.tsv, N em paralelo (padrão 2). Uso: pipeline/rodar_todas.sh [N_PARALELO]
RAIZ="$(cd "$(dirname "$0")/.." && pwd)"; P="${1:-2}"
grep -v '^#' "$RAIZ/pipeline/variantes.tsv" | cut -f1 | xargs -P "$P" -I{} bash -c "'$RAIZ/pipeline/rodar_variante.sh' {} > '$RAIZ/resultados/log_{}.txt' 2>&1 || echo FALHOU {}"
