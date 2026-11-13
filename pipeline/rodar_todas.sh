#!/usr/bin/env bash
# Roda todas as variantes de pipeline/variantes.tsv (atalho para: python reproduzir_tudo.py --modo completo)
cd "$(dirname "$0")/.." && python reproduzir_tudo.py --modo completo "$@"
