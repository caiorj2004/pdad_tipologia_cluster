#!/usr/bin/env bash
# Uso: pipeline/rodar_variante.sh <variante>  -> atalho para: python reproduzir_tudo.py --variantes <variante> --paralelo 1
cd "$(dirname "$0")/.." && python reproduzir_tudo.py --modo completo --variantes "$1" --paralelo 1
