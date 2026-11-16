#!/usr/bin/env python3
"""Embute app/dados.json em app/template.html e grava docs/index.html (página única, sem dependências)."""
from pathlib import Path
R = Path(__file__).resolve().parent.parent
dados = (R/'app'/'dados.json').read_text(encoding='utf-8').replace('</', '<\\/')
html = (R/'app'/'template.html').read_text(encoding='utf-8').replace('/*DADOS*/', dados)
(R/'docs').mkdir(exist_ok=True)
(R/'docs'/'index.html').write_text(html, encoding='utf-8')
(R/'docs'/'.nojekyll').write_text('')
print('ok', round(len(html)/1024), 'KB')
