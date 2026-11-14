# -*- coding: utf-8 -*-
# Figuras: perfis dos grupos com IC 95% (bootstrap) e validação externa. Entrada: perfis_grupos_ic.csv, validacao_externa.csv
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os
NOM = ['A — Urbanização formal\nconsolidada', 'B — Saneamento e\nvias precários', 'C — Urbanização\nparcial'] if os.environ.get('GRUPOS_NOMES') == 'm2' else ['A — Urbanização formal\nconsolidada', 'B — Urbanização\nparcial', 'C — Saneamento e\nvias precários']; COR = ['#440154', '#21918c', '#fde725']
p = pd.read_csv('perfis_grupos_ic.csv')
rot = {'agua_rede': 'Água da rede', 'solucao_sanitaria_estruturada': 'Solução sanitária estruturada', 'rua_pavimentada': 'Rua pavimentada', 'internet_banda_larga_fixa': 'Banda larga fixa',
       'possui_computador': 'Computador', 'possui_automovel': 'Automóvel', 'trabalhou_30d': 'Trabalhou nos últimos 30 dias', 'contribui_previdencia': 'Contribui à previdência'}
fig, axs = plt.subplots(2, 4, figsize=(13, 6.2), sharey=True)
for ax, (v, t) in zip(axs.ravel(), rot.items()):
    for g in range(3):
        r = p[(p.grupo == g) & (p.variavel == v)].iloc[0]
        ax.errorbar(r.estimativa * 100, g, xerr=[[max(0, (r.estimativa - r.ic_inf) * 100)], [max(0, (r.ic_sup - r.estimativa) * 100)]], fmt='o', color=COR[g], ecolor='gray', capsize=3, ms=7)
    ax.set_title(t, fontsize=9); ax.set_yticks(range(3)); ax.set_yticklabels(NOM, fontsize=8); ax.invert_yaxis(); ax.set_xlim(0, 100); ax.grid(axis='x', alpha=.3); ax.set_xlabel('% dos domicílios', fontsize=8)
plt.tight_layout(); plt.savefig('figuras/perfil_grupos_ic.png', dpi=200); plt.close()
v = pd.read_csv('validacao_externa.csv')
rot2 = {'ebia_insegur_qualquer': 'Insegurança alimentar (qualquer nível)', 'ebia_mod_grave': 'Insegurança alimentar moderada ou grave', 'cadunico': 'Inscrito no CadÚnico (algum morador)', 'beneficio_social': 'Recebe ou já recebeu benefício (algum morador)', 'escritura_registrada': 'Escritura registrada em cartório', 'outorga_agua': 'Outorga para uso da água'}
fig, axs = plt.subplots(2, 3, figsize=(13, 6), sharey=True)
for ax, (k, t) in zip(axs.ravel(), rot2.items()):
    for g in range(3):
        r = v[(v.grupo == g) & (v.variavel == k)].iloc[0]
        ax.errorbar(r.estimativa_pct, g, xerr=[[max(0, r.estimativa_pct - r.ic_inf)], [max(0, r.ic_sup - r.estimativa_pct)]], fmt='o', color=COR[g], ecolor='gray', capsize=3, ms=7)
    ax.set_title(t, fontsize=9); ax.set_yticks(range(3)); ax.set_yticklabels(NOM, fontsize=8); ax.invert_yaxis(); ax.grid(axis='x', alpha=.3); ax.set_xlabel('% dos domicílios', fontsize=8)
plt.tight_layout(); plt.savefig('figuras/validacao_externa.png', dpi=200); plt.close()
