"""Gera docs/tabelas_rodadas.md a partir de resultados/comparacao_variantes.csv e das partições de cada variante.
Nenhum número é digitado à mão: o texto de docs/rodadas_stress_test.md cita estas tabelas.
Uso: python pipeline/27_tabelas_rodadas.py   (depois de 22_comparar_variantes.py)"""
import os, pandas as pd, numpy as np
from sklearn.metrics import adjusted_rand_score as ari
BASE = os.environ.get('TCC_BASE', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
R = os.path.join(BASE, 'resultados'); D = os.path.join(BASE, 'docs')
d = pd.read_csv(os.path.join(R, 'comparacao_variantes.csv'), index_col=0)
ordem = ['M0_antigo_K3','M0_antigo_K4','M0_antigo_K5','M1_revisado_K3','M1_revisado_K4','M1_revisado_K5',
         'A_so_ctx','A_so_agua','A_so_esgoto','A_so_renda','A_so_dur','L_sem_ctx','L_sem_agua','L_sem_esgoto','L_sem_renda','L_sem_dur',
         'S_winsor95','S_winsor100','S_clr','E_agua_cisterna_precaria','E_agua_nascente_individual','E_agua_gambiarra_na_rede','E_esgoto_rudimentar_individual',
         'M0_antigo_K5_b500n50','M0_antigo_K5_bootEstratoRegra','M1_revisado_K3_b500n50','M1_revisado_K3_bootEstratoRegra']
ordem = [o for o in ordem if o in d.index]; d = d.loc[ordem]
def f(x, n=3):
    if isinstance(x, str): return x
    if pd.isna(x): return '–'
    if float(x).is_integer() and abs(x) >= 1: return str(int(x))
    return f'{x:.{n}f}'
def tab(titulo, cols, nomes=None, nota=None, n=3):
    nomes = nomes or cols; cols_ok = [c for c in cols if c in d.columns]
    nomes = [nm for c, nm in zip(cols, nomes) if c in d.columns]
    s = f'\n### {titulo}\n\n| variante | ' + ' | '.join(nomes) + ' |\n|---|' + '---|' * len(nomes) + '\n'
    for v in d.index: s += f'| {v} | ' + ' | '.join(f(d.loc[v, c], n) for c in cols_ok) + ' |\n'
    return s + (f'\n{nota}\n' if nota else '')
o = ['# Tabelas das rodadas de estresse (geradas automaticamente)\n',
     'Geradas por `pipeline/27_tabelas_rodadas.py` a partir de `resultados/comparacao_variantes.csv`. Sem interpretação; a leitura está em `docs/rodadas_stress_test.md`. "–" = métrica não se aplica ou não foi calculada.\n']
o.append(tab('T1. Construção e partição', ['V2','K','WINSOR_P','COMP','n_componentes','tamanhos_grupos(seg)','pct_dom_menor_grupo'],
    ['flags','K','winsor','comp.','PCs','segmentos por grupo','% dom. exp. do menor grupo'], n=2))
o.append(tab('T2. Validade interna no K final', ['silhueta','p_embaralhado','p_covariancia','CH','DB','ARI_subamostra80','p_cov_minimo_entre_K','K_do_p_cov_minimo','n_K_com_p_cov<0.05'],
    ['silhueta','p (nulo embaralhado)','p (nulo covariância)','CH','DB','ARI subamostra 80%','menor p_cov entre K=2..8','K desse p','nº de K com p_cov<0,05'],
    'Nulo de covariância: dados sintéticos com a mesma covariância e sem grupos. p alto = a silhueta observada é compatível com ausência de estrutura de grupos.'))
o.append(tab('T3. Bootstrap de domicílios (fluxo inteiro refeito) e nulos', ['boot_ARI_medio','boot_ARI_p5','boot_ARI_obs_medio','boot_ARI_nulo_p95','boot_ARI_p_vs_nulo','metades_ARI_obs_medio','metades_ARI_nulo_p95','encolhimento_ARI','encolhimento_mudam','n_fronteira_boot','pct_dom_fronteira'],
    ['ARI médio (16)','ARI p5','ARI obs. (24)','nulo p95','p vs nulo','metades obs.','metades nulo p95','ARI após encolhimento','segs. que mudam','segs. na fronteira','% dom. na fronteira'],
    '"ARI médio (16)" e "ARI obs. (24)" vêm de rotinas distintas (16: realinhamento húngaro com NINIT_BOOT inicializações; 24: cadeia PCA fixo + K-Means). Comparar cada um com seu próprio nulo.'))
o.append(tab('T4. Algoritmo e semente', ['S1_pct_sementes_no_otimo','S1_ARI_medio_semente_unica','S1_ARI_min_semente_unica','S1_pct_ninit10_no_otimo','S7_ARI_Ward','S7_ARI_GMM_completo','S7_ARI_KMeans_peso_integral','S7_ARI_medio_outros_algoritmos','S7_ARI_min_outros_algoritmos','S4_aparado6_ARI','S4_aparado10_ARI','S5_ARI_argmax_vs_kmeans','S5_pureza_media_Kfinal'],
    ['% sementes no ótimo','ARI semente única (média)','ARI semente única (mín.)','% com 10 inic. no ótimo','ARI Ward','ARI GMM','ARI K-Means peso integral','ARI outros alg. (média)','ARI outros alg. (mín.)','ARI aparado 6%','ARI aparado 10%','ARI arquétipos vs K-Means','pureza arquétipos']))
o.append(tab('T5. Dependência de segmentos, RAs e tamanho', ['S6_ARI_medio_sem_uma_RA','S6b_ARI_medio_sem_uma_RA','S6b_ARI_min_sem_uma_RA','S6b_RA_critica','S6b_pct_RAs_ARI<0.5','S6b_n_componentes_varia','S6b_ARI_medio_sem_uma_subpop','S6b_ARI_min_sem_uma_subpop','S8_ARI_n>=10','S8_ARI_n>=20','S2_spearman_n_dist_Kfinal','S2_p_Kfinal','S3_largura_media_IC_n<12'],
    ['sem uma RA (PCA fixo)','sem uma RA (refaz PCA)','mín.','RA crítica','% RAs com ARI<0,5','nº de PCs varia?','sem uma subpop (refaz PCA)','mín.','ARI corte n≥10','ARI corte n≥20','Spearman n × dist.','p','largura média IC Wilson (n<12)']))
o.append(tab('T6. Confusão territorial', ['AMI_RA_obs','AMI_RA_p_atrib','AMI_seg_p_RA','AMI_subpop_obs','AMI_subpop_p_atrib','AMI_seg_p_subpop','S9_frac_mesmo_grupo_RA','S9_frac_esperada_RA','S9_p_RA','S9_frac_mesmo_grupo_subpop','S9_frac_esperada_subpop','S9_p_subpop','joincount_obs','joincount_esperado','joincount_p'],
    ['AMI RA','p (atrib. domicílios)','p (perm. segmentos)','AMI subpop','p (atrib.)','p (perm.)','pares mesma RA no mesmo grupo','esperado','p','pares mesma subpop no mesmo grupo','esperado','p','join count obs.','esperado','p']))
o.append(tab('T7. Validação externa (Kruskal-Wallis por segmento, p) e ordenada', [c for c in d.columns if c.startswith('VE_KW')] + ['IVS_KW_p','ordenada_n_p<0.05','ordenada_rho_medio_abs'],
    [c.replace('VE_KW_p_','') for c in d.columns if c.startswith('VE_KW')] + ['IVS-DF','nº desfechos ordenados p<0,05','|rho| médio']))
o.append(tab('T9. Validação externa sob incerteza da partição (programa 28)', ['ext28_n_desfechos_p_ref_menor_005','ext28_n_desfechos_mediana_p_menor_005','ext28_pct_replicas_p005_media','ext28_epsilon2_mediano_medio'],
    ['desfechos com p<0,05 (partição de referência)','desfechos com p mediano<0,05 (réplicas)','% médio de réplicas com p<0,05','epsilon² mediano médio'],
    'De 7 desfechos (Kruskal-Wallis entre grupos, 49 segmentos). Desfechos fixos; só a partição varia, com os rótulos de cada réplica do bootstrap. Sem correção para comparações múltiplas.'))
# partições: ARI entre variantes
L = {}
for v in d.index:
    s = pd.read_csv(os.path.join(R, v, 'segmentos_completo.csv')); L[v] = s[f'Cluster_K{int(d.loc[v, "K"])}'].values
M = pd.DataFrame([[ari(L[a], L[b]) for b in d.index] for a in d.index], index=d.index, columns=d.index)
M.round(4).to_csv(os.path.join(R, 'ari_entre_particoes.csv'))
s = '\n### T8. ARI entre a partição de cada variante (no seu K) e as duas de referência\n\n| variante | ARI vs M0_antigo_K5 | ARI vs M1_revisado_K3 |\n|---|---|---|\n'
for v in d.index: s += f'| {v} | {M.loc[v,"M0_antigo_K5"]:.3f} | {M.loc[v,"M1_revisado_K3"]:.3f} |\n'
o.append(s + '\nMatriz completa em `resultados/ari_entre_particoes.csv`. ARI não depende da rotulação dos grupos; compara partições de 49 segmentos.\n')
open(os.path.join(D, 'tabelas_rodadas.md'), 'w', encoding='utf-8').write('\n'.join(o))
print('docs/tabelas_rodadas.md escrito;', len(d), 'variantes')
