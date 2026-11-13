# -*- coding: utf-8 -*-
"""Decomposição territorial no NÍVEL DO DOMICÍLIO, como especificada na avaliação (seção "Decomposição territorial").
Para cada indicador y: efeito fixo da subpopulação + intercepto aleatório de RA + intercepto aleatório de RA x subpopulação (segmento).
 - binários: GLMM logístico (statsmodels BinomialBayesMixedGLM, aproximação variacional); variância residual na escala latente = pi^2/3;
 - contínuos (ln da renda per capita): MixedLM linear (REML) com componentes de variância.
Reporta a fração da variância total devida a subpopulação (fixa), RA, RA x subpopulação (segmento) e residual, com intervalo aproximado (média +- 1,96 DP posterior do log-desvio-padrão; no ajuste
variacional o DP posterior é subestimado, então leia o intervalo como piso de incerteza).
LIMITES: sem pesos amostrais (o estrato é quase coberto por RA x subpopulação); 24 RAs dão pouca informação para a variância de RA. Rodar em resultados/<variante>/ (só usa a Fase 1 do 03). Saída: decomposicao_domiciliar_26.csv"""
import os, warnings, numpy as np, pandas as pd
from pathlib import Path
import statsmodels.formula.api as smf
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
warnings.filterwarnings('ignore')
cfg = dict(kv.split('=', 1) for kv in Path('config.txt').read_text().split())
os.environ.update({'V2': cfg['V2'], 'K_FINAL': cfg['K'], 'WINSOR_P': cfg['WINSOR_P'], 'COMP': cfg['COMP']})
os.environ.setdefault('TCC_BASE', str(Path(__file__).resolve().parent.parent / 'dados' / 'gold')); os.environ.setdefault('TCC_SHP', str(Path(__file__).resolve().parent.parent / 'dados' / 'regioes_administrativas.zip'))
src = open(Path(__file__).resolve().parent / '03_analise_cluster_4d.py', encoding='utf-8').read()
exec(compile(src[:src.index('def monte_carlo_parallel_analysis')], '03_fase1', 'exec'))
d = df_dom.copy(); d['ra'] = d[RA].astype(str); d['sub'] = d[SUB].astype(str); d['seg'] = d['ra'] + '|' + d['sub']
d['ln_renda_pc'] = np.where(d['renda_domiciliar_pc'] >= 0, np.log1p(d['renda_domiciliar_pc']), np.nan)
BIN = [v for v in ['agua_rede', 'internet_banda_larga_fixa', 'solucao_sanitaria_estruturada', 'possui_computador'] if v in d.columns]
linhas = []
def frac(nome, tipo, fe, vra, vseg, vres, ic):
    tot = fe + vra + vseg + vres
    linhas.append({'indicador': nome, 'modelo': tipo, 'n_dom': int(ic['n']), 'subpop(fixa)': fe / tot, 'RA': vra / tot, 'RA x subpop (segmento)': vseg / tot, 'residual': vres / tot,
                   'RA_ic_baixo': ic.get('ra_lo', np.nan), 'RA_ic_alto': ic.get('ra_hi', np.nan), 'segmento_ic_baixo': ic.get('seg_lo', np.nan), 'segmento_ic_alto': ic.get('seg_hi', np.nan)})
for y in BIN:
    x = d[[y, 'ra', 'sub', 'seg']].dropna().copy(); x[y] = x[y].astype(int)
    m = BinomialBayesMixedGLM.from_formula(f'{y} ~ C(sub)', {'ra': '0 + C(ra)', 'seg': '0 + C(seg)'}, x).fit_vb()
    sd_ra, sd_seg = np.exp(m.vcp_mean[0]), np.exp(m.vcp_mean[1]); s_ra, s_seg = m.vcp_sd[0], m.vcp_sd[1]
    X = m.model.exog; fe = np.var(X @ m.fe_mean)
    ic = {'n': len(x), 'ra_lo': np.exp(2 * (m.vcp_mean[0] - 1.96 * s_ra)), 'ra_hi': np.exp(2 * (m.vcp_mean[0] + 1.96 * s_ra)), 'seg_lo': np.exp(2 * (m.vcp_mean[1] - 1.96 * s_seg)), 'seg_hi': np.exp(2 * (m.vcp_mean[1] + 1.96 * s_seg))}
    tot_ref = fe + sd_ra ** 2 + sd_seg ** 2 + np.pi ** 2 / 3
    for k in ('ra_lo', 'ra_hi', 'seg_lo', 'seg_hi'): ic[k] /= tot_ref          # converte para fração (aprox.)
    frac(y, 'logit (VB)', fe, sd_ra ** 2, sd_seg ** 2, np.pi ** 2 / 3, ic)
for y in ['ln_renda_pc']:
    x = d[[y, 'ra', 'sub', 'seg']].dropna().copy()
    for met in ('powell', 'bfgs', 'nm'):
        try: m = smf.mixedlm(f'{y} ~ C(sub)', x, groups=x['ra'], re_formula='1', vc_formula={'seg': '0 + C(seg)'}).fit(reml=True, method=met); break
        except Exception: m = None
    if m is None: continue
    fe = np.var(m.model.exog @ m.fe_params.values); frac(y, 'linear (REML)', fe, float(m.cov_re.iloc[0, 0]), float(m.vcomp[0]), float(m.scale), {'n': len(x)})
R = pd.DataFrame(linhas); R.to_csv('decomposicao_domiciliar_26.csv', index=False)
print('DECOMPOSIÇÃO DOMICILIAR (frações da variância total; binários na escala logit latente):'); print(R.round(3).to_string(index=False))
pd.DataFrame([{'dec26_fracao_RA_media': R['RA'].mean(), 'dec26_fracao_segmento_media': R['RA x subpop (segmento)'].mean(), 'dec26_fracao_subpop_media': R['subpop(fixa)'].mean(), 'dec26_fracao_residual_media': R['residual'].mean()}]).to_csv('decomposicao_domiciliar_26_resumo.csv', index=False)
