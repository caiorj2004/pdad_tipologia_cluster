#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Executa o pipeline do TCC de ponta a ponta com UM comando.

    python reproduzir_tudo.py --modo teste        # verificação rápida da instalação (~5 min)
    python reproduzir_tudo.py --modo principal    # os dois modelos centrais: M0_antigo_K5 (avaliado) e M1_revisado_K3 (revisado) + figuras do TCC
    python reproduzir_tudo.py --modo completo     # todas as variantes de pipeline/variantes.tsv (algumas horas) + comparação

Opções úteis:
    --paralelo N          variantes simultâneas (padrão 2)
    --variantes a,b,c     roda só estas variantes (nomes em pipeline/variantes.tsv)
    --pular-existentes    não refaz variantes que já têm resultados/<variante>/stress_resumo.csv
    --saneamento          antes de tudo, roda o programa 01 (exige dados/silver/Dicionario_aplicado_PDADR_2022.xlsx) e compara com dados/gold/
    --figuras-tcc         roda também os programas 02, 04–13 e 20 (tabelas e figuras do TCC) sobre M1_revisado_K3
    --saida DIR           pasta de saída (padrão: resultados/)

Etapas por variante (rodam em resultados/<variante>/, na ordem):
    03 modelo | 16 bootstrap, nulos, encolhimento | 18 variância, AMI, join count | 19 validação externa | 21 testes de estresse | 23 deixa-uma-RA-fora refazendo o PCA
Ao final: 22 consolida todas as variantes em comparacao_variantes.{csv,md}.
"""
import argparse, os, subprocess, sys, time, shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

RAIZ = Path(__file__).resolve().parent; PIPE = RAIZ / 'pipeline'
ETAPAS = [('03', '03_analise_cluster_4d.py'), ('16', '16_robustez_amostral_e_pesos.py'), ('18', '18_variancia_ami_espacial.py'),
          ('19', '19_validacao_externa.py'), ('21', '21_stress_estrutural.py'), ('23', '23_deixa_uma_RA_refazendo_pca.py')]
FIGURAS_TCC = ['04_comparacao_metodos.py', '05_ward_dendrograma.py', '06_ward_bateria.py', '07_estabilidade_cruzada.py', '08_variaveis_apoio.py',
               '09_kruskal_wallis_escores.py', '10_tabelas_e_figuras.py', '11_mapas_indicadores.py', '12_robustez_tandem.py', '13_ra_grupos_componentes.py',
               '20_figuras_perfil_validacao.py']
PRINCIPAIS = ['M0_antigo_K5', 'M1_revisado_K3']

def ler_variantes():
    v = {}
    for ln in (PIPE / 'variantes.tsv').read_text(encoding='utf-8').splitlines():
        if not ln.strip() or ln.startswith('#'): continue
        c = ln.split('\t'); nome, V2, K, win, comp, nb, ni, fig = c[:8]
        v[nome] = dict(nome=nome, V2='' if V2 == 'nenhuma' else V2, K=K, WINSOR_P=win, COMP=comp, N_BOOT=nb, NINIT_BOOT=ni, figuras=fig == '1')
    return v

def ambiente(cfg):
    e = os.environ.copy()
    e.update(TCC_BASE=str(RAIZ / 'dados' / 'gold'), TCC_SHP=str(RAIZ / 'dados' / 'regioes_administrativas.zip'), TCC_03=str(PIPE / '03_analise_cluster_4d.py'),
             V2=cfg['V2'], K_FINAL=cfg['K'], WINSOR_P=cfg['WINSOR_P'], COMP=cfg['COMP'], N_BOOT=cfg['N_BOOT'], NINIT_BOOT=cfg['NINIT_BOOT'],
             MPLBACKEND='Agg', PYTHONHASHSEED='0', PYTHONIOENCODING='utf-8')
    return e

def rodar(script, cwd, env, log):
    t = time.time()
    with open(log, 'w', encoding='utf-8') as f:
        r = subprocess.run([sys.executable, str(PIPE / script)], cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT)
    return r.returncode, time.time() - t

def variante(cfg, saida, figuras_tcc):
    out = saida / cfg['nome']; shutil.rmtree(out, ignore_errors=True); out.mkdir(parents=True)
    (out / 'config.txt').write_text(f"variante={cfg['nome']} V2={cfg['V2']} K={cfg['K']} WINSOR_P={cfg['WINSOR_P']} COMP={cfg['COMP']} N_BOOT={cfg['N_BOOT']} NINIT_BOOT={cfg['NINIT_BOOT']}\n", encoding='utf-8')
    env = ambiente(cfg); tempos = []
    for cod, script in ETAPAS:
        rc, dt = rodar(script, out, env, out / f'saida_{cod}.txt'); tempos.append((cod, rc, dt))
        if rc != 0: return cfg['nome'], False, tempos
    if figuras_tcc and cfg['nome'] == 'M1_revisado_K3':
        for script in FIGURAS_TCC:
            rc, dt = rodar(script, out, env, out / f"saida_{script[:2]}.txt"); tempos.append((script[:2], rc, dt))
            if rc != 0: return cfg['nome'], False, tempos
    shutil.rmtree(out / 'shape_df_extracao', ignore_errors=True)
    if not (cfg['figuras'] or figuras_tcc and cfg['nome'] == 'M1_revisado_K3'):
        for p in out.rglob('*.png'): p.unlink()
    return cfg['nome'], True, tempos

def saneamento():
    silver = RAIZ / 'dados' / 'silver' / 'Dicionario_aplicado_PDADR_2022.xlsx'
    if not silver.exists():
        print(f'[saneamento] ARQUIVO AUSENTE: {silver}\n  Coloque a planilha-dicionário do IPEDF nessa pasta (veja README) ou rode sem --saneamento.'); return False
    out = RAIZ / 'dados' / 'gold_reproduzido'; out.mkdir(exist_ok=True)
    env = os.environ.copy(); env.update(PDADR_SILVER=str(silver), PDADR_OUT_DIR=str(out), PDADR_EXPORTAR='1', MPLBACKEND='Agg', PYTHONIOENCODING='utf-8')
    rc, dt = rodar('01_saneamento_pdadr_2022.py', out, env, out / 'saida_01.txt'); print(f'[saneamento] programa 01: rc={rc} ({dt:.0f}s)')
    if rc != 0: return False
    import pandas as pd, numpy as np
    ok = True
    for nome in ('pdadr_domicilios.xlsx', 'pdadr_moradores.xlsx'):
        a = pd.read_excel(out / nome); b = pd.read_excel(RAIZ / 'dados' / 'gold' / nome); dif = []
        for c in b.columns:
            if c not in a.columns: dif.append(c + ' (ausente)'); continue
            x, y = a[c], b[c]
            if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y): eq = np.isclose(x.astype(float), y.astype(float), equal_nan=True, rtol=1e-9, atol=1e-9).all()
            else: eq = ((x.isna() & y.isna()) | (x.astype(str) == y.astype(str))).all()
            if not eq: dif.append(c)
        extras = [c for c in a.columns if c not in b.columns]
        print(f'[saneamento] {nome}: colunas diferentes nas comuns: {dif or "nenhuma"} | colunas só na reprodução: {extras or "nenhuma"}'); ok &= not dif
    print('[saneamento] Diferenças esperadas (docs/reproducao_saneamento.md): grafia de RA em E23_1/I08, 1 rótulo em I04_1_1 e 2 colunas de transparência de renda.')
    return ok

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--modo', choices=['teste', 'principal', 'completo'], default='principal')
    ap.add_argument('--paralelo', type=int, default=2); ap.add_argument('--variantes', default='')
    ap.add_argument('--pular-existentes', action='store_true'); ap.add_argument('--saneamento', action='store_true')
    ap.add_argument('--figuras-tcc', action='store_true'); ap.add_argument('--saida', default='')
    a = ap.parse_args(); todas = ler_variantes()
    if a.saneamento: saneamento()
    if a.modo == 'teste':      # smoke test: poucas réplicas, saída separada
        sel = {'M1_revisado_K3': dict(todas['M1_revisado_K3'], N_BOOT='20', NINIT_BOOT='5')}; saida = Path(a.saida or RAIZ / 'resultados_teste')
    else:
        nomes = [n for n in a.variantes.split(',') if n] or (PRINCIPAIS if a.modo == 'principal' else list(todas))
        falta = [n for n in nomes if n not in todas]
        if falta: sys.exit(f'variantes desconhecidas: {falta}; veja pipeline/variantes.tsv')
        sel = {n: todas[n] for n in nomes}; saida = Path(a.saida or RAIZ / 'resultados')
    saida.mkdir(parents=True, exist_ok=True)
    if a.pular_existentes: sel = {n: c for n, c in sel.items() if not (saida / n / 'stress_resumo.csv').exists()}
    fig = a.figuras_tcc or a.modo == 'principal'
    print(f'{len(sel)} variante(s): {", ".join(sel)} | saída: {saida}', flush=True); t0 = time.time(); falhas = []
    with ThreadPoolExecutor(max_workers=max(1, a.paralelo)) as ex:
        for nome, ok, tempos in ex.map(lambda c: variante(c, saida, fig), sel.values()):
            print(f'  {"OK    " if ok else "FALHOU"} {nome}  ' + ' '.join(f'{c}:{dt:.0f}s' for c, rc, dt in tempos), flush=True)
            if not ok: falhas.append(nome)
    if a.modo != 'teste':
        env = os.environ.copy(); env['PYTHONIOENCODING'] = 'utf-8'
        subprocess.run([sys.executable, str(PIPE / '22_comparar_variantes.py')], env=env, cwd=RAIZ)
    print(f'Concluído em {(time.time() - t0) / 60:.1f} min. Falhas: {falhas or "nenhuma"}')
    sys.exit(1 if falhas else 0)

if __name__ == '__main__':
    main()
