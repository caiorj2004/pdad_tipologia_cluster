# Reprodução do saneamento (programa 01) a partir do dicionário aplicado

Entrada: `dados/silver/Dicionario_aplicado_PDADR_2022.xlsx` (abas `PDADR_2022_Domicilios`, 5.075 × 140, e `PDADR 2022 Moradores`, 15.357 × 150), fornecida pelo IPEDF.
Comando: `python reproduzir_tudo.py --saneamento` (ou `PDADR_SILVER=... python pipeline/01_saneamento_pdadr_2022.py`). Tempo: cerca de 1 minuto.

Resultado da comparação com `dados/gold/` (as planilhas usadas em todos os resultados deste repositório), feita em 13/11/2026:

| Base | Formato | Colunas comuns com valores diferentes | Observação |
|---|---|---|---|
| Domicílios | 5.075 × 208 (gold: 206) | **nenhuma** | o programa 01 atual grava 2 colunas a mais de transparência da renda: `n_moradores_renda_ausente` e `renda_dom_parcial` |
| Moradores | 15.357 × 186 (gold: 186) | `E23_1__Local_intencao_novo_domicilio` (8 linhas) e `I08__Local_trabalho_principal` (22): "Recanto das Emas" no lugar de "Recanto Das Emas"; `I04_1_1…_classe` (1 linha): rótulo "Isento / Ausente" no lugar de "Não recebeu" | diferenças de grafia e de rótulo em colunas que **não entram** no agrupamento nem na renda |

Conclusão: o programa 01 reproduz, a partir da planilha do IPEDF, todas as colunas das bases `gold` que alimentam o modelo (inclusive a renda domiciliar per capita). As poucas diferenças são de grafia/rótulo em colunas fora do modelo, mais duas colunas novas. Por isso `dados/gold/` foi mantido como entrada dos resultados já gerados; refazer todos os resultados com a saída atual do programa 01 não deve alterá-los, mas isso não foi verificado rodando o lote de novo.
