# Denominadores das variáveis de moradores usadas em D2

As médias por segmento (`media_ponderada_vetorizada` em `pipeline/03`) usam o peso do morador e ignoram valores ausentes. A função `parse_binaria` aceita apenas "Sim"/"Não" (e 1/0); **"Não se aplica" e "Não sabe" viram valor ausente e saem do denominador**.

Contagem nos microdados (moradores, n = 15.357):

| Variável | Sim | Não | Não se aplica | Não sabe | Denominador efetivo |
|---|---|---|---|---|---|
| `tem_ctps` (I17) | 2.051 | 582 | 12.680 | 44 | **só quem respondeu Sim/Não (2.633 moradores; provavelmente os ocupados, o que o questionário não confirma)**; não é a proporção de todos os moradores |
| `trabalhou_30d` (I05) | 6.044 | 7.222 | 2.034 | 57 | moradores a que a pergunta se aplica (13.266) |
| `contribui_previdencia` (I18) | 3.663 | 8.585 | 2.986 | 123 | moradores a que a pergunta se aplica (12.248) |
| `escolaridade_adulto` | – | – | – | – | moradores com 25 anos ou mais |

Consequência: `tem_ctps` mede a formalização **entre quem tem a pergunta aplicável**, e segmentos com poucos ocupados têm essa média estimada com poucas observações.
