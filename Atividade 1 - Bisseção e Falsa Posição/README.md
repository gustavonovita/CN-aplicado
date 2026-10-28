# Métodos da Bisseção e da Falsa Posição: Escoamento em Canal Aberto

Atividade prática da disciplina **Cálculo Numérico Aplicado**. O objetivo é determinar, com dois métodos numéricos intervalares (**bisseção** e **falsa posição**), a profundidade da lâmina d'água em um canal retangular (equação de Manning) e, a partir dela, a velocidade média do escoamento.

## Sumário

- [Problema físico](#problema-físico)
- [Formulação matemática](#formulação-matemática)
- [Métodos numéricos](#métodos-numéricos)
- [Estrutura do repositório](#estrutura-do-repositório)
- [Como executar](#como-executar)
- [Explicação dos códigos](#explicação-dos-códigos)
- [Resultados](#resultados)
- [Comparação dos métodos](#comparação-dos-métodos)
- [Conclusão](#conclusão)

## Problema físico

Considere o escoamento uniforme de água em um canal aberto retangular de largura $B$, inclinado em relação à horizontal. Para uma vazão volumétrica conhecida $Q$, deseja-se determinar a profundidade $H$ da lâmina d'água e a velocidade média $U$ do escoamento.

**Parâmetros:**

| Símbolo | Descrição | Valor |
|---|---|---|
| $S$ | Inclinação do canal | 0,0002 |
| $B$ | Largura do canal | 20 m |
| $n$ | Coeficiente de rugosidade de Manning | 0,03 |
| $Q$ | Vazão volumétrica | 5 m³/s |

## Formulação matemática

A vazão que atravessa a seção transversal é

$$Q = U A_c = U B H$$

em que $A_c = BH$ é a área molhada. A velocidade média é estimada pela equação semiempírica de Manning:

$$U = \frac{1}{n} R^{2/3} S^{1/2}, \qquad R = \frac{A_c}{P} = \frac{BH}{B + 2H}$$

sendo $R$ o raio hidráulico e $P = B + 2H$ o perímetro molhado. Combinando as expressões, obtém-se uma relação não linear entre $Q$ e $H$:

$$Q = \frac{S^{1/2}}{n}\,\frac{(BH)^{5/3}}{(B+2H)^{2/3}}$$

### Problema de determinação de zero de função

Substituindo os parâmetros, o problema vira $f(H) = 0$, com

$$f(H) = 0{,}471405\,\frac{(20H)^{5/3}}{(20+2H)^{2/3}} - 5$$

O coeficiente vem de $\sqrt{S}/n = \sqrt{0{,}0002}/0{,}03 \approx 0{,}471405$. Fisicamente, a raiz de $f$ é a profundidade para a qual a vazão de Manning é igual à vazão prescrita de 5 m³/s.

Como $f$ é crescente em $H$ (mais profundidade, mais vazão), existe no máximo uma raiz.

### Velocidade média

Depois de obtido $H$, usa-se a relação entre vazão e área:

$$U = \frac{Q}{BH}$$

## Métodos numéricos

Ambos os métodos partem do mesmo intervalo inicial $[H_l, H_u] = [0{,}5;\ 1{,}0]$ m e usam o mesmo critério de parada, $|f(H_r)| < 10^{-6}$, com limite de 100 iterações.

**Verificação da mudança de sinal:**

$$f(0{,}5) \approx -2{,}1254, \qquad f(1{,}0) \approx +3{,}8477 \quad\Rightarrow\quad f(H_l)\,f(H_u) < 0$$

Pelo Teorema do Valor Intermediário, como $f$ é contínua, há ao menos uma raiz no intervalo.

### Método da Bisseção

A nova aproximação é o ponto médio do intervalo:

$$H_r = \frac{H_l + H_u}{2}$$

### Método da Falsa Posição

A nova aproximação é a interseção com o eixo horizontal da reta que passa por $(H_l, f(H_l))$ e $(H_u, f(H_u))$:

$$H_r = H_u - f(H_u)\,\frac{H_u - H_l}{f(H_u) - f(H_l)}$$

### Atualização do intervalo (comum aos dois)

- Se $f(H_l)\,f(H_r) < 0$, a raiz está em $[H_l, H_r]$, então $H_u = H_r$.
- Caso contrário, a raiz está em $[H_r, H_u]$, então $H_l = H_r$.

Assim, o intervalo sempre preserva a mudança de sinal e a raiz nunca sai dele.

## Estrutura do repositório

```
.
├── README.md
├── bissecao.py          # Método da Bisseção
└── falsa_posicao.py     # Método da Falsa Posição
```

## Como executar

Requisitos: Python 3 (apenas biblioteca padrão, sem dependências externas).

```bash
python bissecao.py
python falsa_posicao.py
```

## Explicação dos códigos

### `bissecao.py`

```python
H_u = 1
H_l = 0.5
max_iter = 100
tol = 1e-6

def f(H):
    numerador = (20*H)**(5/3)
    denominador = (20 + 2*H)**(2/3)
    return 0.471405 * (numerador / denominador) - 5

for k in range(1, max_iter + 1):
    H_r = (H_l + H_u) / 2

    print(f"Iteração {k}: Hr = {H_r:.6f}, f(H_r) = {f(H_r):.6f}")

    if abs(f(H_r)) < tol:
        print(f"Raiz encontrada: H = {H_r:.6f} após {k} iterações.")
        break

    if f(H_l) * f(H_r) < 0:
        H_u = H_r
    else:
        H_l = H_r
```

**Como funciona:**

1. Define o intervalo inicial (`H_l`, `H_u`), a tolerância (`tol`) e o número máximo de iterações (`max_iter`), que evita um laço infinito em caso de erro.
2. A função `f(H)` implementa a equação de Manning já substituída com os parâmetros do problema, menos a vazão de 5 m³/s.
3. A cada iteração, calcula o ponto médio `H_r` e imprime `H_r` e `f(H_r)`.
4. Se `|f(H_r)| < tol`, o programa imprime a raiz e o número de iterações e para com `break`.
5. Caso contrário, o teste `f(H_l) * f(H_r) < 0` decide qual metade do intervalo contém a raiz, e o limite correspondente é substituído por `H_r`.

### `falsa_posicao.py`

```python
H_u = 1
H_l = 0.5
tol = 1e-6
max_iter = 100

def f(H):
    numerador = (20*H)**(5/3)
    denominador = (20 + 2*H)**(2/3)
    return 0.471405 * (numerador / denominador) - 5

for k in range(1, max_iter + 1):
    H_r = H_u - (f(H_u) * ((H_u - H_l) / (f(H_u) - f(H_l))))

    if abs(f(H_r)) < tol:
        print(f"Raiz encontrada: H = {H_r:.6f} após {k} iterações.")
        break

    print(f"Iteração {k}: Limite inferior = {H_l:.6f}, Limite superior = {H_u:.6f}, Hr = {H_r:.6f}, f(Hr) = {f(H_r):.6f}")

    if f(H_l) * f(H_r) < 0:
        H_u = H_r
    else:
        H_l = H_r
```

**Como funciona:** a estrutura é idêntica à da bisseção. A única diferença é a linha que calcula `H_r`, que usa a fórmula da falsa posição em vez do ponto médio.

**Observação sobre a saída:** neste código, o teste de convergência vem antes do `print`. Por isso, a iteração em que o critério é satisfeito não aparece na tabela impressa, só na mensagem final (por exemplo, 6 linhas de iteração e "após 7 iterações").

## Resultados

### Resumo

| | Bisseção | Falsa posição |
|---|---|---|
| Profundidade $H$ (m) | 0,702293 | 0,702293 |
| Iterações | 22 | 7 |
| Velocidade média $U$ (m/s) | 0,3560 | 0,3560 |

**Resposta final:** $H \approx 0{,}7023$ m e $U \approx 0{,}356$ m/s.

A velocidade foi calculada com $U = Q/(BH) = 5/(20 \times 0{,}702293) \approx 0{,}35598$ m/s.

**Verificação por Manning:** com $H = 0{,}7023$ m, $A_c = 14{,}05$ m², $P = 21{,}405$ m e $R = 0{,}6562$ m, logo $U = \frac{1}{0{,}03}(0{,}6562)^{2/3}(0{,}0002)^{1/2} \approx 0{,}356$ m/s, o que confere.

### Saída de `bissecao.py`

<details>
<summary>Ver saída completa</summary>

```
Iteração 1: Hr = 0.750000, f(H_r) = 0.562294
Iteração 2: Hr = 0.625000, f(H_r) = -0.863133
Iteração 3: Hr = 0.687500, f(H_r) = -0.169843
Iteração 4: Hr = 0.718750, f(H_r) = 0.191480
Iteração 5: Hr = 0.703125, f(H_r) = 0.009619
Iteração 6: Hr = 0.695312, f(H_r) = -0.080414
Iteração 7: Hr = 0.699219, f(H_r) = -0.035473
Iteração 8: Hr = 0.701172, f(H_r) = -0.012946
Iteração 9: Hr = 0.702148, f(H_r) = -0.001668
Iteração 10: Hr = 0.702637, f(H_r) = 0.003974
Iteração 11: Hr = 0.702393, f(H_r) = 0.001153
Iteração 12: Hr = 0.702271, f(H_r) = -0.000258
Iteração 13: Hr = 0.702332, f(H_r) = 0.000447
Iteração 14: Hr = 0.702301, f(H_r) = 0.000095
Iteração 15: Hr = 0.702286, f(H_r) = -0.000081
Iteração 16: Hr = 0.702293, f(H_r) = 0.000007
Iteração 17: Hr = 0.702290, f(H_r) = -0.000037
Iteração 18: Hr = 0.702291, f(H_r) = -0.000015
Iteração 19: Hr = 0.702292, f(H_r) = -0.000004
Iteração 20: Hr = 0.702293, f(H_r) = 0.000001
Iteração 21: Hr = 0.702293, f(H_r) = -0.000002
Iteração 22: Hr = 0.702293, f(H_r) = -0.000000
Raiz encontrada: H = 0.702293 após 22 iterações.
```

</details>

### Saída de `falsa_posicao.py`

```
Iteração 1: Limite inferior = 0.500000, Limite superior = 1.000000, Hr = 0.677914, f(Hr) = -0.278747
Iteração 2: Limite inferior = 0.677914, Limite superior = 1.000000, Hr = 0.699671, f(Hr) = -0.030257
Iteração 3: Limite inferior = 0.699671, Limite superior = 1.000000, Hr = 0.702015, f(Hr) = -0.003215
Iteração 4: Limite inferior = 0.702015, Limite superior = 1.000000, Hr = 0.702263, f(Hr) = -0.000341
Iteração 5: Limite inferior = 0.702263, Limite superior = 1.000000, Hr = 0.702290, f(Hr) = -0.000036
Iteração 6: Limite inferior = 0.702290, Limite superior = 1.000000, Hr = 0.702292, f(Hr) = -0.000004
Raiz encontrada: H = 0.702293 após 7 iterações.
```

### Como interpretar

- **Bisseção:** o intervalo cai pela metade a cada passo ($0{,}5 \to 0{,}25 \to 0{,}125 \to \dots$). Na iteração 16, $|f|$ chega a $7 \times 10^{-6}$ e na 17 volta a subir ($3{,}7 \times 10^{-5}$). A bisseção não garante que $|f(H_r)|$ diminua a cada passo, só que o intervalo diminui.
- **Falsa posição:** o limite superior fica fixo em $H_u = 1{,}0$ em todas as iterações, e só o limite inferior avança em direção à raiz. Isso é típico do método quando a função é côncava ou convexa no intervalo. Aqui não prejudica, pois $|f(H_r)|$ cai por um fator de cerca de 10 a cada iteração.

## Comparação dos métodos

**Os dois métodos convergiram para a mesma profundidade?**
Sim, $H \approx 0{,}702293$ m. Pequenas diferenças na 8ª casa decimal vêm apenas do critério de parada $|f| < 10^{-6}$.

**Qual precisou de menos iterações?**
A falsa posição: 7 iterações contra 22 da bisseção. Como $f$ é suave e quase linear no intervalo, a interpolação linear aproxima bem a raiz.

**Qual a principal diferença na escolha de $H_r$?**
A bisseção usa apenas a geometria do intervalo (ponto médio) e ignora os valores de $f$. A falsa posição usa os valores de $f$ nos extremos, aproximando a função por uma reta. Se $|f(H_l)|$ é menor que $|f(H_u)|$, o novo ponto fica mais próximo de $H_l$.

**Por que ambos são métodos intervalares?**
Porque trabalham com um intervalo $[H_l, H_u]$ que sempre contém a raiz e o atualizam a cada iteração. Precisam de dois pontos iniciais que "cerquem" a raiz, ao contrário de métodos abertos, como Newton, que partem de um único ponto.

**Qual o papel da condição $f(H_l)\,f(H_u) < 0$?**
Antes de iniciar, garante a existência de raiz no intervalo. Durante as iterações, o teste $f(H_l)\,f(H_r) < 0$ decide qual subintervalo mantém a mudança de sinal. Isso assegura que a raiz nunca saia do intervalo e, portanto, que os métodos convergem.

**A velocidade média é a mesma nos dois casos?**
Sim, $U \approx 0{,}356$ m/s, pois $H$ é o mesmo.

## Conclusão

Os dois métodos encontraram a profundidade $H \approx 0{,}702$ m, com velocidade média $U \approx 0{,}356$ m/s. Ambos têm convergência garantida, por serem intervalares. Para este problema, a falsa posição foi bem mais eficiente (7 iterações contra 22), embora, em geral, ela possa ficar lenta se um dos extremos permanecer fixo e a função tiver forte curvatura. A bisseção, por sua vez, tem velocidade previsível: o erro máximo cai pela metade a cada passo.
