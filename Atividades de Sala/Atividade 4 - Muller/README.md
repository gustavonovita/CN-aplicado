# Profundidade Crítica em Canal Trapezoidal pelo Método de Müller

Atividade prática de **Cálculo Numérico Aplicado**. O programa `Muller.py` determina a
**profundidade crítica** $H_c$ de um canal aberto de seção trapezoidal usando o
**método de Müller**, valida o resultado para $Q = 20\ \text{m}^3/\text{s}$ e, em seguida,
repete o cálculo para **oito vazões diferentes**, gerando um **arquivo de resultados** e um
**gráfico** $H_c \times Q$.

## Sumário

1. [Problema físico](#1-problema-físico)
2. [Formulação matemática](#2-formulação-matemática)
3. [O método de Müller](#3-o-método-de-müller)
4. [Estrutura do repositório](#4-estrutura-do-repositório)
5. [Como executar](#5-como-executar)
6. [Explicação do código](#6-explicação-do-código)
7. [Validação para Q = 20 m³/s](#7-validação-para-q--20-m³s)
8. [Varredura de múltiplas vazões](#8-varredura-de-múltiplas-vazões)
9. [Geração do gráfico Hc × Q](#9-geração-do-gráfico-hc--q)
10. [Interpretação física dos resultados](#10-interpretação-física-dos-resultados)
11. [Limitações e possíveis melhorias](#11-limitações-e-possíveis-melhorias)

---

## 1. Problema físico

Considere água escoando em um canal aberto de seção **trapezoidal**, com:

| Símbolo | Significado | Valor |
|---|---|---|
| $b$ | largura de fundo | 4,0 m |
| $z$ | inclinação lateral $z:1$ (avanço horizontal por unidade de altura) | 1,5 |
| $g$ | aceleração da gravidade | 9,81 m/s² |
| $H$ | profundidade da água | incógnita |
| $Q$ | vazão volumétrica | dada |

A geometria da seção molhada fornece:

$$A(H) = bH + zH^2 \quad \text{(área molhada)}, \qquad T(H) = b + 2zH \quad \text{(largura da superfície livre)}$$

O regime de escoamento é caracterizado pelo **número de Froude**:

$$Fr = \frac{U}{\sqrt{g D_h}}, \qquad U = \frac{Q}{A}, \qquad D_h = \frac{A}{T}$$

O escoamento é **crítico** quando $Fr = 1$, ou seja, quando os efeitos inerciais e
gravitacionais se equilibram. Isso equivale a

$$Fr^2 = \frac{Q^2\,T}{g\,A^3} = 1$$

A profundidade que satisfaz essa condição é a **profundidade crítica** $H_c$.

## 2. Formulação matemática

Substituindo $A(H)$ e $T(H)$ na condição de criticidade e reorganizando:

$$g\,\big[bH + zH^2\big]^3 - Q^2\,(b + 2zH) = 0$$

Como $bH + zH^2 = H(b + zH)$, expandimos $(b+zH)^3$ pelo binômio de Newton e obtemos um
**polinômio de sexta ordem** em $H$:

$$\boxed{\,g z^3 H^6 + 3 g b z^2 H^5 + 3 g b^2 z H^4 + g b^3 H^3 - 2 Q^2 z H - Q^2 b = 0\,}$$

Os coeficientes da forma geral $a_6H^6 + a_5H^5 + \dots + a_0$ são:

| Coeficiente | Expressão | Valor para $Q = 20$ |
|---|---|---|
| $a_6$ | $g z^3$ | 33,10875 |
| $a_5$ | $3 g b z^2$ | 264,87 |
| $a_4$ | $3 g b^2 z$ | 706,32 |
| $a_3$ | $g b^3$ | 627,84 |
| $a_2$ | $0$ | 0 |
| $a_1$ | $-2 Q^2 z$ | −1200 |
| $a_0$ | $-Q^2 b$ | −1600 |

> **Por que sexta ordem?** $A$ já é quadrática em $H$, então $A^3$ é de grau $2 \times 3 = 6$.
> O termo $Q^2 T$ é apenas linear e não altera o grau. O coeficiente de $H^2$ é nulo porque
> a expansão de $H^3(b+zH)^3$ só gera potências de $H^3$ a $H^6$, e $Q^2(b+2zH)$ só gera
> $H^0$ e $H^1$.

> **Admissibilidade física.** Um polinômio de grau 6 tem até seis raízes (reais ou complexas).
> Como $H$ é uma profundidade, só interessa a raiz **real e positiva**. Para $Q = 20$ m³/s,
> as seis raízes são $1{,}17105$; $-1{,}48260$; $-0{,}42794 \pm 1{,}42207i$ e
> $-3{,}41629 \pm 0{,}97465i$, ou seja, apenas **uma** é fisicamente admissível.

## 3. O método de Müller

O método de Müller é uma generalização do método da secante. Em vez de ajustar uma **reta**
a dois pontos, ele ajusta uma **parábola** a **três** pontos $(H_0, f_0)$, $(H_1, f_1)$,
$(H_2, f_2)$, e usa a raiz da parábola mais próxima de $H_2$ como nova aproximação.

**Passo 1 — Avaliar a função nos três pontos:**

$$f_0 = f(H_0), \quad f_1 = f(H_1), \quad f_2 = f(H_2)$$

**Passo 2 — Diferenças e inclinações das secantes:**

$$h_0 = H_1 - H_0, \qquad h_1 = H_2 - H_1, \qquad \delta_0 = \frac{f_1 - f_0}{h_0}, \qquad \delta_1 = \frac{f_2 - f_1}{h_1}$$

**Passo 3 — Coeficientes da parábola** (centrada em $H_2$):

$$a = \frac{\delta_1 - \delta_0}{h_1 + h_0}, \qquad b_M = a\,h_1 + \delta_1, \qquad c = f_2$$

(O índice $M$ em $b_M$ evita confusão com a largura de fundo $b$ do canal.)

**Passo 4 — Discriminante e denominador:**

$$D = \sqrt{b_M^2 - 4ac}, \qquad E = \begin{cases} b_M + D, & |b_M + D| > |b_M - D| \\[2pt] b_M - D, & \text{caso contrário} \end{cases}$$

A escolha do sinal que **maximiza $|E|$** evita o *cancelamento numérico* (subtração de dois
números quase iguais) e seleciona a raiz da parábola mais próxima de $H_2$.

**Passo 5 — Nova aproximação:**

$$\Delta H = -\frac{2c}{E}, \qquad H_3 = H_2 + \Delta H$$

**Passo 6 — Atualização e critério de parada:**

$$H_0 \leftarrow H_1, \quad H_1 \leftarrow H_2, \quad H_2 \leftarrow H_3, \qquad \text{parar quando } |H_3 - H_2| < 10^{-8}\ \text{m}$$

> O discriminante pode ser **negativo** nas iterações intermediárias. Por isso o código usa
> aritmética complexa (`np.emath.sqrt`). Esta é uma vantagem do método: ele pode "passar" pelo
> plano complexo e retornar ao eixo real.

## 4. Estrutura do repositório

```
.
├── Muller.py                # código-fonte completo
├── resultados_Q_Hc.txt      # saída: tabela Q × Hc (gerado ao executar)
├── grafico_Hc_x_Q.png       # saída: gráfico Hc × Q (gerado ao executar)
└── README.md
```

## 5. Como executar

**Requisitos:** Python 3.8+, `numpy` e `matplotlib`.

```bash
pip install numpy matplotlib
python Muller.py
```

Ao rodar, o programa:

1. imprime o histórico de iterações para $Q = 20$ m³/s (validação);
2. calcula $H_c$ para as 8 vazões e imprime cada resultado;
3. grava `resultados_Q_Hc.txt` e `grafico_Hc_x_Q.png` na pasta de execução.

## 6. Explicação do código

### 6.1 Parâmetros globais

```python
b = 4.0        # largura de fundo (m)
z = 1.5        # inclinação lateral
g = 9.81       # gravidade (m/s²)
max_iter = 100 # limite de iterações (evita laço infinito)
tol = 1e-8     # tolerância do critério de parada (m)
```

A vazão **não** é global: ela é passada como parâmetro das funções. Isso é o que permite
reaproveitar o mesmo código para várias vazões.

### 6.2 Função `polinomio_f(H, Q)`

```python
def polinomio_f(H, Q):
    a6 = g * (z**3) * (H**6)
    a5 = 3 * g * b * (z**2) * (H**5)
    a4 = 3 * g * (b**2) * z * (H**4)
    a3 = g * (b**3) * (H**3)
    a1 = -2 * (Q**2) * z * H
    a0 = -(Q**2) * b
    return a1 + a3 + a4 + a5 + a6 + a0
```

Implementa $f(H)$ termo a termo. O termo em $H^2$ (`a2`) não aparece porque é zero. **O termo
`a0` precisa estar na soma do `return`**: sem ele o polinômio ganharia uma raiz espúria e o
método convergiria para um valor errado.

### 6.3 Função `metodo_muller(Q, H0_val, H1_val, H2_val, imprimir_passos=False)`

Recebe a vazão, os três chutes iniciais e uma opção para imprimir o histórico. A cada
iteração `k`:

| Trecho do código | O que faz |
|---|---|
| `f0, f1, f2 = polinomio_f(...)` | avalia $f$ nos três pontos |
| `h0, h1`, `delta_0, delta_1` | diferenças e inclinações das secantes |
| `a`, `b_m`, `c` | coeficientes da parábola interpoladora |
| `D = np.emath.sqrt(b_m**2 - 4*a*c)` | discriminante; devolve complexo se negativo |
| `if abs(b_m + D) > abs(b_m - D)` | escolhe o denominador `E` de maior módulo |
| `Delta_H = -(2*c/E)`; `H3 = H2 + Delta_H` | nova aproximação da raiz |
| `erro = abs(H3 - H2)` | erro entre iterações consecutivas |
| `if erro < tol: return H3.real` | critério de convergência; devolve a parte real |
| `H0, H1, H2 = H1, H2, H3` | desloca a janela de três pontos |

Se `imprimir_passos=True`, cada iteração imprime: número `k`, as três aproximações
($H_0, H_1, H_2$), a nova aproximação $H_3$, o valor $f(H_3)$ e o erro $|H_3 - H_2|$.

O retorno é sempre `H3.real`: se a raiz convergida for real, a parte imaginária é
numericamente nula, e devolver apenas a parte real mantém o resultado como `float` (necessário
para armazenar em lista, gravar em arquivo e plotar).

## 7. Validação para Q = 20 m³/s

Chutes iniciais: $H_0 = 0{,}5$ m, $H_1 = 1{,}0$ m, $H_2 = 1{,}5$ m.

| k | $H_0$ | $H_1$ | $H_2$ | $H_3$ | $f(H_3)$ | $\lvert H_3 - H_2 \rvert$ |
|---|---|---|---|---|---|---|
| 1 | 0,5000 | 1,0000 | 1,5000 | 1,142989 | −2,3804e+02 | 3,5701e−01 |
| 2 | 1,0000 | 1,5000 | 1,1430 | 1,169522 | −1,3450e+01 | 2,6532e−02 |
| 3 | 1,5000 | 1,1430 | 1,1695 | 1,171030 | −1,3306e−01 | 1,5085e−03 |
| 4 | 1,1430 | 1,1695 | 1,1710 | 1,171045 | 5,4932e−06 | 1,5039e−05 |
| 5 | 1,1695 | 1,1710 | 1,1710 | **1,171045** | 0,0000e+00 | 6,2083e−10 |

**Resultado: $H_c = 1{,}171045$ m, em 5 iterações.**

O erro cai de $10^{-1} \to 10^{-2} \to 10^{-3} \to 10^{-5} \to 10^{-10}$, comportamento típico
de convergência **superlinear** (ordem ≈ 1,84, contra ≈ 1,62 da secante).

**Verificação física:** com $H_c = 1{,}171045$ m, $A_c = 6{,}7412$ m² e $T_c = 7{,}5131$ m, o
número de Froude resulta em $Fr^2 = Q^2 T_c / (g A_c^3) \approx 1{,}000000$, confirmando que
a raiz corresponde à condição de escoamento crítico.

## 8. Varredura de múltiplas vazões

Depois de validar o método para $Q = 20$ m³/s, o mesmo código é usado para encontrar a
profundidade crítica de **oito vazões**: 5, 10, 15, 20, 25, 30, 35 e 40 m³/s.

### 8.1 A ideia

Cada vazão $Q$ define um **polinômio de sexta ordem diferente**. Os coeficientes $a_1$ e $a_0$
dependem de $Q^2$:

$$a_1 = -2Q^2 z, \qquad a_0 = -Q^2 b$$

enquanto $a_6, a_5, a_4, a_3$ dependem só da geometria e permanecem fixos. Portanto, para cada
$Q$ é preciso resolver uma equação nova. Isso só é possível sem copiar código porque `Q` é um
**parâmetro** de `polinomio_f` e de `metodo_muller`.

### 8.2 O trecho do código

```python
vazoes = np.array([5, 10, 15, 20, 25, 30, 35, 40], dtype=float)
profundidades_criticas = []

print("--- Calculando Profundidades Críticas (Varredura) ---")
for Q_atual in vazoes:
    Hc = metodo_muller(Q_atual, 0.5, 1.0, 1.5, imprimir_passos=False)
    profundidades_criticas.append(Hc)
    print(f"Vazão Q = {Q_atual:02.0f} m³/s -> Hc = {Hc:.6f} m")
```

### 8.3 Passo a passo

1. **`vazoes = np.array([...], dtype=float)`** — cria um vetor NumPy com as vazões. O
   `dtype=float` garante números de ponto flutuante, evitando divisões inteiras e permitindo
   `Q**2` sem estouro de tipo.
2. **`profundidades_criticas = []`** — lista vazia que vai guardar um $H_c$ por vazão, na mesma
   ordem do vetor `vazoes`.
3. **`for Q_atual in vazoes:`** — percorre as vazões uma a uma.
4. **`metodo_muller(Q_atual, 0.5, 1.0, 1.5, ...)`** — resolve o polinômio daquela vazão. Os
   **três chutes iniciais são os mesmos em todas as vazões** (0,5; 1,0; 1,5 m). Eles ficam
   dentro de uma faixa de profundidades plausíveis, e como a função recebe os chutes por
   argumento, cada chamada começa "do zero", sem herdar valores da vazão anterior.
5. **`imprimir_passos=False`** — a chamada é silenciosa, para não poluir a tela com o histórico
   das 8 vazões. Só o resultado final de cada vazão é impresso.
6. **`profundidades_criticas.append(Hc)`** — armazena o resultado.
7. **`print(f"... {Q_atual:02.0f} ... {Hc:.6f}")`** — formata a vazão com 2 dígitos e zero
   à esquerda (`05`, `10`, ...) e $H_c$ com 6 casas decimais.

### 8.4 Salvando os resultados em arquivo

```python
dados_exportacao = np.column_stack((vazoes, profundidades_criticas))
np.savetxt("resultados_Q_Hc.txt", dados_exportacao,
           fmt="%.2f %.6f", header="Q(m3/s) Hc(m)")
```

- **`np.column_stack((vazoes, profundidades_criticas))`** junta os dois vetores como **colunas**
  de uma matriz $8 \times 2$: coluna 0 = $Q$, coluna 1 = $H_c$.
- **`np.savetxt(...)`** grava a matriz em texto:
  - `fmt="%.2f %.6f"` define o formato de cada linha: a 1ª coluna com 2 casas decimais e a 2ª
    com 6, separadas por espaço;
  - `header="Q(m3/s) Hc(m)"` escreve um cabeçalho. O NumPy acrescenta automaticamente `# ` no
    início, então a linha vira um comentário.

Conteúdo de `resultados_Q_Hc.txt`:

```
# Q(m3/s) Hc(m)
5.00 0.507112
10.00 0.776569
15.00 0.989536
20.00 1.171045
25.00 1.331607
30.00 1.476888
35.00 1.610373
40.00 1.734394
```

Para ler o arquivo de volta em Python (o `#` é ignorado automaticamente):

```python
import numpy as np
Q, Hc = np.loadtxt("resultados_Q_Hc.txt", unpack=True)
```

### 8.5 Resultados

| $Q$ (m³/s) | $H_c$ (m) | $\Delta H_c$ em relação à vazão anterior (m) | $Fr^2$ ¹ |
|---|---|---|---|
| 5 | 0,507112 | — | 0,9999996 |
| 10 | 0,776569 | +0,269457 | 1,0000010 |
| 15 | 0,989536 | +0,212967 | 1,0000004 |
| 20 | 1,171045 | +0,181509 | 1,0000014 |
| 25 | 1,331607 | +0,160562 | 1,0000005 |
| 30 | 1,476888 | +0,145281 | 0,9999991 |
| 35 | 1,610373 | +0,133485 | 0,9999995 |
| 40 | 1,734394 | +0,124021 | 0,9999997 |

¹ Calculado com os valores de $H_c$ arredondados em 6 casas decimais (os mesmos do arquivo). A
pequena diferença em relação a 1 (da ordem de $10^{-6}$) vem desse arredondamento, não do
método.

## 9. Geração do gráfico Hc × Q

### 9.1 O trecho do código

```python
plt.figure(figsize=(8, 5))
plt.plot(vazoes, profundidades_criticas, marker='o', linestyle='-', color='b')
plt.title("Variação da Profundidade Crítica com a Vazão")
plt.xlabel("Vazão Volumétrica Q (m³/s)")
plt.ylabel("Profundidade Crítica Hc (m)")
plt.grid(True)
plt.savefig("grafico_Hc_x_Q.png")
```

### 9.2 Linha por linha

| Linha | O que faz |
|---|---|
| `plt.figure(figsize=(8, 5))` | cria uma nova figura de 8 × 5 polegadas (largura × altura) |
| `plt.plot(vazoes, profundidades_criticas, ...)` | desenha a curva: eixo **x** = vazões, eixo **y** = profundidades críticas |
| `marker='o'` | coloca um círculo em cada ponto calculado (8 pontos), mostrando onde estão os dados reais |
| `linestyle='-'` | liga os pontos com uma linha contínua |
| `color='b'` | cor azul (*blue*) |
| `plt.title(...)` | título do gráfico |
| `plt.xlabel(...)`, `plt.ylabel(...)` | rótulos dos eixos, com as unidades |
| `plt.grid(True)` | adiciona uma grade para facilitar a leitura dos valores |
| `plt.savefig("grafico_Hc_x_Q.png")` | grava a figura em um arquivo PNG na pasta de execução |

`vazoes` (vetor NumPy) e `profundidades_criticas` (lista) têm o mesmo comprimento, 8, como
exige o `plt.plot`. O que está no eixo x é a **variável independente** ($Q$, que escolhemos), e
no eixo y a **variável dependente** ($H_c$, resultado do método de Müller).

### 9.3 Resultado

![Gráfico Hc x Q](grafico_Hc_x_Q.png)

### 9.4 Observações sobre o `plt.savefig`

- O código **salva** o gráfico mas não o **exibe** na tela. Para visualizar durante a execução,
  adicione `plt.show()` **depois** do `plt.savefig(...)`. Se `plt.show()` vier antes, em muitos
  ambientes a figura é fechada ao fechar a janela e o arquivo salvo sai vazio.
- Para um PNG de maior resolução (por exemplo, para relatório), use
  `plt.savefig("grafico_Hc_x_Q.png", dpi=300, bbox_inches="tight")`.
- O arquivo é gravado na **pasta de trabalho atual** (de onde o comando `python` foi
  executado), e não necessariamente na pasta do script.

## 10. Interpretação física dos resultados

**Como $H_c$ varia com $Q$?** A profundidade crítica **cresce** com a vazão, mas com
**concavidade para baixo**: cada aumento de 5 m³/s eleva $H_c$ cada vez menos (+0,269 m entre
5 e 10 m³/s, mas apenas +0,124 m entre 35 e 40 m³/s).

**Isso é fisicamente coerente?** Sim.

- Mais vazão exige mais profundidade para que o escoamento atinja o equilíbrio $Fr = 1$.
- O crescimento é **sublinear** porque a seção trapezoidal se **alarga** com a profundidade
  ($T = b + 2zH$). Parte da vazão extra é acomodada pelo aumento da largura da superfície, e
  não só da profundidade.
- Em um canal retangular, $H_c = (Q^2/(g b^2))^{1/3} \propto Q^{2/3}$, que também é côncava.
  O trapézio mostra o mesmo comportamento qualitativo.

**Sensibilidade aos chutes iniciais.** O método de Müller é *local*: chutes razoáveis
($0 < H < 3$ m) convergem para a raiz física $1{,}17105$ m, enquanto chutes muito distantes ou
negativos podem convergir para as outras raízes do polinômio (negativas ou complexas), que não
têm significado físico. Por isso é recomendável verificar sempre que $H_c > 0$, $H_c$ é real e
$Fr \approx 1$.

## 11. Limitações e possíveis melhorias

- **Divisão por zero:** se duas aproximações consecutivas coincidirem ($h_0 = 0$ ou $h_1 = 0$),
  ou se $E = 0$, ocorre `ZeroDivisionError`. Pode-se proteger com um `if` antes das divisões.
- **Sem aviso de não convergência:** se `max_iter` for atingido, a função retorna o último
  valor sem avisar. Pode-se imprimir uma mensagem ou lançar uma exceção nesse caso.
- **Verificação de $Fr$ e da admissibilidade:** o programa não calcula automaticamente
  $Fr^2$ nem testa $H_c > 0$ e $H_c \in \mathbb{R}$. Esses cálculos podem ser adicionados após cada
  convergência.
- **Contagem de iterações:** `metodo_muller` poderia retornar também o número de iterações
  `k` para ser registrado junto com cada resultado.
- **Teste de chutes iniciais diferentes** (ex.: $(0{,}1;\,0{,}5;\,1{,}0)$ e $(1{,}0;\,2{,}0;\,3{,}0)$)
  para verificar se o método converge sempre para a mesma raiz.

## Autor

<!-- Seu nome / GitHub / disciplina / instituição -->
