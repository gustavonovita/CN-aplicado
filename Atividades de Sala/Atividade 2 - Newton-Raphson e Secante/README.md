# Temperatura de Equilíbrio de uma Superfície Aquecida
### Métodos de Newton–Raphson e da Secante em Python

Atividade prática de **Cálculo Numérico Aplicado**. O objetivo é determinar numericamente a temperatura de equilíbrio de uma superfície que recebe uma potência térmica constante e perde calor para o ambiente por **convecção** e por **radiação**, usando dois métodos iterativos de busca de raízes:

- **Newton–Raphson**, que usa a derivada analítica da função;
- **Secante**, que usa apenas valores da função (como se ela fosse uma "caixa-preta").

---

## Sumário

1. [Problema físico](#1-problema-físico)
2. [Formulação matemática](#2-formulação-matemática)
3. [Método de Newton–Raphson](#3-método-de-newtonraphson)
4. [Método da Secante](#4-método-da-secante)
5. [Estrutura do repositório](#5-estrutura-do-repositório)
6. [Como executar](#6-como-executar)
7. [Resultados](#7-resultados)
8. [Comparação dos métodos](#8-comparação-dos-métodos)
9. [Observações sobre as implementações](#9-observações-sobre-as-implementações)

---

## 1. Problema físico

Uma pequena superfície de área $A$ recebe uma potência térmica constante $P$. Ela perde calor para o ambiente ao mesmo tempo por convecção e por radiação. Em regime permanente, tudo o que entra é igual ao que sai:

$$P = hA\,(T - T_\infty) + \varepsilon\sigma A\,(T^4 - T_\infty^4)$$

em que $T$ é a temperatura de equilíbrio da superfície (a incógnita).

| Símbolo | Significado | Valor | Variável no código |
|---|---|---|---|
| $A$ | Área da superfície | 0,01 m² | `A` |
| $h$ | Coeficiente de convecção | 10 W/(m²·K) | `h` |
| $\varepsilon$ | Emissividade da superfície | 0,8 | `epsilon` |
| $T_\infty$ | Temperatura do ambiente | 300 K | `T_amb` |
| $P$ | Potência fornecida | 20 W | `P` |
| $\sigma$ | Constante de Stefan–Boltzmann | 5,67 × 10⁻⁸ W/(m²·K⁴) | `sigma` |

## 2. Formulação matemática

### Problema como zero de função

Passando todos os termos para um lado, o problema vira "encontrar $T$ tal que $f(T) = 0$":

$$f(T) = hA\,(T - T_\infty) + \varepsilon\sigma A\,(T^4 - T_\infty^4) - P$$

$f(T)$ é o **resíduo do balanço térmico**, em watts:

- $f(T) > 0$: a superfície perde mais calor do que recebe (temperatura acima do equilíbrio);
- $f(T) < 0$: a superfície recebe mais calor do que perde (temperatura abaixo do equilíbrio);
- $f(T) = 0$: equilíbrio.

### Derivada analítica

Derivando termo a termo em relação a $T$:

$$f'(T) = hA + 4\,\varepsilon\sigma A\,T^3$$

Como $f'(T) > 0$ para qualquer $T > 0$, a função é estritamente crescente e existe **uma única raiz** positiva. Por isso os métodos não correm o risco de encontrar uma raiz "errada" nem de encontrar derivada nula no caminho.

### Critério de parada

Os dois métodos usam o mesmo critério:

$$|T_{k+1} - T_k| < 10^{-6}$$

com um número máximo de iterações (`max_iter = 100`) para evitar que o programa fique preso em um laço infinito caso não haja convergência.

## 3. Método de Newton–Raphson

**Ideia geométrica:** em cada passo, traça-se a reta **tangente** à curva $f$ no ponto atual e usa-se o ponto onde essa reta cruza o eixo $T$ como a nova aproximação.

$$T_{k+1} = T_k - \frac{f(T_k)}{f'(T_k)}$$

Arquivo: [`newton_raphson.py`](newton_raphson.py)

Pontos principais do código:

- `f(T)` implementa o resíduo do balanço térmico.
- `df(T)` implementa a derivada analítica $f'(T) = hA + 4\varepsilon\sigma A T^3$.
- O laço `for k in range(1, max_iter + 1)` faz as iterações. Na primeira, `Tk` recebe o chute inicial `T0`; nas seguintes, recebe o valor calculado na iteração anterior (`Tk_novo`).
- A atualização é `Tk_novo = Tk - f(Tk) / df(Tk)`.
- O erro é `abs(Tk_novo - Tk)`, ou seja, quanto a temperatura mudou naquele passo. Se for menor que `tol`, o programa imprime a raiz e encerra com `break`.
- Chute inicial usado no código: `T0 = 300` K.

## 4. Método da Secante

**Ideia geométrica:** é o Newton sem derivada. Em vez da tangente, usa-se a reta **secante**, que passa pelos dois últimos pontos calculados. A inclinação dessa reta aproxima a derivada:

$$f'(T_k) \approx \frac{f(T_k) - f(T_{k-1})}{T_k - T_{k-1}}$$

Substituindo na fórmula de Newton:

$$T_{k+1} = T_k - f(T_k)\,\frac{T_k - T_{k-1}}{f(T_k) - f(T_{k-1})}$$

Arquivo: [`secante.py`](secante.py)

Pontos principais do código:

- Aqui só `f(T)` é usada. **A derivada não aparece em nenhum momento**, como exige a atividade (a função é tratada como caixa-preta).
- A secante precisa de dois pontos iniciais. Eles foram escolhidos a partir da tabela de resíduos fornecida na atividade, em lados opostos da raiz: `T0 = 400` K, com $f(400) = -2{,}062$ W, e `T1 = 425` K, com $f(425) = +3{,}625$ W. A mudança de sinal indica que a raiz está entre eles.
- Em cada iteração, `Tk_novo = T1 - f(T1) * ((T1 - T0) / (f(T1) - f(T0)))`.
- Depois de testar a convergência, o código "desliza" a janela de pontos: `T0 = T1` e `T1 = Tk_novo`.

**Tabela de referência da atividade** (resíduos de um "experimento"):

| T (K) | f(T) (W) |
|---|---|
| 350 | −11,867325 |
| 375 | −7,204043 |
| 400 | −2,062000 |
| 425 | 3,624717 |
| 450 | 9,926275 |
| 475 | 16,917092 |
| 500 | 24,675840 |

## 5. Estrutura do repositório

```
.
├── newton_raphson.py   # Parte 1: método de Newton–Raphson
├── secante.py          # Parte 2: método da Secante
└── README.md           # Este arquivo
```

Os dois scripts são independentes: cada um define seus próprios parâmetros e funções e pode ser executado separadamente.

## 6. Como executar

**Requisitos:** Python 3.6 ou superior (usa f-strings). Não há dependências externas.

```bash
git clone https://github.com/<seu-usuario>/<nome-do-repositorio>.git
cd <nome-do-repositorio>

python3 newton_raphson.py
python3 secante.py
```

A cada iteração, os programas imprimem o número da iteração, a temperatura, o valor de $f$ e o erro usado no critério de parada. Ao final, imprimem a temperatura de equilíbrio e o número de iterações.

## 7. Resultados

### Newton–Raphson (`T0 = 300` K)

Saída real do programa:

```
Iteração: 1, T1 = 300.000000, f(T1) = -20.000000, erro = 134.238278
Iteração: 2, T2 = 434.238278, f(T2) = 5.877852, erro = 23.647113
Iteração: 3, T3 = 410.591165, f(T3) = 0.276694, erro = 1.226523
Iteração: 4, T4 = 409.364642, f(T4) = 0.000689, erro = 0.003069
Iteração: 5, T5 = 409.361573, f(T5) = 0.000000, erro = 0.000000
Raiz encontrada: T = 409.361573 após 5 iterações.
```

| k | $T_k$ (K) | $f(T_k)$ (W) | erro $\lvert T_{k+1}-T_k\rvert$ |
|---|---|---|---|
| 1 | 300,000000 | −20,000000 | 134,238278 |
| 2 | 434,238278 | 5,877852 | 23,647113 |
| 3 | 410,591165 | 0,276694 | 1,226523 |
| 4 | 409,364642 | 0,000689 | 0,003069 |
| 5 | 409,361573 | 0,000000 | 0,000000 |

Observe como o erro despenca rapidamente perto da raiz (134 → 23,6 → 1,2 → 0,003 → ~0). Isso é a **convergência quadrática** do método de Newton: perto da solução, o número de dígitos corretos aproximadamente dobra a cada iteração.

### Secante (`T0 = 400` K, `T1 = 425` K)

Saída real do programa:

```
Iteração: 1, Temperatura = 409.064984, f(Tk) = 3.624717, Erro = 15.935016
Iteração: 2, Temperatura = 409.352211, f(Tk) = -0.066534, Erro = 0.287227
Iteração: 3, Temperatura = 409.361579, f(Tk) = -0.002101, Erro = 0.009368
Iteração: 4, Temperatura = 409.361573, f(Tk) = 0.000001, Erro = 0.000006
Iteração: 5, Temperatura = 409.361573, f(Tk) = -0.000000, Erro = 0.000000
Temperatura de Equilíbrio encontrada = 409.3615729509906
```

| k | Nova temperatura (K) | erro $\lvert T_{k+1}-T_k\rvert$ |
|---|---|---|
| 1 | 409,064984 | 15,935016 |
| 2 | 409,352211 | 0,287227 |
| 3 | 409,361579 | 0,009368 |
| 4 | 409,361573 | 0,000006 |
| 5 | 409,361573 | 0,000000 |

### Resultado final

$$\boxed{T_{\text{equilíbrio}} \approx 409{,}3616\ \text{K} \;(\approx 136{,}2\ ^\circ\text{C})}$$

**Verificação física:** nessa temperatura, a convecção dissipa ≈ 10,9 W e a radiação ≈ 9,1 W, somando os 20 W fornecidos. Os dois mecanismos têm peso comparável, por isso a radiação (termo em $T^4$) não pode ser desprezada.

## 8. Comparação dos métodos

**a) Os dois métodos convergiram para a mesma temperatura?**
Sim. Ambos chegaram a **409,361573 K**. Isso era esperado, pois a raiz é única ($f$ é estritamente crescente).

**b) Qual necessitou de menos iterações?**
Os dois necessitaram de **5 iterações**, ou seja, empate nesta execução. Vale notar que o Newton partiu de um chute mais distante da raiz (300 K), enquanto a secante partiu de dois pontos que já cercavam a raiz (400 K e 425 K). Em ordem de convergência, o Newton é **quadrático** (ordem 2) e a secante é **superlinear** (ordem ≈ 1,618, a razão áurea). Em custo por iteração, o Newton avalia $f$ e $f'$ a cada passo, enquanto a secante avalia apenas **uma** nova $f$, reaproveitando o valor anterior.

**c) Por que o Newton–Raphson necessita do conhecimento de $f'(T)$?**
Porque o passo $T_{k+1} = T_k - f(T_k)/f'(T_k)$ vem da linearização de $f$ em torno de $T_k$ (série de Taylor de 1ª ordem): $f(T) \approx f(T_k) + f'(T_k)(T - T_k)$. A derivada é a inclinação da reta tangente e determina o quanto "andar" para zerar a função. Sem $f'$, não há como calcular o passo.

**d) Qual a principal vantagem da secante quando $f(T)$ vem de um experimento, simulação ou programa externo?**
Ela precisa **apenas de valores de $f(T)$**, que é exatamente o que uma caixa-preta fornece. Não exige a expressão analítica da derivada, que pode ser desconhecida, trabalhosa de obter ou simplesmente inexistente (um experimento não devolve derivadas).

**e) Em que sentido a secante é uma modificação do Newton–Raphson?**
Os dois métodos têm a mesma estrutura: "aproximar $f$ por uma reta perto da raiz e ver onde ela cruza o eixo". A diferença está em **qual reta** se usa. O Newton usa a **tangente** (derivada exata em um ponto). A secante usa a **reta por dois pontos**, trocando a derivada por uma diferença finita. Quando $T_{k-1}$ se aproxima de $T_k$, a inclinação da secante tende à da tangente, e a secante tende ao Newton. Em resumo: **a secante é o Newton com $f'$ substituída por uma diferença finita de dois valores consecutivos de $f$**.

## 9. Observações sobre as implementações

- **Como ler a saída do Newton:** na iteração $k$, o programa imprime $T_k$ e $f(T_k)$ do ponto *atual*, e o erro $\lvert T_{k+1}-T_k\rvert$ do passo que leva ao próximo ponto. A temperatura final reportada é `Tk_novo`, o último valor calculado.
- **Como ler a saída da secante:** na iteração $k$, a coluna "Temperatura" mostra o valor *recém-calculado* (`Tk_novo`), enquanto a coluna `f(Tk)` mostra $f(T_1)$, ou seja, o valor de $f$ no ponto usado para dar o passo (o ponto anterior à nova temperatura). Por isso a primeira linha mostra $f(425) = 3{,}624717$, valor que coincide com a Tabela 1. Para imprimir o resíduo da nova aproximação, basta trocar `f(T1)` por `f(Tk_novo)` na instrução `print`.
- **Chute inicial do Newton:** o enunciado sugere $T_0 = 350$ K, e o código usa $T_0 = 300$ K. Ambos funcionam. Com $T_0 = 350$ K, o método também converge em 5 iterações para o mesmo valor (409,361573 K), com passos iniciais menores (erro 66,75 → 7,28 → 0,11 → 0,000025 → ~0). Com $T_0 = 300$ K, o primeiro passo é maior (134 K), pois a derivada é pequena em temperaturas baixas, mas o método se recupera rápido.
- **Segurança na secante:** a fórmula divide por $f(T_1) - f(T_0)$. Neste problema isso nunca é zero (a função é estritamente crescente), mas em problemas gerais convém verificar essa divisão por zero.
- **Critério de parada:** $\lvert T_{k+1}-T_k\rvert < 10^{-6}$ mede o tamanho do último passo, e não o erro verdadeiro em relação à raiz exata. Funciona bem aqui por causa da convergência rápida, e olhar também o valor de $f(T_k)$ (que chega a ≈ 0) ajuda a confirmar que a solução é boa.

---

## Autor

<!-- Preencha com seus dados -->
**Nome:** _seu nome_
**Disciplina:** Cálculo Numérico Aplicado
**Instituição:** _sua instituição_
