# Condução de Calor em uma Placa 2D com o Método de Gauss–Seidel

Atividade prática da disciplina **Cálculo Numérico Aplicado**.

Este repositório contém um programa em **Python** que calcula a distribuição de temperatura, em regime permanente, no interior de uma placa quadrada. Ele resolve o sistema linear que surge da discretização da equação de Laplace usando o **método iterativo de Gauss–Seidel**, sem nunca montar a matriz de coeficientes.

Este README foi escrito como material de consulta: explica o problema físico, toda a teoria numérica envolvida, o código linha a linha, os resultados e como interpretá-los.

---

## Sumário

1. [Visão geral](#1-visão-geral)
2. [Estrutura do repositório e como executar](#2-estrutura-do-repositório-e-como-executar)
3. [O problema físico](#3-o-problema-físico)
4. [Fundamentação teórica](#4-fundamentação-teórica)
   - 4.1 [A equação de Laplace](#41-a-equação-de-laplace)
   - 4.2 [Discretização por diferenças finitas](#42-discretização-por-diferenças-finitas)
   - 4.3 [O sistema linear resultante](#43-o-sistema-linear-resultante)
   - 4.4 [Métodos iterativos: Jacobi e Gauss–Seidel](#44-métodos-iterativos-jacobi-e-gaussseidel)
   - 4.5 [Gauss–Seidel aplicado à placa](#45-gaussseidel-aplicado-à-placa)
   - 4.6 [Convergência](#46-convergência)
   - 4.7 [Critério de parada × resíduo](#47-critério-de-parada--resíduo)
5. [Explicação detalhada do código](#5-explicação-detalhada-do-código)
6. [Resultados](#6-resultados)
7. [Análise dos resultados](#7-análise-dos-resultados)
8. [Questões teóricas respondidas](#8-questões-teóricas-respondidas)
9. [Possíveis melhorias e adaptações](#9-possíveis-melhorias-e-adaptações)
10. [Glossário](#10-glossário)

---

## 1. Visão geral

| Item | Descrição |
|------|-----------|
| **Problema** | Temperatura em regime permanente numa placa quadrada com bordas a temperaturas fixas |
| **Equação** | Equação de Laplace em 2D: $\partial^2 T/\partial x^2 + \partial^2 T/\partial y^2 = 0$ |
| **Discretização** | Diferenças finitas centrais, malha uniforme ($\Delta x = \Delta y$) |
| **Método numérico** | Gauss–Seidel (iterativo, *matrix-free*) |
| **Malhas testadas** | $n = 5$, $n = 10$ e $n = 20$ nós internos por direção |
| **Critério de parada** | Maior variação de temperatura entre duas varreduras $< 10^{-6}$ °C |
| **Bibliotecas** | NumPy (matriz) e Matplotlib (gráfico) |

---

## 2. Estrutura do repositório e como executar

```
.
├── gauss_seidel.py            # código-fonte (o programa completo)
├── grafico_calor_placa.png    # gráfico gerado pelo programa (n = 20)
└── README.md                  # este arquivo
```

### Requisitos

- Python 3.8 ou superior
- NumPy e Matplotlib

```bash
pip install numpy matplotlib
```

### Execução

```bash
python gauss_seidel.py
```

O programa imprime no terminal, para cada malha, o número de iterações, o erro final e o resíduo máximo (e a temperatura do nó central para $n = 5$). Ao final, grava o arquivo `grafico_calor_placa.png` na pasta atual.

### Saída esperada no terminal

```
--- Simulação para n = 5 ---
Iterações realizadas: 56
Erro final (critério de parada): 8.4928e-07
Temperatura do nó central T(4,4) [Python T[3,3]]: 40.0000 ºC
Resíduo máximo no sistema linear (max |ri,j|): 1.6986e-06

--- Simulação para n = 10 ---
Iterações realizadas: 178
Erro final: 9.6961e-07
Resíduo máximo: 1.8607e-06

--- Simulação para n = 20 ---
Iterações realizadas: 598
Erro final: 9.7862e-07
Resíduo máximo: 1.9354e-06

O gráfico de preenchimento térmico foi salvo como 'grafico_calor_placa.png'.
```

---

## 3. O problema físico

Uma placa quadrada de lado $L = 1\ \text{m}$ conduz calor em **regime permanente** (as temperaturas não mudam mais com o tempo) e **sem geração interna** de calor (não há fontes nem sumidouros dentro da placa). As temperaturas nas bordas são mantidas constantes:

| Borda | Temperatura |
|-------|-------------|
| Inferior | 100 °C |
| Superior | 20 °C |
| Esquerda | 20 °C |
| Direita | 20 °C |

```
              T = 20 °C
        ┌─────────────────────┐
        │                     │
 T = 20 │       placa         │ T = 20
  °C    │    (nós internos    │  °C
        │     desconhecidos)  │
        │                     │
        └─────────────────────┘
             T = 100 °C
```

**O que se quer calcular:** a temperatura em cada ponto (nó) **interno** da placa. Os valores das bordas já são conhecidos, pois são impostos pelo problema.

**Intuição física:** o calor entra pela borda inferior quente e flui para cima e para os lados, escoando pelas três bordas frias. Em regime permanente, a placa fica mais quente perto da borda de baixo e mais fria perto das outras três.

---

## 4. Fundamentação teórica

### 4.1 A equação de Laplace

A condução de calor em regime permanente, sem fontes internas e com condutividade térmica constante, é governada pela **equação de Laplace** em duas dimensões:

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

Ela diz que a curvatura da temperatura numa direção é compensada pela curvatura na outra, de modo que o campo de temperatura é "o mais suave possível". Uma consequência importante é o **princípio do máximo**: a temperatura em qualquer ponto interno fica **entre** o menor e o maior valor das bordas. Neste problema, todas as temperaturas internas ficam entre 20 °C e 100 °C.

Como a equação é uma equação diferencial parcial (EDP) com solução analítica difícil para casos gerais de contorno, usa-se um método numérico: **discretiza-se** a placa numa malha de pontos e transforma-se a EDP em um sistema de equações algébricas.

### 4.2 Discretização por diferenças finitas

#### A malha

A placa é dividida numa malha uniforme com **$n$ nós internos em cada direção**. Incluindo as bordas, há $n + 2$ pontos por direção. O espaçamento entre pontos é

$$h = \Delta x = \Delta y = \frac{L}{n+1}$$

Cada ponto é identificado por um par de índices $(i, j)$, onde $i$ percorre a direção horizontal ($x$) e $j$ a vertical ($y$).

```
n = 3 (5 x 5 pontos no total):

   j=4   B   B   B   B   B      B = nó de borda (valor conhecido)
   j=3   B   o   o   o   B      o = nó interno (incógnita)
   j=2   B   o   o   o   B
   j=1   B   o   o   o   B
   j=0   B   B   B   B   B  <- borda inferior (100 °C)
        i=0 i=1 i=2 i=3 i=4
```

#### Dedução da fórmula (série de Taylor)

Expandindo $T$ em série de Taylor em torno de um ponto $x$, para um passo $h$ à direita e à esquerda:

$$T(x+h) = T(x) + h\,T'(x) + \frac{h^2}{2}T''(x) + \frac{h^3}{6}T'''(x) + O(h^4)$$

$$T(x-h) = T(x) - h\,T'(x) + \frac{h^2}{2}T''(x) - \frac{h^3}{6}T'''(x) + O(h^4)$$

Somando as duas expressões, os termos de ordem ímpar se cancelam:

$$T(x+h) + T(x-h) = 2T(x) + h^2\,T''(x) + O(h^4)$$

Isolando a derivada segunda, obtém-se a **diferença central**:

$$\frac{\partial^2 T}{\partial x^2} \approx \frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{h^2}$$

O erro de truncamento dessa aproximação é de ordem $O(h^2)$. De forma análoga, para a direção $y$:

$$\frac{\partial^2 T}{\partial y^2} \approx \frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{h^2}$$

#### A equação em cada nó interno

Substituindo na equação de Laplace, com $\Delta x = \Delta y = h$:

$$\frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{h^2} + \frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{h^2} = 0$$

Multiplicando por $h^2$ e agrupando:

$$4T_{i,j} - T_{i-1,j} - T_{i+1,j} - T_{i,j-1} - T_{i,j+1} = 0$$

Isolando $T_{i,j}$:

$$\boxed{T_{i,j} = \frac{T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1}}{4}}$$

**Interpretação:** a temperatura de cada nó interno é a **média aritmética dos quatro vizinhos** (esquerda, direita, abaixo, acima). Essa é a "propriedade do valor médio" das funções harmônicas, reproduzida no nível discreto. O valor de $h$ desaparece da fórmula final, por isso o programa não precisa dele.

### 4.3 O sistema linear resultante

Escrevendo a equação acima para **todos** os $n^2$ nós internos, obtém-se um sistema linear $A\,\mathbf{T} = \mathbf{b}$ com $n^2$ equações e $n^2$ incógnitas:

- **Matriz $A$:** tem $4$ na diagonal principal e $-1$ nas posições que correspondem aos vizinhos de cada nó. Cada linha tem no máximo 5 elementos não nulos. Ela é **esparsa**, **simétrica**, **definida positiva** e tem estrutura **bloco-tridiagonal**.
- **Vetor $\mathbf{b}$:** é zero, exceto nas equações dos nós vizinhos às bordas, onde o valor conhecido da borda é passado para o lado direito. Por exemplo, para um nó adjacente à borda inferior, o termo $T_{i,0} = 100$ vai para $\mathbf{b}$.

Tamanho do sistema:

| $n$ | Incógnitas ($n^2$) | Tamanho de $A$ (denso) | Elementos não nulos de $A$ ($n^2 + 4n(n-1)$) |
|-----|--------------------|------------------------|-----------------------------------------------|
| 5 | 25 | 25 × 25 = 625 | 105 |
| 10 | 100 | 100 × 100 = 10 000 | 460 |
| 20 | 400 | 400 × 400 = 160 000 | 1 920 |

Montar e guardar $A$ seria desperdício: para $n = 20$, cerca de 98,8 % dos elementos são zero. Por isso o enunciado pede para **não montar a matriz** e usar a fórmula da média diretamente. Métodos iterativos como o Gauss–Seidel permitem isso, porque só precisam saber **quem são os vizinhos** de cada nó.

### 4.4 Métodos iterativos: Jacobi e Gauss–Seidel

Em vez de resolver $A\mathbf{T} = \mathbf{b}$ de uma vez (por exemplo, por eliminação de Gauss), um método iterativo parte de um **chute inicial** e o melhora repetidamente, gerando uma sequência $\mathbf{T}^{(0)}, \mathbf{T}^{(1)}, \mathbf{T}^{(2)}, \dots$ que se aproxima da solução.

Para um sistema geral $A\mathbf{x} = \mathbf{b}$, isola-se a $i$-ésima variável na $i$-ésima equação.

#### Método de Jacobi

$$x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j \ne i} a_{ij}\,x_j^{(k)}\right)$$

Todas as variáveis da iteração $k+1$ são calculadas **somente com valores da iteração $k$**. É preciso guardar duas cópias do vetor (a antiga e a nova).

#### Método de Gauss–Seidel

$$x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j < i} a_{ij}\,x_j^{(k+1)} - \sum_{j > i} a_{ij}\,x_j^{(k)}\right)$$

Aqui, ao calcular $x_i^{(k+1)}$, as variáveis $x_j$ com $j < i$ **já foram atualizadas** nesta mesma varredura, e seus valores novos são usados imediatamente. Só as variáveis com $j > i$ ainda estão na versão antiga.

| | Jacobi | Gauss–Seidel |
|---|--------|--------------|
| Usa valores novos imediatamente? | Não | Sim |
| Memória | Duas cópias do vetor | Uma só (atualização no lugar) |
| Dependência da ordem dos nós | Nenhuma | Sim (a ordem altera os valores intermediários) |
| Velocidade de convergência | Mais lento | Em geral ~2× mais rápido neste problema |
| Paralelização | Trivial | Mais difícil |

### 4.5 Gauss–Seidel aplicado à placa

Neste problema, $a_{ii} = 4$ e os coeficientes dos quatro vizinhos são $-1$. Substituindo na fórmula de Gauss–Seidel, tudo se reduz à média dos vizinhos (os termos de contorno entram naturalmente, porque os nós de borda estão na matriz $T$).

A ordem de varredura adotada no código é:

```
j (vertical)    de 1 até n      <- laço externo
  i (horizontal) de 1 até n     <- laço interno
```

Ou seja, a varredura anda **da esquerda para a direita, de baixo para cima**. Com essa ordem, quando o nó $(i,j)$ é calculado:

| Vizinho | Situação nesta varredura |
|---------|--------------------------|
| $T_{i-1,j}$ (esquerda) | **já atualizado** (novo) |
| $T_{i,j-1}$ (abaixo) | **já atualizado** (novo) |
| $T_{i+1,j}$ (direita) | ainda antigo |
| $T_{i,j+1}$ (acima) | ainda antigo |

```
            acima  (antigo)
               ↑
esquerda → [ T(i,j) ] ← direita
 (novo)        ↑        (antigo)
            abaixo (novo)
```

Isso é exatamente o que o código faz: ele sobrescreve `T[i,j]` imediatamente, e os próximos nós "enxergam" o valor novo.

### 4.6 Convergência

**Quando o método converge?** Para matrizes simétricas definidas positivas, como a deste problema, o método de Gauss–Seidel **sempre converge**, qualquer que seja o chute inicial. A matriz também é diagonalmente dominante (de forma fraca), outra condição clássica favorável.

**Quão rápido?** A velocidade é governada pelo **raio espectral** $\rho$ da matriz de iteração: o erro cai, a cada iteração, por um fator aproximadamente igual a $\rho$. Para a equação de Laplace numa malha quadrada com $h = 1/(n+1)$, os raios espectrais são conhecidos:

$$\rho_{\text{Jacobi}} = \cos(\pi h) \qquad\qquad \rho_{\text{Gauss–Seidel}} = \cos^2(\pi h) = \rho_{\text{Jacobi}}^2$$

Isso explica por que Gauss–Seidel precisa de cerca de **metade** das iterações de Jacobi: cada iteração dele vale por duas de Jacobi.

Quando $h \to 0$ (malha refinada), $\cos^2(\pi h) \approx 1 - \pi^2 h^2$, ou seja, **$\rho \to 1$**: o método fica cada vez mais lento. Para reduzir o erro de um fator $\varepsilon$ são necessárias aproximadamente

$$k \approx \frac{\ln \varepsilon}{\ln \rho_{\text{GS}}} \sim \frac{\ln(1/\varepsilon)}{\pi^2 h^2} \quad\Longrightarrow\quad k = O(n^2)$$

iterações. Comparação com os resultados do programa (com $\varepsilon = 10^{-6}$):

| $n$ | $h = 1/(n+1)$ | $\rho_{\text{GS}} = \cos^2(\pi h)$ | Estimativa $\ln(10^{-6})/\ln\rho$ | Iterações reais |
|-----|---------------|-------------------------------------|-----------------------------------|------------------|
| 5 | 1/6 | 0,7500 | ≈ 48 | 56 |
| 10 | 1/11 | 0,9206 | ≈ 167 | 178 |
| 20 | 1/21 | 0,9778 | ≈ 615 | 598 |

A estimativa não é exata, pois o critério de parada mede a diferença entre iterações (não o erro verdadeiro) e o chute inicial influencia as primeiras iterações, mas a **ordem de grandeza e a tendência** batem muito bem.

**Custo total:** cada varredura custa $O(n^2)$ operações (uma por nó) e são necessárias $O(n^2)$ varreduras, então o custo total cresce como $O(n^4)$. Ao dobrar $n$, o tempo de execução aumenta, assintoticamente, cerca de 16 vezes.

### 4.7 Critério de parada × resíduo

O programa usa **dois** indicadores diferentes. Vale entender a diferença:

**1. Erro de parada (variação entre varreduras)**

$$\text{erro} = \max_{i,j}\ \lvert T_{i,j}^{(k+1)} - T_{i,j}^{(k)} \rvert$$

É o que decide quando o laço `while` termina. Mede **quanto a solução ainda está se movendo**. É barato (calculado durante a própria varredura), mas não é uma medida direta de quão correta é a solução.

**2. Resíduo**

$$r_{i,j} = 4T_{i,j} - \left(T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1}\right)$$

Mede **o quanto cada equação do sistema linear deixa de ser satisfeita** pela solução obtida. Se a solução fosse exata, todos os $r_{i,j}$ seriam zero. É calculado uma única vez, no final, e funciona como verificação independente.

O resíduo é naturalmente da ordem de **quatro vezes** a variação de um nó (pois, em módulo, $r = 4\,\Delta T$ no momento da atualização). Por isso, com `tol = 1e-6`, o resíduo final fica na casa de $10^{-6}$ (aqui, entre $1{,}7\times10^{-6}$ e $1{,}9\times10^{-6}$).

> **Cuidado:** um erro de parada pequeno nem sempre garante erro pequeno na solução verdadeira, principalmente quando $\rho$ é próximo de 1 (malhas finas, convergência lenta). Por isso a conferência pelo resíduo é uma boa prática.

---

## 5. Explicação detalhada do código

O programa tem duas partes: a função `resolver_gauss_seidel`, que faz todo o cálculo, e o bloco `if __name__ == "__main__":`, que a executa para várias malhas e apresenta os resultados.

### 5.1 Importações

```python
import numpy as np
import matplotlib.pyplot as plt
```

- `numpy` fornece o array 2D (`np.full`) que guarda as temperaturas.
- `matplotlib.pyplot` gera o gráfico no final.

### 5.2 Assinatura da função

```python
def resolver_gauss_seidel(n, tol=1e-6, max_iter=100000):
```

| Parâmetro | Significado | Valor padrão |
|-----------|-------------|--------------|
| `n` | Número de nós internos em cada direção | (obrigatório) |
| `tol` | Tolerância do critério de parada, em °C | $10^{-6}$ |
| `max_iter` | Número máximo de iterações (segurança contra laço infinito) | 100000 |

A função **retorna** quatro valores: a matriz `T` com as temperaturas, o número de iterações, o erro final e o resíduo máximo.

### 5.3 Passo 1: malha e condições de contorno

```python
T = np.full((n + 2, n + 2), 20.0)
T[1:n+1, 0] = 100.0
```

**Linha 1:** cria uma matriz $(n+2)\times(n+2)$ preenchida com 20,0. O "+2" existe para incluir as duas bordas de cada direção. Isso resolve dois problemas de uma vez:

- **Chute inicial:** todos os nós internos começam em 20 °C.
- **Bordas laterais e superior:** já ficam com o valor correto (20 °C), sem precisar de nenhum comando adicional.

**Linha 2:** fixa a borda inferior em 100 °C. O fatiamento `T[1:n+1, 0]` seleciona:

- primeiro índice `1:n+1`, ou seja, $i = 1, 2, \dots, n$ (o Python exclui o último, então `n+1` não entra);
- segundo índice `0`, ou seja, $j = 0$, a linha de baixo.

**Por que os cantos ficam fora?** `T[0,0]` e `T[n+1,0]` continuam em 20 °C porque o fatiamento começa em 1 e termina em $n$. Os cantos nunca aparecem como vizinho de nenhum nó interno (os vizinhos são só esquerda, direita, abaixo e acima), então o valor deles não influencia o resultado.

#### Tabela de correspondência de índices (Fortran × Python)

O enunciado usa índices de 1 a $n+2$ (convenção do Fortran/MATLAB). O Python usa de 0 a $n+1$. Por isso, **todo índice do Python é o do enunciado menos 1**:

| Conceito | Enunciado (1-based) | Python (0-based) |
|----------|--------------------|-------------------|
| Primeira/última linha ou coluna | 1 e $n+2$ | 0 e $n+1$ |
| Nós internos | 2 até $n+1$ | 1 até $n$ → `range(1, n+1)` |
| Borda inferior | $j = 1$ | $j = 0$ |
| Nó central para $n=5$ | `T(4,4)` | `T[3,3]` |

#### Convenção dos eixos

`T[i, j]`: **`i` é a posição horizontal ($x$) e `j` é a posição vertical ($y$)**. Isso é diferente da convenção usual de matrizes (linha, coluna), onde o primeiro índice seria a linha (vertical). Essa escolha segue o enunciado, mas exige atenção na hora de desenhar (ver 5.7).

### 5.4 Passo 2: laço iterativo

```python
iteracao = 0
erro = float('inf')

while erro > tol and iteracao < max_iter:
    erro = 0.0
    for j in range(1, n + 1):
        for i in range(1, n + 1):
            antigo = T[i, j]
            T[i, j] = (T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1]) / 4.0
            diferenca = abs(T[i, j] - antigo)
            if diferenca > erro:
                erro = diferenca
    iteracao += 1
```

**Inicialização.** `erro = float('inf')` garante que o `while` entre na primeira vez (infinito é maior que qualquer tolerância). `iteracao` conta as varreduras.

**Condição do `while`.** O laço continua enquanto o erro for maior que a tolerância **e** o número máximo de iterações não tiver sido atingido. Ele para por convergência ou por segurança, o que ocorrer primeiro.

**`erro = 0.0` no início de cada varredura.** Zera o "maior erro" para calcular o máximo apenas desta varredura.

**Laços `for`.** `j` é o laço externo (vertical) e `i` o interno (horizontal), como visto em 4.5. Ambos usam `range(1, n+1)`, percorrendo só os nós internos, **sem tocar nas bordas**.

**Corpo do laço:**

1. `antigo = T[i, j]`: guarda o valor antes de modificá-lo, para medir a variação.
2. `T[i, j] = (...) / 4.0`: aplica a fórmula da média. Essa é a linha que **implementa o método de Gauss–Seidel**. Como a escrita é feita na própria matriz `T`, `T[i-1, j]` e `T[i, j-1]` já têm o valor novo quando são lidos.
3. `diferenca = abs(T[i, j] - antigo)`: quanto esse nó mudou.
4. `if diferenca > erro: erro = diferenca`: mantém o **máximo** das diferenças (equivale a `erro = max(erro, diferenca)`).

**`iteracao += 1`** ao fim de cada varredura completa.

> **Por que não usar operações vetorizadas do NumPy?** Algo como `T[1:-1,1:-1] = (T[:-2,1:-1] + T[2:,1:-1] + ...) / 4` calcularia todos os nós **a partir dos valores antigos**, o que é o método de **Jacobi**, não o de Gauss–Seidel. Os laços explícitos são necessários para que cada nó use imediatamente os valores recém-atualizados.

### 5.5 Conferência: resíduo máximo

```python
max_res = 0.0
for j in range(1, n + 1):
    for i in range(1, n + 1):
        r = 4.0 * T[i, j] - (T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1])
        if abs(r) > max_res:
            max_res = abs(r)

return T, iteracao, erro, max_res
```

Depois que o laço iterativo termina, percorre-se novamente os nós internos calculando o resíduo $r_{i,j}$ (ver 4.7) e guardando o maior valor absoluto. Isso valida se a solução obtida de fato satisfaz as equações discretizadas. Por fim, a função devolve os quatro resultados.

### 5.6 As questões teóricas (comentários no código)

As respostas às questões rápidas do enunciado estão escritas como comentários entre a função e o bloco principal. Versões completas na [seção 8](#8-questões-teóricas-respondidas).

### 5.7 Bloco principal e apresentação dos resultados

```python
if __name__ == "__main__":
```

Garante que o bloco só rode quando o arquivo é executado diretamente (`python gauss_seidel.py`), e **não** quando é importado como módulo por outro programa. Assim a função `resolver_gauss_seidel` pode ser reutilizada em outros scripts sem disparar as simulações.

**Execução para $n = 5$:**

```python
T_5, iter_5, erro_5, res_5 = resolver_gauss_seidel(5)
centro_idx_5 = 3
print(f"... {T_5[centro_idx_5, centro_idx_5]:.4f} ºC")
```

A função devolve uma tupla de quatro elementos, que é desempacotada em quatro variáveis. Para $n = 5$ a matriz tem $7 \times 7$ posições (índices 0 a 6) e o nó central é o de índice $3$, que corresponde a `T(4,4)` do enunciado. O formato `:.4f` imprime 4 casas decimais, e `:.4e` imprime em notação científica (usado para erros, que são números muito pequenos).

**Execuções para $n = 10$ e $n = 20$:** repetem o mesmo processo e imprimem iterações, erro e resíduo. Como o número de nós internos é par, o centro geométrico da placa **cai entre dois nós**, e não existe um único nó central. Por isso o código não imprime a temperatura central para essas malhas.

**Gráfico (item 4a do enunciado):**

```python
plt.figure(figsize=(7, 6))
plt.imshow(T_20.T, origin='lower', cmap='jet', extent=[0, 1, 0, 1])
plt.colorbar(label='Temperatura (°C)')
plt.title('Distribuição de Temperatura na Placa (Gauss-Seidel, n=20)')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('grafico_calor_placa.png')
```

| Comando | O que faz |
|---------|-----------|
| `plt.figure(figsize=(7, 6))` | Cria a figura com 7 × 6 polegadas |
| `T_20.T` | Transpõe a matriz. O `imshow` espera `[linha, coluna] = [y, x]`, mas aqui `T[i, j] = T[x, y]`, então é preciso trocar |
| `origin='lower'` | Põe a linha 0 (borda a 100 °C) **embaixo** do desenho. O padrão do `imshow` seria em cima |
| `cmap='jet'` | Mapa de cores azul (frio) → vermelho (quente), igual ao da Figura 1 do enunciado |
| `extent=[0, 1, 0, 1]` | Rotula os eixos de 0 a 1 m (lado da placa) |
| `plt.colorbar(...)` | Mostra a escala de cores com o rótulo da temperatura |
| `plt.savefig(...)` | Salva a imagem em arquivo PNG |

Sem a **transposição** e o **`origin='lower'`**, o gráfico sairia girado ou de cabeça para baixo, com a borda quente no lado errado.

### 5.8 Glossário das variáveis do código

| Variável | Tipo | Significado |
|----------|------|-------------|
| `n` | `int` | Número de nós internos por direção |
| `tol` | `float` | Tolerância do critério de parada (°C) |
| `max_iter` | `int` | Limite de iterações |
| `T` | `ndarray` $(n+2)\times(n+2)$ | Temperaturas de todos os pontos, incluindo bordas |
| `iteracao` | `int` | Contador de varreduras |
| `erro` | `float` | Maior variação de temperatura na última varredura |
| `antigo` | `float` | Valor de `T[i,j]` antes da atualização |
| `diferenca` | `float` | $\lvert T_{novo} - T_{antigo}\rvert$ de um nó |
| `max_res` | `float` | Maior valor absoluto do resíduo $r_{i,j}$ |
| `r` | `float` | Resíduo de um nó |

---

## 6. Resultados

### 6.1 Resumo das simulações

| $n$ | Nós internos | Iterações | Erro final (máx. ΔT) | Resíduo máximo | Convergiu? |
|-----|--------------|-----------|----------------------|----------------|------------|
| 5 | 25 | 56 | $8{,}49\times10^{-7}$ | $1{,}70\times10^{-6}$ | Sim |
| 10 | 100 | 178 | $9{,}70\times10^{-7}$ | $1{,}86\times10^{-6}$ | Sim |
| 20 | 400 | 598 | $9{,}79\times10^{-7}$ | $1{,}94\times10^{-6}$ | Sim |

Em todos os casos o erro final ficou abaixo de $10^{-6}$ antes de atingir `max_iter = 100000`, ou seja, houve convergência. O programa não imprime a palavra "convergiu", mas isso se verifica comparando o erro final com `tol` (o laço só termina com `erro < tol` quando a convergência foi atingida).

### 6.2 Temperatura no nó central ($n = 5$)

$$T(4,4) = \texttt{T[3,3]} \approx 40{,}0000\ °\text{C}$$

### 6.3 Matriz de temperaturas para $n = 5$

Topo = borda superior (20 °C); base = borda aquecida (100 °C). Colunas da esquerda para a direita ($i = 0 \dots 6$):

```
20.0000  20.0000  20.0000  20.0000  20.0000  20.0000  20.0000    <- j = 6 (borda superior)
20.0000  22.5050  24.2774  24.9075  24.2774  22.5051  20.0000    <- j = 5
20.0000  25.7428  29.6970  31.0754  29.6970  25.7428  20.0000    <- j = 4
20.0000  30.7692  37.6923  40.0000  37.6923  30.7692  20.0000    <- j = 3
20.0000  39.6418  50.3030  53.5400  50.3030  39.6418  20.0000    <- j = 2
20.0000  57.4949  70.3380  73.5540  70.3380  57.4949  20.0000    <- j = 1
20.0000 100.0000 100.0000 100.0000 100.0000 100.0000  20.0000    <- j = 0 (borda inferior)
```

### 6.4 Gráfico gerado ($n = 20$)

![Distribuição de temperatura na placa, n = 20](grafico_calor_placa.png)

A borda inferior é a região mais quente (vermelho) e a temperatura decai suavemente para cima e em direção às paredes laterais (azul). Perto dos cantos inferiores o gradiente é muito forte, pois ali a borda a 100 °C encontra as paredes a 20 °C.

---

## 7. Análise dos resultados

### 7.1 Por que $T(4,4) = 40$ °C exatamente?

Há uma justificativa analítica elegante. Subtraindo 20 °C de todas as temperaturas, o problema fica com uma borda a 80 °C e as outras três a 0 °C (a equação de Laplace continua valendo, pois ela é linear e uma constante a satisfaz).

Agora imagine os **quatro** problemas em que cada um dos lados, separadamente, é o único a 80 °C. Somando as quatro soluções, obtém-se o problema com **todas** as bordas a 80 °C, cuja solução é simplesmente 80 °C em toda a placa (pelo princípio do máximo). No **centro**, pela simetria do quadrado (girar a placa de 90° troca um problema pelo outro), as quatro soluções valem o mesmo. Logo, cada uma vale $80/4 = 20$ °C no centro. Voltando à escala original:

$$T_{centro} = 20 + 20 = 40\ °\text{C} \qquad\Big(\text{equivalente a } \tfrac{100 + 20 + 20 + 20}{4}\Big)$$

O resultado numérico (39,999997 °C, impresso como 40,0000) confirma que o programa está correto. Isso vale **somente** para $n$ ímpar, quando existe um nó exatamente no centro da placa.

### 7.2 Outras verificações de sanidade

| Verificação | Resultado |
|-------------|-----------|
| **Princípio do máximo** (todo valor interno entre 20 e 100 °C) | Satisfeito |
| **Simetria** em relação ao eixo vertical ($T_{i,j} = T_{n+1-i,j}$) | Satisfeita (ver matriz 6.3) |
| **Monotonicidade vertical** (a temperatura diminui de baixo para cima na coluna central) | Satisfeita: 100 → 73,55 → 53,54 → 40,00 → 31,08 → 24,91 → 20 |
| **Resíduo** pequeno | $\sim10^{-6}$ em todas as malhas |

Pequenas diferenças na quarta casa decimal entre pontos simétricos (por exemplo, 22,5050 e 22,5051) vêm de o método parar antes da convergência total, combinado com a ordem de varredura, que não é simétrica.

### 7.3 Efeito do refinamento da malha

| Mudança | Razão entre iterações |
|---------|-----------------------|
| $n = 5 \to 10$ | $178/56 \approx 3{,}18$ |
| $n = 10 \to 20$ | $598/178 \approx 3{,}36$ |

Dobrar $n$ multiplica o número de iterações por um valor que **tende a 4**, em acordo com a previsão teórica $k = O(n^2)$ da seção 4.6.

### 7.4 Sobre o centro nas malhas $n = 10$ e $n = 20$

Nessas malhas o centro geométrico fica entre quatro nós. Dois nós diagonalmente opostos em torno do centro têm, no resultado do programa:

| $n$ | Nó (índice Python) | Temperatura |
|-----|--------------------|-------------|
| 10 | `T[5,5]` | 42,9999 °C |
| 10 | `T[6,6]` | 37,0000 °C |
| 20 | `T[10,10]` | 41,5845 °C |
| 20 | `T[11,11]` | 38,4155 °C |

Um deles fica acima e o outro abaixo de 40 °C (o de baixo, mais perto da borda quente, é o mais quente), e a média de cada par é 40 °C nos resultados numéricos. O afastamento em relação a 40 °C diminui com o refinamento (cerca de 3 °C para $n = 10$ e 1,6 °C para $n = 20$), aproximando-se do valor exato no centro da placa.

---

## 8. Questões teóricas respondidas

### a) Por que as temperaturas das bordas não são atualizadas?

Porque as bordas são **condições de contorno de Dirichlet**: o problema físico impõe o valor da temperatura nelas (100 °C embaixo, 20 °C nas demais). Elas são **dados do problema**, não incógnitas. Só os nós internos têm uma equação discretizada (a média dos vizinhos) para ser resolvida. Os valores das bordas entram apenas como termos conhecidos nas médias dos nós adjacentes. Se fossem atualizadas, o código deixaria de representar uma parede mantida a temperatura constante e a solução perderia o sentido físico.

No código, isso é garantido pelos laços `range(1, n+1)`, que nunca visitam os índices $0$ e $n+1$.

### b) O que diferencia a atualização de Gauss–Seidel da de Jacobi?

No **Jacobi**, toda a nova varredura é calculada usando **apenas valores da iteração anterior**, o que exige manter uma cópia da matriz. No **Gauss–Seidel**, cada valor recém-calculado é usado **imediatamente** nos cálculos seguintes da mesma varredura (neste código, `T[i-1,j]` e `T[i,j-1]` já são novos quando `T[i,j]` é calculado). Isso:

- dispensa a segunda matriz (menos memória);
- acelera a convergência: o raio espectral do Gauss–Seidel é o **quadrado** do de Jacobi, então precisa de cerca de **metade das iterações**.

### c) Como o número de iterações muda quando a malha é refinada?

Ele **aumenta bastante**, aproximadamente **4 vezes** a cada vez que $n$ dobra (56 → 178 → 598 neste trabalho). Há duas razões complementares:

1. **Propagação da informação:** a informação da borda só se espalha de vizinho para vizinho, um nó por varredura. Com mais nós entre a borda e o centro, são necessárias mais varreduras para ela chegar lá.
2. **Raio espectral:** matematicamente, $\rho_{GS} = \cos^2(\pi h) \approx 1 - \pi^2 h^2 \to 1$ quando $h \to 0$. Os erros "suaves" (de baixa frequência) são amortecidos muito lentamente.

Além disso, cada iteração também custa mais (mais nós por varredura), então o tempo total cresce ainda mais rápido, como $O(n^4)$.

---

## 9. Possíveis melhorias e adaptações

### Adaptar o problema

| Para mudar... | O que alterar |
|---------------|---------------|
| Tolerância | Parâmetro `tol` na chamada, por exemplo `resolver_gauss_seidel(20, tol=1e-8)` |
| Malha | Parâmetro `n` |
| Temperatura da borda quente | Valor em `T[1:n+1, 0] = 100.0` |
| Temperatura das bordas frias / chute inicial | Valor em `np.full(..., 20.0)` |
| Outra borda quente (por exemplo, a superior) | Trocar a linha de contorno por `T[1:n+1, n+1] = 100.0` |

### Melhorias possíveis no programa

- **Exibir o gráfico na tela:** adicionar `plt.show()` ao final (hoje o gráfico só é salvo em arquivo).
- **Imprimir se houve convergência:** por exemplo, `print("Convergiu" if erro < tol else "Não convergiu")`.
- **Salvar a matriz em arquivo de dados:** `np.savetxt("temperaturas_n5.dat", T_5.T[::-1], fmt="%9.4f")` (a transposição e a inversão colocam a borda quente embaixo, como no desenho da placa).
- **Gráficos para as três malhas**, e não só para $n = 20$, em uma figura com três subgráficos.
- **Comparar com Jacobi:** implementar a versão vetorizada com NumPy e verificar que ela usa aproximadamente o dobro de iterações.
- **Acelerar o Gauss–Seidel:** usar sobrerrelaxação (**SOR**) com $\omega$ ótimo, que reduz o número de iterações de $O(n^2)$ para $O(n)$; ou compilar os laços com Numba.
- **Rótulos geométricos mais precisos no gráfico:** o `imshow` desenha cada ponto como um "pixel" quadrado. Para uma representação rigorosa das coordenadas ($x_i = i\,h$), pode-se usar `plt.contourf`.

---

## 10. Glossário

| Termo | Significado |
|-------|-------------|
| **Regime permanente** | Estado em que as temperaturas não variam mais com o tempo |
| **Equação de Laplace** | $\nabla^2 T = 0$; descreve o equilíbrio térmico sem fontes |
| **Diferenças finitas** | Técnica que troca derivadas por diferenças entre valores de pontos vizinhos |
| **Nó** | Ponto da malha onde se calcula a temperatura |
| **Nó interno / de borda** | Nó incógnita / nó de valor conhecido |
| **Condição de contorno de Dirichlet** | Condição que fixa o *valor* da função na borda |
| **Sistema linear esparso** | Sistema $A\mathbf{x} = \mathbf{b}$ cuja matriz tem quase todos os elementos nulos |
| **Método iterativo** | Método que gera uma sequência de aproximações que converge para a solução |
| **Jacobi** | Método iterativo que só usa valores da iteração anterior |
| **Gauss–Seidel** | Método iterativo que usa os valores novos assim que disponíveis |
| **Varredura** | Uma passada completa por todos os nós internos |
| **Raio espectral** ($\rho$) | Maior autovalor em módulo da matriz de iteração; controla a velocidade de convergência |
| **Resíduo** | $r = b - A\mathbf{x}$; mede quanto o sistema deixa de ser satisfeito |
| **Convergência** | Aproximação da sequência iterativa à solução do sistema |
| **Vetorização** | Operar sobre arrays inteiros de uma vez, sem laços explícitos |
| **Princípio do máximo** | Os valores internos de uma solução de Laplace ficam entre o mínimo e o máximo das bordas |

---

## Autor

<!-- Preencha com seus dados -->
**Nome:** _seu nome aqui_
**Disciplina:** Cálculo Numérico Aplicado
**Instituição:** _sua instituição aqui_
