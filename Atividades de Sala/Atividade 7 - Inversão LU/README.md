# Inversão Matricial por Decomposição LU
### Aplicação à Flexibilidade de uma Estrutura de Molas (5 graus de liberdade)

Atividade prática de **Cálculo Numérico Aplicado**. O programa calcula a **matriz de flexibilidade** `C = K⁻¹` de uma cadeia de cinco nós ligados por seis molas, usando **decomposição LU implementada manualmente**, sem funções prontas de fatoração ou inversão (`np.linalg.inv`, `np.linalg.solve`, `scipy.linalg.lu` etc.).

Este documento serve como material de consulta e está dividido em duas grandes partes: a **teoria** (física e matemática) e a **programação** (explicação do código bloco a bloco).

---

## Sumário

1. [Como executar](#1-como-executar)
2. [Visão geral do problema](#2-visão-geral-do-problema)
3. [Teoria física](#3-teoria-física)
4. [Teoria matemática: decomposição LU](#4-teoria-matemática-decomposição-lu)
5. [Teoria matemática: inversão por LU](#5-teoria-matemática-inversão-por-lu)
6. [Explicação do código, bloco a bloco](#6-explicação-do-código-bloco-a-bloco)
7. [Exemplo numérico acompanhado passo a passo](#7-exemplo-numérico-acompanhado-passo-a-passo)
8. [Saída do programa e resultados](#8-saída-do-programa-e-resultados)
9. [Interpretação física dos resultados](#9-interpretação-física-dos-resultados)
10. [Questões de discussão respondidas](#10-questões-de-discussão-respondidas)
11. [Limitações e possíveis extensões](#11-limitações-e-possíveis-extensões)
12. [Glossário](#12-glossário)

---

## 1. Como executar

```bash
pip install numpy
python inversao_lu.py
```

**Requisitos:** Python 3 e NumPy.

> **Nota sobre o arquivo:** o código foi escrito em um arquivo sem extensão (`Inversão_LU`). Para o GitHub e para executar normalmente, renomeie para `inversao_lu.py` (sem acento e com `.py`).
>
> Os comentários do código contêm marcações `[cite: 3]`, resíduos de cópia do enunciado. Elas não afetam a execução, mas convém apagá-las para deixar o código limpo.

---

## 2. Visão geral do problema

Uma estrutura mecânica linear é formada por **5 nós móveis** (massas) ligados em série por **6 molas**. As molas das extremidades estão presas a paredes fixas.

```
|parede|--k0--(1)--k1--(2)--k2--(3)--k3--(4)--k4--(5)--k5--|parede|
              u1       u2       u3       u4       u5
```

- `u1 … u5`: deslocamentos horizontais dos nós.
- `F1 … F5`: forças externas aplicadas em cada nó.
- `k0 … k5`: rigidezes das molas, em N/m.

**Objetivo:** em vez de resolver o sistema para um único carregamento, determinar a **matriz de flexibilidade** `C`, que permite calcular os deslocamentos para **qualquer** carregamento por meio de `u = C·F`.

**Dados usados:**

| Mola | k0 | k1 | k2 | k3 | k4 | k5 |
|---|---|---|---|---|---|---|
| Rigidez (N/m) | 1200 | 800 | 1500 | 1000 | 1800 | 900 |

---

## 3. Teoria física

### 3.1 Lei de Hooke

Uma mola de rigidez `k` que sofre alongamento `Δ` exerce uma força restauradora de módulo `k·Δ`. A relação é linear, o que torna o problema linear.

### 3.2 Equilíbrio de um nó

No nó `i`, duas molas atuam: a da esquerda (`k_{i-1}`) e a da direita (`k_i`).

- O alongamento da mola da esquerda é `uᵢ − uᵢ₋₁`.
- O alongamento da mola da direita é `uᵢ₊₁ − uᵢ`.

No equilíbrio estático, a força externa equilibra as forças das molas:

$$F_i = k_{i-1}(u_i - u_{i-1}) - k_i(u_{i+1} - u_i)$$

Reorganizando:

$$F_i = -k_{i-1}\,u_{i-1} + (k_{i-1}+k_i)\,u_i - k_i\,u_{i+1}$$

Cada equação envolve apenas o nó `i` e seus vizinhos. Isso explica a estrutura **tridiagonal** da matriz.

### 3.3 Matriz de rigidez K

Escrevendo as 5 equações juntas, obtém-se `F = K·u`:

$$K=\begin{bmatrix}k_0+k_1&-k_1&0&0&0\\-k_1&k_1+k_2&-k_2&0&0\\0&-k_2&k_2+k_3&-k_3&0\\0&0&-k_3&k_3+k_4&-k_4\\0&0&0&-k_4&k_4+k_5\end{bmatrix}$$

Com os valores numéricos:

$$K=\begin{bmatrix}2000&-800&0&0&0\\-800&2300&-1500&0&0\\0&-1500&2500&-1000&0\\0&0&-1000&2800&-1800\\0&0&0&-1800&2700\end{bmatrix}\ \text{N/m}$$

**Como ler a matriz:**

- **Diagonal `K[i,i]`:** soma das rigidezes das duas molas ligadas ao nó `i`. É a força necessária para deslocar só o nó `i` em 1 m mantendo os outros parados.
- **Fora da diagonal `K[i,i+1]`:** `−k_i`, o acoplamento causado pela mola entre os nós `i` e `i+1`.
- **Zeros:** nós não vizinhos não interagem diretamente.

### 3.4 Propriedades importantes de K

1. **Simétrica:** `K = Kᵀ`. A mola entre `i` e `j` exerce em `i` a mesma intensidade de força que exerce em `j`.
2. **Tridiagonal:** só há elementos na diagonal principal e nas duas vizinhas.
3. **Definida positiva:** `uᵀ K u > 0` para todo `u ≠ 0`. De fato:

$$u^T K u = k_0u_1^2 + \sum_{i=1}^{4}k_i(u_{i+1}-u_i)^2 + k_5u_5^2$$

Essa soma é o **dobro da energia elástica** armazenada nas molas. Como todo `kᵢ > 0`, ela só é zero se todos os alongamentos forem nulos, ou seja, se `u = 0`. Isso depende das paredes (ver a [questão g](#g-se-todas-as-restrições-às-extremidades-fossem-removidas-seria-possível-inverter-a-matriz-de-rigidez)).

4. **Diagonalmente dominante** (com dominância estrita nas linhas das pontas).

Consequência prática de ser simétrica definida positiva: a decomposição LU **existe e é estável sem pivoteamento**, pois todos os pivôs são positivos.

### 3.5 Matriz de flexibilidade C

Se `K` é inversível:

$$u = K^{-1}F \quad\Rightarrow\quad C = K^{-1},\qquad u = C\,F$$

| Matriz | Relaciona | Unidade | Ideia |
|---|---|---|---|
| `K` (rigidez) | deslocamentos → forças | N/m | "quão dura" é a estrutura |
| `C` (flexibilidade / complacência) | forças → deslocamentos | m/N | "quão macia" é a estrutura |

### 3.6 Significado das colunas e dos elementos de C

Multiplicar `C` por um vetor da base canônica `eⱼ` seleciona a coluna `j`:

$$C\,e_j = c_j$$

Logo, **a coluna `j` de `C` é o vetor de deslocamentos de todos os nós quando se aplica uma força de 1 N apenas no nó `j`.**

Portanto, `Cᵢⱼ` é o deslocamento do nó `i` causado por uma força unitária no nó `j`.

Como o sistema é linear, vale a **superposição**: para um carregamento qualquer,

$$u = C\,F = \sum_j F_j\,c_j$$

ou seja, o deslocamento total é a soma das respostas a cada força, ponderadas pelo seu valor.

### 3.7 Simetria de C (reciprocidade de Maxwell–Betti)

Como `K = Kᵀ`, a inversa também é simétrica:

$$C^T = (K^{-1})^T = (K^T)^{-1} = K^{-1} = C$$

Fisicamente, `Cᵢⱼ = Cⱼᵢ` significa que **o deslocamento do nó `i` devido a uma força em `j` é igual ao deslocamento do nó `j` devido a uma força em `i`**. Esse é o teorema da reciprocidade, válido para estruturas lineares e elásticas.

---

## 4. Teoria matemática: decomposição LU

### 4.1 Ideia

Fatorar `K` como o produto de duas matrizes triangulares:

$$K = L\,U$$

- `L` (*lower*): triangular **inferior**, com **diagonal unitária** (convenção usada neste código).
- `U` (*upper*): triangular **superior**.

Sistemas triangulares são muito baratos de resolver, por substituição. Fatorar uma vez e reaproveitar é a vantagem do método.

### 4.2 De onde vem a decomposição: eliminação de Gauss

A eliminação de Gauss transforma `K` em uma matriz triangular superior subtraindo múltiplos de uma linha das linhas abaixo dela.

No passo `k`, para cada linha `i > k`:

$$\ell_{ik} = \frac{U_{ik}}{U_{kk}}, \qquad \text{linha}_i \leftarrow \text{linha}_i - \ell_{ik}\cdot \text{linha}_k$$

- `U_kk` é o **pivô**.
- `ℓᵢₖ` é o **multiplicador** (ou fator).
- Essa operação zera o elemento `U[i,k]`.

Cada passo equivale a multiplicar `K` à esquerda por uma matriz triangular inferior `M_k`:

$$U = M_{n-1}\cdots M_2 M_1\,K \quad\Rightarrow\quad K = M_1^{-1}M_2^{-1}\cdots M_{n-1}^{-1}\,U = L\,U$$

O resultado notável é que, ao multiplicar as inversas, os multiplicadores **não se misturam**: `L` é simplesmente a matriz com 1 na diagonal e os multiplicadores `ℓᵢₖ` nas posições abaixo da diagonal. Por isso o código só precisa guardar `fator` em `L[i,k]`.

### 4.3 Condição de existência

A LU **sem pivoteamento** funciona quando todos os pivôs `U_kk` são diferentes de zero. Para matrizes **simétricas definidas positivas** (caso de `K`), isso é garantido, e os pivôs são todos positivos. Por isso o código não faz troca de linhas.

O **determinante** vem de graça:

$$\det K = \prod_{k=1}^{n} U_{kk}$$

(pois `det L = 1`). Se algum pivô for zero, `K` é singular.

### 4.4 Custo

Para uma matriz densa `n × n`, a fatoração custa cerca de `(2/3)n³` operações de ponto flutuante.

---

## 5. Teoria matemática: inversão por LU

### 5.1 Inversa como `n` sistemas lineares

Por definição, `K·C = I`. Olhando coluna a coluna, e chamando de `cⱼ` a coluna `j` de `C` e `eⱼ` a coluna `j` de `I`:

$$K\,c_j = e_j, \qquad j = 1, 2, \dots, n$$

Ou seja, calcular a inversa equivale a resolver `n` sistemas lineares com a **mesma matriz** `K` e lados direitos diferentes.

### 5.2 Uso da fatoração

Substituindo `K = LU`:

$$L\,U\,c_j = e_j$$

Define-se uma variável auxiliar `yⱼ = U·cⱼ` e resolvem-se dois sistemas triangulares, em sequência:

1. **Substituição progressiva** (`L·yⱼ = eⱼ`): de cima para baixo.
2. **Substituição regressiva** (`U·cⱼ = yⱼ`): de baixo para cima.

### 5.3 Substituição progressiva

Como `L` é triangular inferior, a primeira equação envolve só `y₁`, a segunda `y₁` e `y₂`, e assim por diante:

$$y_i = \frac{e_i - \sum_{j=1}^{i-1} L_{ij}\,y_j}{L_{ii}}, \qquad i = 1,\dots,n$$

(Aqui `L_ii = 1`, então a divisão não altera nada, mas ela é mantida por generalidade.)

### 5.4 Substituição regressiva

Como `U` é triangular superior, a última equação envolve só `xₙ`:

$$x_n = \frac{y_n}{U_{nn}}, \qquad x_i = \frac{y_i - \sum_{j=i+1}^{n} U_{ij}\,x_j}{U_{ii}}, \quad i = n-1,\dots,1$$

O vetor `x` resultante é a coluna `j` de `C`.

### 5.5 Por que a LU é feita uma única vez

A fatoração depende **somente de `K`**. O que muda de uma coluna para outra é apenas o lado direito `eⱼ`. Logo, `L` e `U` são reutilizadas, e só as substituições são repetidas.

| Etapa | Custo |
|---|---|
| Fatoração LU (uma vez) | ≈ `(2/3)n³` |
| Um par de substituições (por coluna) | ≈ `2n²` |
| `n` colunas | ≈ `2n³` |
| **Total** | **da ordem de `n³`** |

Se a eliminação fosse refeita para cada coluna, o custo total seria da ordem de `n⁴`.

### 5.6 Verificações numéricas

- **`LU ≈ K`:** confirma que a fatoração está correta.
- **`KC ≈ I`:** confirma que `C` é a inversa. Mede-se `erromax = max|KC − I|`.
- **`C ≈ Cᵀ`:** confirma a simetria. Mede-se `errosim = max|Cᵢⱼ − Cⱼᵢ|`.

Os erros não são exatamente zero devido ao **arredondamento de ponto flutuante** (precisão da ordem de `10⁻¹⁶` em dupla precisão). Erros dessa magnitude indicam um resultado correto.

---

## 6. Explicação do código, bloco a bloco

### 6.1 Importação e dados

```python
import numpy as np

n = 5

k0 = 1200.0
k1 = 800.0
k2 = 1500.0
k3 = 1000.0
k4 = 1800.0
k5 = 900.0
```

- `numpy` é usado para armazenar e imprimir matrizes. As operações do método (LU, substituições) são feitas por laços explícitos.
- `n` é o número de graus de liberdade.
- Os valores são escritos como `float` (`1200.0`). Isso evita divisão inteira e deixa explícito o tipo numérico.

### 6.2 Construção de K

```python
K = np.array([
    [k0 + k1, -k1,      0.0,     0.0,     0.0],
    [-k1,     k1 + k2, -k2,      0.0,     0.0],
    [0.0,    -k2,       k2 + k3, -k3,     0.0],
    [0.0,     0.0,     -k3,      k3 + k4, -k4],
    [0.0,     0.0,      0.0,    -k4,      k4 + k5]
])
```

Esta é a matriz da seção 3.3. Cada elemento da diagonal é a soma das duas molas vizinhas, e cada elemento ao lado da diagonal é o negativo da mola compartilhada. O programa imprime `K` para conferência antes de seguir.

### 6.3 Decomposição LU

```python
L = np.eye(n)
U = K.copy()

for k in range(n - 1):
    for i in range(k + 1, n):
        fator = U[i, k] / U[k, k]
        L[i, k] = fator
        for j in range(k + 1, n):
            U[i, j] = U[i, j] - fator * U[k, j]
        U[i, k] = 0.0
```

| Linha | O que faz |
|---|---|
| `L = np.eye(n)` | Começa `L` como identidade (diagonal 1, resto 0). Os multiplicadores serão preenchidos abaixo da diagonal. |
| `U = K.copy()` | Começa `U` como cópia de `K`. O `.copy()` é essencial: sem ele, `U` e `K` seriam o **mesmo objeto** e modificar `U` destruiria `K`, que é necessário depois para as verificações. |
| `for k in range(n-1)` | Percorre os pivôs. O último pivô não precisa eliminar nada, por isso o laço vai só até o índice `n−2`. |
| `for i in range(k+1, n)` | Percorre as linhas abaixo do pivô, que precisam ser zeradas na coluna `k`. |
| `fator = U[i,k]/U[k,k]` | Calcula o multiplicador `ℓᵢₖ`. |
| `L[i,k] = fator` | Guarda o multiplicador em `L`. |
| `for j in range(k+1, n)` | Atualiza as colunas à direita do pivô. |
| `U[i,j] = U[i,j] - fator*U[k,j]` | Subtrai da linha `i` o múltiplo da linha do pivô. |
| `U[i,k] = 0.0` | Zera explicitamente o elemento eliminado, deixando `U` rigorosamente triangular. |

**Observação:** o laço interno começa em `k+1`, e não em `k`, porque `U[i,k]` é tratado separadamente na última linha (e já é usado no cálculo de `fator` antes de ser zerado).

**Verificação:**

```python
erro_LU = np.max(np.abs(np.dot(L, U) - K))
```

- `np.dot(L, U)` multiplica as matrizes (é um produto, não uma inversão ou fatoração).
- `np.abs(...)` toma o módulo de cada elemento da diferença.
- `np.max(...)` pega o maior. Esse é o erro máximo elemento a elemento. Aqui o resultado é `0.0`.

### 6.4 Montagem da inversa

```python
C = np.zeros((n, n))

for coluna in range(n):
    e = np.zeros(n)
    e[coluna] = 1.0
    ...
```

- `C` é inicializada com zeros e será preenchida coluna a coluna.
- `e` é a coluna `coluna` da identidade: zeros com um `1.0` na posição `coluna`.

#### Substituição progressiva (`L y = e`)

```python
y = np.zeros(n)
for i in range(n):
    soma = e[i]
    for j in range(i):
        soma = soma - L[i, j] * y[j]
    y[i] = soma / L[i, i]
```

- Para cada linha `i`, parte-se de `e[i]` e subtrai-se a contribuição dos `y[j]` já conhecidos (`j < i`).
- `range(i)` percorre `j = 0 … i−1`, ou seja, só a parte abaixo da diagonal de `L`.
- A divisão por `L[i,i]` (que vale 1) segue a fórmula geral.
- Como `y[j]` com `j < i` já foi calculado, a ordem de cima para baixo é obrigatória.

#### Substituição regressiva (`U x = y`)

```python
x = np.zeros(n)
x[n - 1] = y[n - 1] / U[n - 1, n - 1]
for i in range(n - 2, -1, -1):
    soma = y[i]
    for j in range(i + 1, n):
        soma = soma - U[i, j] * x[j]
    x[i] = soma / U[i, i]
```

- A última incógnita sai direto: `x[n-1] = y[n-1]/U[n-1,n-1]`.
- `range(n-2, -1, -1)` conta de trás para frente: `3, 2, 1, 0` para `n = 5`.
- Para cada linha `i`, subtrai-se a contribuição dos `x[j]` já conhecidos (`j > i`) e divide-se pelo pivô `U[i,i]`.

#### Armazenamento

```python
for i in range(n):
    C[i, coluna] = x[i]
```

Copia `x` para a coluna `coluna` de `C`. (Equivale a `C[:, coluna] = x`.)

> **Ponto central:** a decomposição LU aparece **antes** do laço `for coluna in range(n)`. Ela é feita uma única vez, e dentro do laço só se resolvem os dois sistemas triangulares.

### 6.5 Verificações

```python
M = np.dot(K, C)
I = np.eye(n)
erromax = np.max(np.abs(M - I))

errosim = np.max(np.abs(C - C.T))
```

- `M = KC` deve ser a identidade. `erromax` mede o desvio.
- `C.T` é a transposta de `C`. `errosim` mede quanto `C` difere de sua transposta, que seria zero se fosse perfeitamente simétrica.

### 6.6 Aplicação: deslocamentos para um carregamento

```python
F = np.array([0.0, 100.0, 0.0, -50.0, 25.0])
u = np.dot(C, F)
```

Depois de ter `C`, calcular os deslocamentos para um novo carregamento é **apenas um produto matriz-vetor**. Não é preciso refazer nenhuma fatoração.

### 6.7 Impressão dos resultados

O programa usa `np.round(L, 4)` e `np.round(C, 8)` apenas para exibir menos casas decimais; os valores armazenados continuam com precisão total. Os deslocamentos são impressos com `f"{u[idx]:.6e}"`, que usa notação científica com 6 casas decimais.

---

## 7. Exemplo numérico acompanhado passo a passo

### 7.1 Fatoração LU

**Passo k = 0 (pivô 2000):**

- `fator = −800 / 2000 = −0,4` e `L[1,0] = −0,4`
- `U[1,1] = 2300 − (−0,4)(−800) = 2300 − 320 = 1980`
- `U[1,2] = −1500` (não muda, pois `U[0,2] = 0`)

**Passo k = 1 (pivô 1980):**

- `fator = −1500 / 1980 ≈ −0,757576` e `L[2,1] ≈ −0,7576`
- `U[2,2] = 2500 − (−0,757576)(−1500) = 2500 − 1136,364 = 1363,636`

**Passo k = 2 (pivô 1363,636):**

- `fator = −1000 / 1363,636 ≈ −0,733333` e `L[3,2] ≈ −0,7333`
- `U[3,3] = 2800 − (−0,733333)(−1000) = 2800 − 733,333 = 2066,667`

**Passo k = 3 (pivô 2066,667):**

- `fator = −1800 / 2066,667 ≈ −0,870968` e `L[4,3] ≈ −0,8710`
- `U[4,4] = 2700 − (−0,870968)(−1800) = 2700 − 1567,742 = 1132,258`

**Resultado:**

$$L=\begin{bmatrix}1&0&0&0&0\\-0{,}4&1&0&0&0\\0&-0{,}7576&1&0&0\\0&0&-0{,}7333&1&0\\0&0&0&-0{,}8710&1\end{bmatrix}\qquad U=\begin{bmatrix}2000&-800&0&0&0\\0&1980&-1500&0&0\\0&0&1363{,}64&-1000&0\\0&0&0&2066{,}67&-1800\\0&0&0&0&1132{,}26\end{bmatrix}$$

Todos os pivôs são positivos, o que confirma que `K` é definida positiva e que não é necessário pivoteamento.

Como `K` é tridiagonal, `L` e `U` são **bidiagonais**: o preenchimento de zeros se preserva.

**Determinante:**

$$\det K = 2000 \times 1980 \times 1363{,}64 \times 2066{,}67 \times 1132{,}26 \approx 1{,}26\times 10^{16} \neq 0$$

### 7.2 Primeira coluna da inversa (j = 1, `e₁ = [1,0,0,0,0]ᵀ`)

**Substituição progressiva (`L y = e₁`):**

| i | Cálculo | y |
|---|---|---|
| 1 | `1` | `1` |
| 2 | `0 − (−0,4)(1)` | `0,4` |
| 3 | `0 − (−0,7576)(0,4)` | `0,30303` |
| 4 | `0 − (−0,7333)(0,30303)` | `0,22222` |
| 5 | `0 − (−0,8710)(0,22222)` | `0,19355` |

**Substituição regressiva (`U x = y`):**

| i | Cálculo | x (m/N) |
|---|---|---|
| 5 | `0,19355 / 1132,258` | `1,7094 × 10⁻⁴` |
| 4 | `(0,22222 + 1800 · 1,7094×10⁻⁴) / 2066,667` | `2,5641 × 10⁻⁴` |
| 3 | `(0,30303 + 1000 · 2,5641×10⁻⁴) / 1363,636` | `4,1026 × 10⁻⁴` |
| 2 | `(0,4 + 1500 · 4,1026×10⁻⁴) / 1980` | `5,1282 × 10⁻⁴` |
| 1 | `(1 + 800 · 5,1282×10⁻⁴) / 2000` | `7,0513 × 10⁻⁴` |

Esse vetor é a **primeira coluna de `C`**, igual ao que o programa imprime.

---

## 8. Saída do programa e resultados

Saída real da execução do código:

```
--- Matriz de Rigidez K ---
[[ 2000.  -800.     0.     0.     0.]
 [ -800.  2300. -1500.     0.     0.]
 [    0. -1500.  2500. -1000.     0.]
 [    0.     0. -1000.  2800. -1800.]
 [    0.     0.     0. -1800.  2700.]]

--- Matriz L ---
[[ 1.      0.      0.      0.      0.    ]
 [-0.4     1.      0.      0.      0.    ]
 [ 0.     -0.7576  1.      0.      0.    ]
 [ 0.      0.     -0.7333  1.      0.    ]
 [ 0.      0.      0.     -0.871   1.    ]]

--- Matriz U ---
[[ 2000.      -800.         0.         0.         0.    ]
 [    0.      1980.     -1500.         0.         0.    ]
 [    0.         0.      1363.6364 -1000.         0.    ]
 [    0.         0.         0.      2066.6667 -1800.    ]
 [    0.         0.         0.         0.      1132.2581]]

Verificação da fatoração (Erro max |LU - K|): 0.0000e+00

--- Matriz de Flexibilidade C (K^-1) ---
[[0.00070513 0.00051282 0.00041026 0.00025641 0.00017094]
 [0.00051282 0.00128205 0.00102564 0.00064103 0.00042735]
 [0.00041026 0.00102564 0.00135385 0.00084615 0.0005641 ]
 [0.00025641 0.00064103 0.00084615 0.00115385 0.00076923]
 [0.00017094 0.00042735 0.0005641  0.00076923 0.00088319]]

Erro máximo da identidade (|KC - I|): 5.5511e-16
Erro de simetria (|C_ij - C_ji|): 1.0842e-19

--- Deslocamentos nodais (u = CF) ---
u1 = 4.273504e-02 m
u2 = 1.068376e-01 m
u3 = 7.435897e-02 m
u4 = 2.564103e-02 m
u5 = 2.635328e-02 m
```

### 8.1 Resumo das verificações

| Verificação | Resultado | Conclusão |
|---|---|---|
| `max｜LU − K｜` | `0` | Fatoração exata |
| `erromax = max｜KC − I｜` | `5,55 × 10⁻¹⁶` | Inversa correta (erro de arredondamento) |
| `errosim = max｜C − Cᵀ｜` | `1,08 × 10⁻¹⁹` | `C` simétrica |

### 8.2 Deslocamentos para `F = [0, 100, 0, −50, 25]ᵀ N`

| Nó | Força (N) | Deslocamento (mm) |
|---|---|---|
| 1 | 0 | 42,735 |
| 2 | 100 | 106,838 |
| 3 | 0 | 74,359 |
| 4 | −50 | 25,641 |
| 5 | 25 | 26,353 |

---

## 9. Interpretação física dos resultados

- **Todos os elementos de `C` são positivos.** Uma força para a direita em qualquer nó desloca todos os nós para a direita, pois os nós são arrastados pelas molas.
- **Cada coluna tem seu máximo na diagonal** (o nó onde a força é aplicada é o que mais se desloca, e os demais decaem conforme se afastam).
- **Maior flexibilidade em `C₃₃ ≈ 1,354 mm/N`:** o nó central é o mais "macio", pois está mais longe das paredes e ligado por molas relativamente moles (`k₂ = 1500` e `k₃ = 1000`).
- **Menor flexibilidade em `C₁₁ ≈ 0,705 mm/N`:** o nó 1 está preso à parede por uma mola rígida (`k₀ = 1200`).
- **Nó 3 sem força direta, mas com 74 mm de deslocamento:** é arrastado pelas forças aplicadas nos nós 2, 4 e 5. É a superposição em ação.
- **Nó 4 com força `−50 N` (para a esquerda), mas deslocando para a direita:** o efeito dos `+100 N` no nó 2 e dos `+25 N` no nó 5 supera o do nó 4.
- **Teste de consistência:** `K·u` reproduz `F`, o que fecha o ciclo `F → u → F`.

---

## 10. Questões de discussão respondidas

### a) Qual é a diferença física entre K e C?

`K` relaciona **deslocamentos a forças** (N/m): mede a resistência da estrutura. `C` relaciona **forças a deslocamentos** (m/N): mede o quanto a estrutura se deforma. Uma estrutura rígida tem `K` grande e `C` pequeno.

### b) O que representa fisicamente a coluna j de C?

Os deslocamentos de todos os nós quando se aplica uma força de 1 N somente no nó `j`. É a "resposta estática ao impulso unitário" da estrutura.

### c) Por que não é necessário repetir a decomposição LU para cada coluna?

Porque `L` e `U` dependem só de `K`, e `K` é a mesma para todas as colunas. O que muda é o lado direito `eⱼ`, tratado pelas substituições.

### d) Qual é a vantagem computacional de reutilizar L e U?

A fatoração (`~n³`) é feita uma vez, e cada coluna adicional custa só `~n²`. O custo total é da ordem de `n³`, e não `n⁴`.

### e) Por que C deve ser simétrica?

Porque `K` é simétrica (a mola exerce forças de mesma intensidade nos dois extremos) e a inversa de uma matriz simétrica é simétrica. Fisicamente, é a reciprocidade de Maxwell–Betti: o efeito de `j` em `i` é igual ao efeito de `i` em `j`.

### f) Se uma das molas se tornar muito pouco rígida, o que se espera nos elementos de C?

Os elementos de `C` **aumentam**: pequenas forças causam grandes deslocamentos, pois a ligação entre partes da estrutura enfraquece. Numericamente, `K` fica **mal condicionada** (o determinante tende a zero) e a inversão perde precisão. No limite `k → 0` com perda total do apoio, `C` diverge.

### g) Se todas as restrições às extremidades fossem removidas, seria possível inverter a matriz de rigidez?

**Não.**

- **Matematicamente:** sem `k₀` e `k₅`, a primeira e a última diagonais ficam `k₁` e `k₄` e **cada linha de `K` soma zero**. Logo `K·[1,1,1,1,1]ᵀ = 0`, com vetor não nulo no núcleo. `K` é singular (`det K = 0`). Na LU, o último pivô `Uₙₙ` sai nulo e a divisão `x[n-1] = y[n-1]/U[n-1,n-1]` falha.
- **Fisicamente:** a estrutura passa a poder se mover como **corpo rígido**: todos os nós podem se deslocar o mesmo valor sem esticar nenhuma mola. Os deslocamentos não são únicos, e se a soma das forças for diferente de zero não existe equilíbrio estático. As paredes eliminam esse movimento.

### h) Depois de conhecida C, qual é a vantagem de utilizá-la para muitos carregamentos?

Cada novo carregamento custa apenas um produto matriz-vetor (`~n²`), sem refatorar nem resolver sistemas. Além disso, `C` permite usar a superposição: qualquer resposta é combinação linear das colunas. Isso é útil em projeto estrutural, com muitas combinações de carga.

---

## 11. Limitações e possíveis extensões

**Limitações do código:**

- **Sem pivoteamento.** Correto para este `K` (simétrica definida positiva), mas se a matriz fosse geral, um pivô zero ou muito pequeno causaria erro ou instabilidade. A solução é a **LU com pivoteamento parcial** (`PK = LU`).
- **Não explora a estrutura tridiagonal.** Há muitas multiplicações por zero. Um algoritmo específico (de Thomas) resolveria cada sistema em `O(n)`.
- **Não verifica singularidade.** Se `U[k,k] == 0`, ocorre divisão por zero. Seria prudente testar antes.
- **Valores fixos no código.** As rigidezes e `F` estão escritos diretamente. Uma versão mais geral usaria uma função `montar_K(k)` para `n` arbitrário.
- **Inverter para um único carregamento não compensa.** Se só interessa um `u`, basta resolver `K u = F` com LU, sem calcular `C` inteira.

**Possíveis extensões:**

1. Transformar os blocos em funções (`decomposicao_lu`, `subst_progressiva`, `subst_regressiva`).
2. Aceitar `n` e as rigidezes como entrada.
3. Implementar pivoteamento parcial.
4. Calcular o número de condição de `K` e estudar o efeito de molas muito moles.
5. Comparar com `np.linalg.inv` apenas como validação (não como solução).
6. Incluir a parte 1 (sistema 2×2) para comparação com o cálculo manual.

---

## 12. Glossário

| Termo | Significado |
|---|---|
| **Grau de liberdade** | Cada deslocamento independente do sistema (aqui, `u₁…u₅`) |
| **Matriz de rigidez K** | Relaciona deslocamentos às forças (`F = Ku`) |
| **Matriz de flexibilidade (complacência) C** | Inversa de `K`; relaciona forças aos deslocamentos (`u = CF`) |
| **Decomposição LU** | Fatoração `K = LU` em triangular inferior e superior |
| **Pivô** | Elemento diagonal `U_kk` usado para eliminar os elementos abaixo dele |
| **Multiplicador (fator)** | Razão `ℓᵢₖ = Uᵢₖ / U_kk`, guardada em `L` |
| **Substituição progressiva** | Resolução de um sistema triangular inferior, de cima para baixo |
| **Substituição regressiva** | Resolução de um sistema triangular superior, de baixo para cima |
| **Matriz simétrica** | `A = Aᵀ` |
| **Definida positiva** | `xᵀAx > 0` para todo `x ≠ 0` |
| **Matriz singular** | Matriz sem inversa (`det = 0`) |
| **Tridiagonal** | Matriz com elementos não nulos apenas na diagonal principal e nas duas adjacentes |
| **Condicionamento** | Medida da sensibilidade da solução a erros nos dados |
| **Reciprocidade de Maxwell–Betti** | Propriedade que torna `C` simétrica |
| **Superposição** | Em sistemas lineares, a resposta total é a soma das respostas individuais |

---

## Estrutura do repositório

```
.
├── inversao_lu.py
└── README.md
```

## Autor

*(adicione seu nome, curso e instituição)*
