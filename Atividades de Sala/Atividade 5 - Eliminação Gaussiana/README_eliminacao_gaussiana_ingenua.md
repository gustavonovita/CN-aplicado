# Eliminação Gaussiana (versão ingênua) em Python

Implementação, em Python + NumPy, do método de **eliminação gaussiana sem pivoteamento** para resolver um sistema linear denso $A\,\mathbf{c}=\mathbf{b}$ de dimensão $5\times5$. O programa também faz a **substituição regressiva**, **conta as operações de ponto flutuante (FLOPs)** e **verifica a solução** pelo resíduo.

Este README serve de material de consulta: explica a **teoria** (o que o algoritmo faz e por quê) e a **programação** (o que cada trecho do código faz), com o exemplo numérico completo, passo a passo.

---

## Sumário

1. [O problema](#1-o-problema)
2. [Teoria](#2-teoria)
   - [2.1 Sistema linear e forma matricial](#21-sistema-linear-e-forma-matricial)
   - [2.2 Ideia da eliminação gaussiana](#22-ideia-da-eliminação-gaussiana)
   - [2.3 Fórmulas da eliminação progressiva](#23-fórmulas-da-eliminação-progressiva)
   - [2.4 Substituição regressiva](#24-substituição-regressiva)
   - [2.5 O pivô e quando o método falha](#25-o-pivô-e-quando-o-método-falha)
   - [2.6 Relação com a fatoração LU e o determinante](#26-relação-com-a-fatoração-lu-e-o-determinante)
   - [2.7 Custo computacional (FLOPs)](#27-custo-computacional-flops)
   - [2.8 Verificação: resíduo e normas](#28-verificação-resíduo-e-normas)
   - [2.9 Erros de arredondamento e condicionamento](#29-erros-de-arredondamento-e-condicionamento)
3. [Exemplo resolvido passo a passo](#3-exemplo-resolvido-passo-a-passo)
4. [O código explicado](#4-o-código-explicado)
5. [Como executar e o que esperar da saída](#5-como-executar-e-o-que-esperar-da-saída)
6. [Limitações e melhorias possíveis](#6-limitações-e-melhorias-possíveis)
7. [Resumo rápido (cola)](#7-resumo-rápido-cola)

---

## 1. O problema

O código resolve o sistema:

$$
\underbrace{\begin{bmatrix}
 6 &  0 & -1 &  0 &  0\\
-3 &  3 &  0 &  0 &  0\\
 0 & -1 &  9 &  0 &  0\\
 0 & -1 & -8 & 11 & -2\\
-3 & -1 &  0 &  0 &  4
\end{bmatrix}}_{A}
\underbrace{\begin{bmatrix}c_1\\c_2\\c_3\\c_4\\c_5\end{bmatrix}}_{\mathbf c}
=
\underbrace{\begin{bmatrix}50\\0\\160\\0\\0\end{bmatrix}}_{\mathbf b}
$$

- As incógnitas são **concentrações** $c_1,\dots,c_5$, impressas em **mg/m³** pelo programa.
- A forma da matriz (diagonal positiva, termos negativos fora da diagonal, lado direito com "cargas" em algumas equações) é a típica de um **balanço de massa em regime permanente** entre elementos conectados (por exemplo, uma série de reatores). Cada linha de $A$ é a equação de balanço de um elemento.
- Os coeficientes vêm do enunciado da atividade. Aqui o foco é o **método numérico**, que serve para qualquer sistema $n\times n$.

**Solução (calculada pelo programa):**

| Incógnita | Valor (mg/m³) |
|---|---|
| $c_1$ | 11,5094 |
| $c_2$ | 11,5094 |
| $c_3$ | 19,0566 |
| $c_4$ | 16,9983 |
| $c_5$ | 11,5094 |

---

## 2. Teoria

### 2.1 Sistema linear e forma matricial

Um sistema de $n$ equações lineares com $n$ incógnitas

$$
a_{i1}x_1+a_{i2}x_2+\dots+a_{in}x_n=b_i,\qquad i=1,\dots,n
$$

escreve-se como $A\mathbf{x}=\mathbf{b}$, onde $A$ é a matriz dos coeficientes $n\times n$, $\mathbf x$ é o vetor de incógnitas e $\mathbf b$ é o vetor dos termos independentes. Neste código o vetor de incógnitas se chama `c` (concentrações).

Uma solução única existe quando $\det A\neq0$ (matriz **não singular**).

### 2.2 Ideia da eliminação gaussiana

Resolver um sistema **triangular superior** é fácil, porque a última equação tem uma só incógnita, a penúltima tem duas (uma já conhecida), e assim por diante. A eliminação gaussiana transforma o sistema original em um triangular equivalente, usando apenas operações que **não alteram a solução**:

> Substituir a linha $i$ por (linha $i$) $-\ m\cdot$ (linha $k$), onde $m$ é um número qualquer.

O algoritmo tem duas fases:

1. **Eliminação progressiva** (*forward elimination*): zera todos os elementos abaixo da diagonal, coluna por coluna, e o resultado é $U\mathbf x=\mathbf{b}'$ com $U$ triangular superior.
2. **Substituição regressiva** (*back substitution*): resolve $U\mathbf x=\mathbf b'$ de baixo para cima.

```
 [ x x x x x | x ]      [ x x x x x | x ]
 [ x x x x x | x ]      [ 0 x x x x | x ]
 [ x x x x x | x ]  →   [ 0 0 x x x | x ]   →  substituição regressiva
 [ x x x x x | x ]      [ 0 0 0 x x | x ]
 [ x x x x x | x ]      [ 0 0 0 0 x | x ]
```

### 2.3 Fórmulas da eliminação progressiva

Na etapa $k$ ($k=1,\dots,n-1$), a linha $k$ é a **linha pivô** e o elemento $a_{kk}$ é o **pivô**. Para cada linha abaixo ($i=k+1,\dots,n$):

**Fator multiplicador**
$$
m_{ik}=\frac{a_{ik}}{a_{kk}}
$$

**Atualização da linha $i$**
$$
a_{ij}\leftarrow a_{ij}-m_{ik}\,a_{kj}\quad(j=k+1,\dots,n),\qquad
b_i\leftarrow b_i-m_{ik}\,b_k
$$

Com esse $m_{ik}$, o novo $a_{ik}$ seria $a_{ik}-m_{ik}a_{kk}=0$ exatamente. Por isso o código **não calcula** a coluna $k$ e só atribui `A[i,k] = 0.0`.

### 2.4 Substituição regressiva

Com $U$ triangular superior (elementos $u_{ij}$) e lado direito $\mathbf b'$:

$$
x_n=\frac{b'_n}{u_{nn}},\qquad
x_i=\frac{b'_i-\sum_{j=i+1}^{n}u_{ij}\,x_j}{u_{ii}},\quad i=n-1,\dots,1
$$

### 2.5 O pivô e quando o método falha

O algoritmo **divide pelo pivô** $a_{kk}$. Isso cria dois problemas:

- **Pivô nulo:** a divisão é impossível. O método para, mesmo que o sistema tenha solução (às vezes basta trocar linhas de lugar).
- **Pivô muito pequeno (em módulo):** o multiplicador $m_{ik}$ fica enorme, os números crescem e os erros de arredondamento se amplificam, ou seja, o algoritmo fica **numericamente instável**.

A solução é o **pivoteamento parcial** (escolher como pivô o maior $|a_{ik}|$ da coluna e trocar as linhas), assunto de uma próxima etapa. A versão deste repositório é **ingênua**: usa sempre o elemento da diagonal.

Neste problema funciona bem: a matriz é **diagonalmente dominante** (em fraca em algumas linhas), e os pivôs obtidos são $6,\;3,\;53/6,\;11,\;4$, todos longe de zero.

### 2.6 Relação com a fatoração LU e o determinante

A eliminação gaussiana equivale a fatorar $A=LU$, onde:

- $U$ é a matriz triangular superior final;
- $L$ é triangular inferior com 1 na diagonal e os multiplicadores $m_{ik}$ abaixo dela.

A atualização de $\mathbf b$ durante a eliminação equivale a resolver $L\mathbf y=\mathbf b$, e a substituição regressiva resolve $U\mathbf x=\mathbf y$.

**Determinante:** $\det A=\prod_k u_{kk}$ (produto dos pivôs, sem troca de linhas). Aqui:

$$
\det A=6\cdot3\cdot\tfrac{53}{6}\cdot11\cdot4=6996
$$

Isso confere com o valor obtido por `numpy.linalg.det`. Como é diferente de zero, a solução existe e é única.

### 2.7 Custo computacional (FLOPs)

O código conta separadamente **multiplicações/divisões** e **somas/subtrações** (sem contar comparações, atribuições e laços). Para uma matriz densa $n\times n$:

| Fase | Mult/Div | Soma/Sub |
|---|---|---|
| Eliminação | $\displaystyle\sum_{m=1}^{n-1}m(m+2)=\frac{n(n-1)(2n+5)}{6}$ | $\displaystyle\sum_{m=1}^{n-1}m(m+1)=\frac{(n-1)n(n+1)}{3}$ |
| Substituição regressiva | $\dfrac{n(n+1)}{2}$ | $\dfrac{n(n-1)}{2}$ |

onde $m=n-k$ é o número de linhas abaixo do pivô. A contagem por linha eliminada $i$ é: 1 divisão (fator) + $(n-k)$ multiplicações (laço em $j$) + 1 multiplicação (para $b$), e $(n-k)+1$ subtrações.

Para $n=5$:

| | Mult/Div | Soma/Sub | Total |
|---|---|---|---|
| Eliminação | 50 | 40 | 90 |
| Substituição | 15 | 10 | 25 |
| **Total** | **65** | **50** | **115** |

Para $n$ grande, o custo domina em $\approx\frac{n^3}{3}$ multiplicações e $\frac{n^3}{3}$ somas, isto é, $\mathcal O(n^3)$, enquanto a substituição regressiva custa apenas $\mathcal O(n^2)$.

> O código é **genérico**: executa todas as operações mesmo quando o multiplicador é 0 (por exemplo, $m_{31}=0$ na primeira etapa). Por isso a contagem é sempre 65 + 50 para $n=5$, independentemente da esparsidade da matriz.

### 2.8 Verificação: resíduo e normas

Como a solução exata não é fornecida ao algoritmo, a verificação é feita pelo **vetor resíduo**, calculado com a matriz e o vetor **originais**:

$$
\mathbf r=A_{orig}\,\mathbf c_{calc}-\mathbf b_{orig},\qquad
\|\mathbf r\|_\infty=\max_i|r_i|
$$

Se $\mathbf c_{calc}$ fosse exato, $\mathbf r=\mathbf 0$. Na prática, o resíduo fica na ordem do arredondamento de máquina ($\approx10^{-15}$ a $10^{-14}$ em `float64`).

É importante usar $A_{orig}$ e $b_{orig}$ (e não os modificados pela eliminação). Por isso o código guarda cópias.

### 2.9 Erros de arredondamento e condicionamento

- Em `float64`, o epsilon de máquina é $\varepsilon\approx2{,}2\times10^{-16}$. Cada operação carrega um erro relativo desse tamanho.
- O **número de condição** $\kappa(A)=\|A\|\,\|A^{-1}\|$ mede quanto o problema amplifica erros. Vale a relação

$$
\frac{\|\mathbf e\|}{\|\mathbf x\|}\;\lesssim\;\kappa(A)\,\frac{\|\mathbf r\|}{\|\mathbf b\|}
$$

onde $\mathbf e$ é o erro da solução. Aqui $\kappa_\infty(A)\approx11{,}6$ (**bem condicionado**), então um resíduo de $\sim10^{-14}$ implica erro da solução também minúsculo. Em sistemas mal condicionados, um resíduo pequeno **não** garante solução precisa.

- **Condicionamento** é propriedade do problema; **estabilidade** é propriedade do algoritmo. A eliminação sem pivoteamento pode ser instável mesmo num problema bem condicionado, se um pivô for pequeno (não é o caso aqui).

---

## 3. Exemplo resolvido passo a passo

Aqui uso a numeração matemática (linhas e colunas de 1 a 5). O código usa `k = 0..3`, ou seja, **etapa matemática $k$ = `k` do código + 1**.

### Etapa 1: pivô $a_{11}=6$

| Linha | Multiplicador $m_{i1}$ | Resultado |
|---|---|---|
| 2 | $-3/6=-0{,}5$ | $a_{22}=3,\;a_{23}=-0{,}5,\;b_2=25$ |
| 3 | $0/6=0$ | linha inalterada |
| 4 | $0/6=0$ | linha inalterada |
| 5 | $-3/6=-0{,}5$ | $a_{52}=-1,\;a_{53}=-0{,}5,\;b_5=25$ |

Exemplo da conta da linha 2: $a_{23}=0-(-0{,}5)(-1)=-0{,}5$ e $b_2=0-(-0{,}5)(50)=25$.

### Etapa 2: pivô $a_{22}=3$

| Linha | $m_{i2}$ | Resultado |
|---|---|---|
| 3 | $-1/3$ | $a_{33}=9-\frac{1}{6}=8{,}8333,\;b_3=160+\frac{25}{3}=168{,}3333$ |
| 4 | $-1/3$ | $a_{43}=-8-\frac16=-8{,}1667,\;b_4=8{,}3333$ |
| 5 | $-1/3$ | $a_{53}=-0{,}5-\frac16=-0{,}6667,\;b_5=33{,}3333$ |

### Etapa 3: pivô $a_{33}=53/6\approx8{,}8333$

| Linha | $m_{i3}$ | Resultado |
|---|---|---|
| 4 | $-\dfrac{49}{53}\approx-0{,}9245$ | $a_{44}=11,\;a_{45}=-2,\;b_4=163{,}9623$ |
| 5 | $-\dfrac{4}{53}\approx-0{,}0755$ | $a_{55}=4,\;b_5=46{,}0377$ |

### Etapa 4: pivô $a_{44}=11$

$m_{54}=0/11=0$, então nada muda (o código ainda executa as operações e as conta como FLOPs).

### Matriz triangular obtida

$$
U=\begin{bmatrix}
6&0&-1&0&0\\
0&3&-0{,}5&0&0\\
0&0&8{,}8333&0&0\\
0&0&0&11&-2\\
0&0&0&0&4
\end{bmatrix},\qquad
\mathbf b'=\begin{bmatrix}50\\25\\168{,}3333\\163{,}9623\\46{,}0377\end{bmatrix}
$$

### Substituição regressiva

$$
\begin{aligned}
c_5&=\tfrac{46{,}0377}{4}=11{,}5094\\
c_4&=\tfrac{163{,}9623+2\cdot11{,}5094}{11}=16{,}9983\\
c_3&=\tfrac{168{,}3333}{8{,}8333}=19{,}0566\\
c_2&=\tfrac{25+0{,}5\cdot19{,}0566}{3}=11{,}5094\\
c_1&=\tfrac{50+19{,}0566}{6}=11{,}5094
\end{aligned}
$$

### Conferência manual nas equações originais

| Equação | Conta | Resultado |
|---|---|---|
| 1 | $6(11{,}5094)-19{,}0566$ | $50{,}00$ ✓ |
| 2 | $-3(11{,}5094)+3(11{,}5094)$ | $0$ ✓ |
| 3 | $-11{,}5094+9(19{,}0566)$ | $160{,}00$ ✓ |
| 5 | $-3(11{,}5094)-11{,}5094+4(11{,}5094)$ | $0$ ✓ |

---

## 4. O código explicado

O código segue as "partes" do roteiro da atividade (comentários `Parte 1`, `Parte 3`, `Parte 4`, `Parte 6`).

### 4.1 Importação e dados (Parte 1)

```python
import numpy as np

n = 5

A_orig = np.array([
    [ 6.0,  0.0, -1.0,  0.0,  0.0],
    [-3.0,  3.0,  0.0,  0.0,  0.0],
    ...
])
b_orig = np.array([50.0, 0.0, 160.0, 0.0, 0.0])
```

- `np.array` cria uma matriz (2D) e um vetor (1D) de `float64`.
- Os valores são escritos com `.0` para garantir ponto flutuante. Com inteiros a divisão `A[i,k]/A[k,k]` também daria float em Python 3, mas é boa prática deixar explícito.
- `n` é a dimensão do sistema, usada em todos os laços.

### 4.2 Cópias de trabalho

```python
A = A_orig.copy()
b = b_orig.copy()
```

A eliminação **modifica** `A` e `b`. Em NumPy, `A = A_orig` não copiaria os dados, apenas criaria outro nome para o mesmo objeto. Por isso `.copy()` é essencial: `A_orig` e `b_orig` permanecem intactos para calcular o resíduo no fim.

### 4.3 Contadores de FLOPs

```python
cont_mult_div = 0
cont_soma_sub = 0
```

Dois inteiros, incrementados logo após cada operação aritmética correspondente. Não se contam comparações (`if abs(...) < ...`), atribuições (`A[i,k] = 0.0`) nem a verificação de pivô.

### 4.4 Eliminação progressiva (Parte 3)

```python
for k in range(n - 1):
    if abs(A[k, k]) < 1e-12:
        print(f"Erro: Pivô nulo detectado na etapa k={k}. Execução interrompida.")
        break

    for i in range(k + 1, n):
        fator = A[i, k] / A[k, k]
        cont_mult_div += 1

        for j in range(k + 1, n):
            A[i, j] = A[i, j] - fator * A[k, j]
            cont_mult_div += 1
            cont_soma_sub += 1

        b[i] = b[i] - fator * b[k]
        cont_mult_div += 1
        cont_soma_sub += 1

        A[i, k] = 0.0
```

**Mapeamento com a teoria (indexação 0 do Python):**

| Código | Significado |
|---|---|
| `range(n - 1)` → `k = 0,…,3` | etapas de eliminação: as colunas 1 a $n-1$ (a última coluna não tem nada abaixo da diagonal) |
| `A[k, k]` | pivô $a_{kk}$ |
| `range(k + 1, n)` → `i` | linhas abaixo do pivô |
| `fator` | multiplicador $m_{ik}=a_{ik}/a_{kk}$ |
| laço em `j` de `k+1` a `n-1` | atualiza a linha $i$ da coluna $k+1$ em diante |
| `b[i] = b[i] - fator*b[k]` | aplica a mesma operação de linha ao termo independente |
| `A[i, k] = 0.0` | grava o zero que o cálculo produziria (veja abaixo) |

**Detalhes de programação importantes:**

1. **Ordem das operações:** `fator` precisa ser calculado **antes** de zerar `A[i,k]`, porque usa esse valor. Por isso `A[i,k] = 0.0` fica no **fim** do laço de `i`.
2. **Por que `j` começa em `k+1`:** a coluna $k$ não é calculada por subtração. Ela vale zero por construção, e isso economiza operações e evita resíduos de arredondamento (como $10^{-17}$) abaixo da diagonal.
3. **Verificação de pivô:** `abs(A[k,k]) < 1e-12` protege contra divisão por (quase) zero, com tolerância **absoluta** fixa de $10^{-12}$.
4. **Contadores:** por linha eliminada, são 1 divisão + $(n-k-1)$ multiplicações no laço de `j` + 1 multiplicação em `b`, e $(n-k-1)+1$ subtrações.

**Impressão a cada etapa:**

```python
print(f"\nApós pivô k = {k}:")
print(np.round(A, 4))
print(np.round(b, 4))
```

`np.round(A, 4)` apenas **formata** a saída com 4 casas, e `A` continua com precisão total internamente.

### 4.5 Substituição regressiva (Parte 4)

```python
c = np.zeros(n)

c[n - 1] = b[n - 1] / A[n - 1, n - 1]
cont_mult_div += 1

for i in range(n - 2, -1, -1):
    soma = b[i]
    for j in range(i + 1, n):
        soma = soma - A[i, j] * c[j]
        cont_mult_div += 1
        cont_soma_sub += 1

    c[i] = soma / A[i, i]
    cont_mult_div += 1
```

- `np.zeros(n)` cria o vetor solução, que é preenchido de trás para frente.
- `range(n - 2, -1, -1)` percorre `i = n-2, n-3, …, 0`.
- `soma` começa em `b[i]` e vai subtraindo as contribuições das incógnitas já conhecidas (`A[i,j]*c[j]`, com `j > i`). O que sobra é dividido pelo elemento diagonal.
- Corresponde diretamente às fórmulas da seção [2.4](#24-substituição-regressiva).

### 4.6 Saída da solução e do custo

```python
for idx in range(n):
    print(f"c{idx+1} = {c[idx]:.4f} mg/m³")
```

- `f"..."` é uma *f-string*: o que está entre `{}` é calculado e inserido.
- `{c[idx]:.4f}` imprime em ponto fixo com 4 casas decimais.
- `idx+1` converte o índice 0 do Python para o nome matemático $c_1,\dots,c_5$.

O custo total é `cont_mult_div + cont_soma_sub` (115 para $n=5$).

### 4.7 Verificação (Parte 6)

```python
r = np.dot(A_orig, c) - b_orig
norma_inf_r = np.max(np.abs(r))
print(f"Norma infinito do resíduo (||r||_inf): {norma_inf_r:.4e}")
```

- `np.dot(A_orig, c)` é o produto matriz-vetor $A_{orig}\mathbf c$.
- `np.abs(r)` aplica módulo em cada componente, e `np.max` pega o maior, que é a norma infinito.
- `:.4e` imprime em notação científica (ex.: `7.1054e-15`).
- Usa-se `A_orig`/`b_orig` para provar que a solução satisfaz o sistema **original**.

---

## 5. Como executar e o que esperar da saída

Requisitos: Python 3 e NumPy.

```bash
pip install numpy
python eliminacao_gaussiana.py
```

Saída do programa, na ordem:

1. **Sistema original** (`A` e `b`).
2. **Quatro blocos** "Após pivô k = 0…3", cada um com `A` e `b` após a etapa. O primeiro bloco mostra as linhas 2 e 5 modificadas, e o último mostra a matriz triangular completa.
3. **Matriz triangular superior resultante.**
4. **Solução**, com os valores abaixo:

```
c1 = 11.5094 mg/m³
c2 = 11.5094 mg/m³
c3 = 19.0566 mg/m³
c4 = 16.9983 mg/m³
c5 = 11.5094 mg/m³
```

5. **Custo computacional:**

```
Multiplicações/Divisões: 65
Somas/Subtrações: 50
Total de FLOPs: 115
```

6. **Verificação:** $\|\mathbf r\|_\infty\approx7{,}1\times10^{-15}$ (só nas componentes 4 e 5; as demais são exatamente 0). É resíduo de arredondamento, e a solução está correta.

**Conferência independente:** `numpy.linalg.solve(A_orig, b_orig)` retorna os mesmos valores. O determinante é 6996 e $\kappa_\infty(A)\approx11{,}6$.

---

## 6. Limitações e melhorias possíveis

| Ponto | Comentário |
|---|---|
| **Sem pivoteamento** | Funciona aqui porque os pivôs são grandes. Para pivôs pequenos ou nulos, é preciso pivoteamento parcial (busca do maior $|a_{ik}|$ na coluna e troca de linhas). |
| **`break` não encerra o programa** | Se um pivô for pequeno, o `break` só sai do laço da eliminação, e a substituição regressiva roda mesmo assim com uma matriz inconsistente, produzindo resultado sem sentido. O ideal é parar o programa (`sys.exit()`, `raise` ou uma variável de controle). |
| **Último pivô não é testado** | O laço vai até `k = n-2`; `A[n-1, n-1]` também é divisor na substituição e deveria ser verificado. |
| **Tolerância absoluta fixa (`1e-12`)** | Depende da escala dos dados. Uma tolerância **relativa**, como $n\,\varepsilon\max|A_{ij}|$, é mais robusta. |
| **Comentários `[cite: 3]`** | Sobras do texto-fonte nos comentários do código. Convém apagá-los antes de publicar. |
| **Não explora esparsidade** | Intencional (o enunciado pede um solucionador genérico), mas o custo $\mathcal O(n^3)$ seria evitável para matrizes esparsas. |
| **Dados fixos no código** | Para reaproveitar, vale encapsular em funções (`eliminacao(A, b)`, `subst_regressiva(U, b)`) e ler `A`, `b` como parâmetros. |

Exemplo de verificação mais segura do último pivô:

```python
if abs(A[n - 1, n - 1]) < 1e-12:
    raise ZeroDivisionError("Pivô nulo na última etapa")
```

---

## 7. Resumo rápido (cola)

- **Objetivo:** resolver $A\mathbf c=\mathbf b$ transformando $A$ em triangular superior e resolvendo de baixo para cima.
- **Eliminação:** para cada pivô $a_{kk}$, para cada linha $i>k$: $m=a_{ik}/a_{kk}$; linha$_i\leftarrow$ linha$_i-m\cdot$linha$_k$ (inclusive em $b$).
- **Substituição regressiva:** $x_i=(b_i-\sum_{j>i}u_{ij}x_j)/u_{ii}$.
- **Pivôs deste problema:** $6,\,3,\,53/6,\,11,\,4$ → $\det A=6996$.
- **Solução:** $c_1=c_2=c_5\approx11{,}5094$; $c_3\approx19{,}0566$; $c_4\approx16{,}9983$ mg/m³.
- **Custo ($n=5$):** 65 mult/div + 50 soma/sub = **115 FLOPs**.
- **Verificação:** $\|\mathbf r\|_\infty\approx7\times10^{-15}$; $\kappa_\infty(A)\approx11{,}6$ (bem condicionado).
- **Cuidado:** sem pivoteamento, pivôs pequenos ou nulos quebram o método. Esse é o assunto da versão com pivoteamento parcial.
