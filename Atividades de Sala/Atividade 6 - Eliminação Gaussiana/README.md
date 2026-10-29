# Eliminação Gaussiana com Pivoteamento Parcial: estrutura de barras

Implementação em Python (NumPy) da **eliminação gaussiana com pivoteamento parcial** aplicada ao cálculo dos esforços axiais de uma treliça plana de 5 barras. O projeto faz parte da disciplina **Cálculo Numérico Aplicado** e tem como foco distinguir dois conceitos:

| Conceito | O que é | Quem "tem" essa propriedade |
|---|---|---|
| **Condicionamento** | Sensibilidade da solução a pequenas perturbações nos dados | O **problema** (a geometria das barras) |
| **Estabilidade numérica** | Quanto os erros de arredondamento crescem durante o cálculo | O **algoritmo** (a eliminação gaussiana) |

Este documento foi escrito para consulta: explica a teoria, depois o código linha por linha, depois os resultados e por fim as limitações e melhorias possíveis.

---

## Sumário

1. [Como executar](#1-como-executar)
2. [Problema físico](#2-problema-físico)
3. [Teoria matemática](#3-teoria-matemática)
4. [Explicação do código](#4-explicação-do-código)
5. [Resultados obtidos e interpretação](#5-resultados-obtidos-e-interpretação)
6. [Limitações e melhorias](#6-limitações-e-melhorias)
7. [Glossário](#7-glossário)
8. [Referências](#8-referências)

---

## 1. Como executar

**Requisitos:** Python 3.8 ou superior e NumPy.

```bash
pip install numpy
python3 eliminacao_gaussiana.py
```

> Salve o arquivo com extensão `.py`. O arquivo original está sem extensão (`Eliminação_Gaussiana`), e acentos em nomes de arquivo costumam causar problemas em terminais e no `import`. Sugestão de nome: `eliminacao_gaussiana.py`.

**Estrutura do arquivo:**

```
eliminacao_gaussiana.py
├── criar_sistema_barras()               # monta A, b e F_ref
├── eliminacao_gaussiana_pivoteamento()  # resolve A·F = b
├── rodar_cenario()                      # executa e avalia um cenário
└── bloco principal (__main__)           # cenários A e C + teste de sensibilidade
```

---

## 2. Problema físico

### 2.1 Treliça de 5 barras

Uma estrutura plana é formada por barras articuladas, sob carga estática. Cada barra transmite **apenas esforço axial**:

- $F_i > 0$: barra **tracionada**;
- $F_i < 0$: barra **comprimida**.

Há três nós (A, B, C):

- **Nó A:** barras 1 e 2 ligam A a dois apoios fixos, formando ângulos $\theta_1$ e $\theta_2$ com a vertical; a barra 3 liga A a B na horizontal.
- **Nó B:** a barra 4 liga B a um apoio fixo, a um ângulo $\phi$ abaixo da horizontal; a barra 5 liga B a C na vertical.
- **Nó C:** possui um guia que fornece reação horizontal, mas **não** fornece reação vertical.

### 2.2 Equilíbrio de forças

Em cada nó, a soma das forças em x e em y deve ser zero. Isso gera as cinco equações:

$$
\begin{aligned}
F_1\sin\theta_1 + F_2\sin\theta_2 + F_3 &= -P_{Ax} &&\text{(nó A, horizontal)}\\
F_1\cos\theta_1 + F_2\cos\theta_2 &= -P_{Ay} &&\text{(nó A, vertical)}\\
-F_3 + F_4\cos\phi &= -P_{Bx} &&\text{(nó B, horizontal)}\\
-F_4\sin\phi - F_5 &= -P_{By} &&\text{(nó B, vertical)}\\
F_5 &= -P_{Cy} &&\text{(nó C, vertical)}
\end{aligned}
$$

A reação horizontal do guia em C **não é incógnita**, por isso o nó C contribui com uma única equação. São 5 incógnitas ($F_1,\dots,F_5$) e 5 equações.

### 2.3 Forma matricial

$$
\underbrace{\begin{bmatrix}
\sin\theta_1 & \sin\theta_2 & 1 & 0 & 0\\
\cos\theta_1 & \cos\theta_2 & 0 & 0 & 0\\
0 & 0 & -1 & \cos\phi & 0\\
0 & 0 & 0 & -\sin\phi & -1\\
0 & 0 & 0 & 0 & 1
\end{bmatrix}}_{A}
\underbrace{\begin{bmatrix}F_1\\F_2\\F_3\\F_4\\F_5\end{bmatrix}}_{F}
=
\underbrace{-\begin{bmatrix}P_{Ax}\\P_{Ay}\\P_{Bx}\\P_{By}\\P_{Cy}\end{bmatrix}}_{b}
$$

A **ordem das equações é intencional**: quando $\theta_1$ é muito pequeno, o elemento $A_{11}=\sin\theta_1$ (o primeiro pivô) também é muito pequeno. Isso cria um caso de teste para o pivoteamento.

### 2.4 Por que barras quase paralelas são problemáticas

As duas primeiras colunas da matriz, $(\sin\theta_1,\cos\theta_1)$ e $(\sin\theta_2,\cos\theta_2)$, são os vetores unitários das direções das barras 1 e 2. Se as barras ficam quase paralelas, as colunas ficam quase linearmente dependentes e a matriz quase singular.

Fisicamente, duas barras alinhadas só conseguem resistir a forças na direção delas. Para equilibrar qualquer força com componente transversal, os esforços nas barras precisam crescer muito. Por isso uma pequena variação na carga pode causar uma grande variação nos esforços.

No bloco 2×2 formado pelas barras 1 e 2:

$$
\det\begin{bmatrix}\sin\theta_1 & \sin\theta_2\\ \cos\theta_1 & \cos\theta_2\end{bmatrix}
= \sin\theta_1\cos\theta_2 - \sin\theta_2\cos\theta_1 = \sin(\theta_1-\theta_2)
$$

A última igualdade é a identidade trigonométrica do seno da diferença. Quando $\theta_1-\theta_2\to 0$, o determinante vai a zero.

---

## 3. Teoria matemática

### 3.1 Eliminação gaussiana

Para resolver $A x = b$, o método transforma o sistema em um sistema **triangular superior** equivalente, que é fácil de resolver. Tem duas fases.

**Fase 1: eliminação progressiva.** Para cada coluna $k = 1,\dots,n-1$, zera-se tudo abaixo da diagonal subtraindo múltiplos da linha $k$ das linhas seguintes:

$$
m_{ik} = \frac{a_{ik}}{a_{kk}}, \qquad
a_{ij} \leftarrow a_{ij} - m_{ik}\,a_{kj}, \qquad
b_i \leftarrow b_i - m_{ik}\,b_k
\qquad (i,j>k)
$$

- $a_{kk}$ é o **pivô** da etapa $k$.
- $m_{ik}$ é o **fator de eliminação** (multiplicador).

Subtrair um múltiplo de uma equação de outra **não altera a solução** do sistema. Ao final, a matriz $U$ é triangular superior.

**Fase 2: substituição regressiva.** Resolve-se de baixo para cima:

$$
x_n = \frac{b_n}{u_{nn}}, \qquad
x_i = \frac{b_i - \sum_{j=i+1}^{n} u_{ij}\,x_j}{u_{ii}} \quad (i = n-1,\dots,1)
$$

**Custo:** aproximadamente $\tfrac{2}{3}n^3$ operações de ponto flutuante na eliminação, ou seja, $O(n^3)$. A substituição regressiva custa só $O(n^2)$. Para $n=5$ são 65 multiplicações/divisões e 50 somas/subtrações (considerando a matriz densa, sem aproveitar zeros).

### 3.2 Por que o pivô importa

O fator $m_{ik}=a_{ik}/a_{kk}$ tem o pivô no denominador. Se $|a_{kk}|$ for muito menor que os elementos abaixo dele, $|m_{ik}|$ será enorme. Os números grandes subtraem números de tamanho normal, e as informações originais ficam "engolidas" pelo arredondamento.

**Exemplo clássico** (aritmética de ponto flutuante dupla, $\varepsilon = 10^{-20}$):

$$
\begin{bmatrix}\varepsilon & 1\\ 1 & 1\end{bmatrix}
\begin{bmatrix}x_1\\x_2\end{bmatrix} =
\begin{bmatrix}1\\2\end{bmatrix}
\quad\Rightarrow\quad x_1\approx 1,\; x_2\approx 1
$$

- **Sem pivoteamento:** $m = 1/\varepsilon = 10^{20}$; $u_{22} = 1 - 10^{20} \approx -10^{20}$; $b_2' = 2 - 10^{20}\approx -10^{20}$. Então $x_2 = 1$ e $x_1 = (1 - 1)/\varepsilon = 0$. **Errado.** O "1" da primeira equação se perdeu na soma.
- **Com pivoteamento** (trocar as linhas): $m=\varepsilon$, e as operações ficam bem comportadas. Resultado: $x_1 = x_2 = 1$.

### 3.3 Pivoteamento parcial

Na etapa $k$, entre as linhas $k,k+1,\dots,n$, escolhe-se a linha $p$ com o **maior valor absoluto** na coluna $k$:

$$
p = \arg\max_{k\le i\le n}\,\lvert a_{ik}\rvert
$$

e troca-se a linha $k$ com a linha $p$. Isso garante $\lvert m_{ik}\rvert \le 1$ para todo $i>k$, o que mantém os números sob controle.

Pontos importantes:

- A troca deve ser feita na **linha inteira**, em todas as colunas de $A$ **e** no elemento correspondente de $b$. Trocar só parte da linha destruiria a equação.
- Trocar linhas **apenas reordena as equações**. O sistema físico e a solução exata não mudam.
- Matricialmente, é equivalente a resolver $PA\,x = Pb$, onde $P$ é uma matriz de permutação (resultando na fatoração $PA = LU$).
- É "parcial" porque se procura o pivô só na coluna (trocando linhas). O pivoteamento *total* procuraria também nas colunas (trocando colunas), o que é mais caro e raramente necessário.

### 3.4 Pivô numericamente nulo

Se o maior $\lvert a_{pk}\rvert$ da coluna for (numericamente) zero, a matriz é singular ou quase singular, e a eliminação não pode continuar. O código interrompe a execução e avisa.

Em ponto flutuante, "zero" quase nunca é exatamente zero (aparecem valores como $10^{-17}$ por arredondamento). Por isso se compara com um **limiar** (`limite_pivo`). O código usa $10^{-12}$ (ver a [seção 6](#6-limitações-e-melhorias) sobre as limitações de um limiar fixo).

### 3.5 Normas e a norma infinito

Para medir o "tamanho" de um vetor usa-se uma norma. A **norma infinito** é o maior valor absoluto das componentes:

$$
\lVert v\rVert_\infty = \max_i \lvert v_i\rvert
$$

Para matrizes, a norma infinito induzida é a maior soma de valores absolutos por linha:

$$
\lVert A\rVert_\infty = \max_i \sum_j \lvert a_{ij}\rvert
$$

### 3.6 Número de condição

O **número de condição** de uma matriz inversível é

$$
\kappa(A) = \lVert A\rVert\,\lVert A^{-1}\rVert \;\ge 1
$$

Ele limita quanto uma perturbação relativa em $b$ pode ser amplificada na solução. Se $A x = b$ e $A(x+\delta x) = b+\delta b$, então $\delta x = A^{-1}\delta b$, e como $\lVert b\rVert\le\lVert A\rVert\lVert x\rVert$:

$$
\boxed{\frac{\lVert\delta x\rVert}{\lVert x\rVert}\;\le\;\kappa(A)\,\frac{\lVert\delta b\rVert}{\lVert b\rVert}}
$$

Interpretação:

- $\kappa\approx 1$: **bem condicionado**. Perturbações relativas não são amplificadas.
- $\kappa \gg 1$: **mal condicionado**. Perturbações pequenas podem ser amplificadas por até $\kappa$ vezes (no pior caso).
- Uma regra prática: perde-se cerca de $\log_{10}\kappa$ algarismos decimais de precisão na solução.

Os valores deste projeto (norma infinito):

| Cenário | θ₁ | θ₂ | κ∞(A) |
|---|---|---|---|
| A (afastadas) | 15° | 65° | **11,65** |
| C (quase paralelas) | 0,01° | 0,02° | **4,58·10⁴** |

**Determinante pequeno não é o mesmo que mal condicionado.** O determinante depende da escala. $A=10^{-3}I$ tem determinante minúsculo, mas $\kappa=1$. Trocar a unidade de uma equação (N por kN) multiplica o determinante, e a sensibilidade física não muda. O número de condição é invariante a essas escalas, por isso é a medida correta.

**O condicionamento é o *pior caso*.** Dependendo da direção da perturbação, a amplificação real pode ser bem menor que $\kappa$. Veremos isso na [seção 5.3](#53-teste-de-sensibilidade).

### 3.7 Erro × resíduo

Seja $\hat F$ a solução calculada e $F$ a solução exata.

| Grandeza | Definição | O que mede |
|---|---|---|
| **Erro** | $e = \hat F - F$ | Distância até a solução exata (só se conhece se a solução exata for conhecida) |
| **Resíduo** | $r = A\hat F - b$ | Quão bem $\hat F$ satisfaz as equações (sempre calculável) |

A relação entre eles envolve o condicionamento:

$$
\frac{\lVert e\rVert}{\lVert F\rVert}\;\le\;\kappa(A)\,\frac{\lVert r\rVert}{\lVert b\rVert}
$$

**Consequência:** em um problema mal condicionado, um resíduo pequeno **não garante** um erro pequeno. No cenário C o resíduo é ~$10^{-15}$ e o erro ~$10^{-11}$, a diferença é de quatro ordens de grandeza.

### 3.8 Estabilidade × condicionamento

- Um **algoritmo estável** (backward stable) entrega a solução *exata* de um problema ligeiramente perturbado, o que se manifesta como resíduo pequeno. A eliminação gaussiana com pivoteamento parcial é estável na prática.
- O erro final é aproximadamente: **(condicionamento do problema) × (erro de arredondamento do algoritmo)**.
- Um bom algoritmo não "cura" um problema mal condicionado. Um mau algoritmo (eliminação sem pivô, com pivôs minúsculos) pode estragar até um problema bem condicionado.

Por isso **a troca de linhas não elimina a sensibilidade física** das barras quase paralelas: $\kappa(PA)=\kappa(A)$, já que permutar linhas não altera as normas de $A$ nem de $A^{-1}$.

### 3.9 O que acontece se θ₁ = θ₂

As colunas 1 e 2 ficam idênticas, a matriz é singular ($\det=0$) e o sistema não tem solução única (ou nenhuma solução, ou infinitas). Fisicamente, a estrutura é um mecanismo. No programa, o pivô da etapa 2 seria nulo e a execução seria interrompida.

---

## 4. Explicação do código

> **Convenção de índices:** o pseudocódigo do enunciado é em estilo Fortran (índices começam em 1). Python começa em 0. Então `for k in range(n-1)` equivale a `do k = 1, n-1`.

### 4.1 Importação

```python
import numpy as np
```

O NumPy fornece arrays, funções trigonométricas, produto matriz-vetor e operações vetorizadas.

### 4.2 `criar_sistema_barras(theta1_deg, theta2_deg, phi_deg=45.0)`

Monta o sistema do problema.

```python
th1 = np.radians(theta1_deg)
th2 = np.radians(theta2_deg)
phi = np.radians(phi_deg)
```

As funções `np.sin` e `np.cos` trabalham em **radianos**. Os ângulos chegam em graus e são convertidos com `np.radians`, que equivale a multiplicar por $\pi/180$. Esquecer essa conversão é um erro muito comum.

```python
A = np.array([
    [np.sin(th1), np.sin(th2),  1.0,        0.0,  0.0],
    [np.cos(th1), np.cos(th2),  0.0,        0.0,  0.0],
    [        0.0,         0.0, -1.0, np.cos(phi),  0.0],
    [        0.0,         0.0,  0.0,-np.sin(phi), -1.0],
    [        0.0,         0.0,  0.0,        0.0,  1.0]
])
```

É exatamente a matriz $A$ da [seção 2.3](#23-forma-matricial). Cada linha é uma equação de equilíbrio e cada coluna corresponde a uma barra. Usar `1.0` e `0.0` (e não inteiros) garante que o array seja de ponto flutuante.

```python
F_ref = np.array([100.0, -80.0, 20.0, 35.0, 15.0])
b = np.dot(A, F_ref)
```

Define os esforços de referência $F_{ref}$ e calcula $b = A F_{ref}$. Assim se obtém um conjunto de cargas **fisicamente consistentes** com uma solução conhecida. O algoritmo só recebe `A` e `b`, e `F_ref` serve **apenas para conferir** a resposta. Retorna `A, b, F_ref`.

### 4.3 `eliminacao_gaussiana_pivoteamento(A_orig, b_orig, limite_pivo=1e-12)`

O solucionador. Recebe a matriz, o vetor e o limiar de pivô, e devolve o vetor solução (ou `None` em caso de falha).

**Cópias de segurança**

```python
A = A_orig.copy()
b = b_orig.copy()
n = len(b)
```

O algoritmo modifica `A` e `b` durante a eliminação. Trabalhar em cópias preserva `A_orig` e `b_orig`, que são necessários depois para calcular o resíduo $r = A_{orig}\hat F - b_{orig}$. Sem o `.copy()`, os arrays originais seriam destruídos (em Python, atribuir um array não o duplica, só cria outro nome para ele).

**Laço externo: uma etapa por coluna**

```python
for k in range(n - 1):
```

Percorre as colunas $0,\dots,n-2$. A última coluna não precisa de eliminação, pois não há linhas abaixo da última.

**Etapa 1: encontrar a linha do pivô**

```python
p = k
for i in range(k + 1, n):
    if abs(A[i, k]) > abs(A[p, k]):
        p = i
```

Começa assumindo que a linha pivô é a própria $k$. Percorre as linhas abaixo e, sempre que encontra um elemento de maior valor absoluto na coluna $k$, atualiza `p`. Ao final, `p` é o índice do maior $\lvert a_{ik}\rvert$ entre as linhas $k,\dots,n-1$. Esta busca ocorre **uma vez por etapa $k$**, fora do laço de eliminação, como exige o enunciado.

**Etapa 2: teste de pivô nulo**

```python
if abs(A[p, k]) <= limite_pivo:
    print(f"Erro: Pivô nulo/muito pequeno na etapa {k}. Execução interrompida.")
    return None
```

Se o maior elemento da coluna é praticamente zero, a matriz é singular na precisão adotada e a eliminação não pode prosseguir. A função informa o erro e retorna `None`.

**Etapa 3: troca de linhas**

```python
if p != k:
    A[[k, p]] = A[[p, k]]
    b[[k, p]] = b[[p, k]]
```

Só troca se for necessário. `A[[k, p]] = A[[p, k]]` usa *indexação avançada* (lista de índices) do NumPy. O lado direito `A[[p, k]]` cria uma **cópia** das duas linhas já invertidas, e depois a atribuição as grava nas posições `k` e `p`. A troca é feita em todas as colunas, e o vetor `b` é trocado junto.

> **Cuidado:** a forma "intuitiva" `A[k], A[p] = A[p], A[k]` **não funciona** em NumPy. `A[k]` e `A[p]` são *views* (janelas para a mesma memória), e a primeira atribuição sobrescreve dados que a segunda ainda precisa ler, o que acaba duplicando uma linha e perdendo a outra. A indexação com lista evita esse problema.

**Etapa 4: eliminação progressiva**

```python
for i in range(k + 1, n):
    fator = A[i, k] / A[k, k]
    for j in range(k + 1, n):
        A[i, j] = A[i, j] - fator * A[k, j]
    b[i] = b[i] - fator * b[k]
    A[i, k] = 0.0
```

Para cada linha $i$ abaixo do pivô:

1. `fator` é o multiplicador $m_{ik}=a_{ik}/a_{kk}$. Com o pivoteamento, $\lvert \text{fator}\rvert\le 1$.
2. O laço em `j` subtrai `fator × linha k` da linha $i$ para as colunas à direita do pivô ($j>k$).
3. O mesmo é feito no vetor `b`. Esta linha precisa estar **fora** do laço de `j`, pois `b[i]` é atualizado uma única vez por linha.
4. `A[i, k] = 0.0` registra explicitamente o zero. Não se calcula `A[i,k] - fator*A[k,k]`, que daria zero a menos de arredondamento, e a coluna `j=k` é pulada para economizar uma operação por linha.

**Teste do último pivô**

```python
if abs(A[n-1, n-1]) <= limite_pivo:
    print("Erro: Pivô nulo na última linha. Execução interrompida.")
    return None
```

O laço externo vai só até $n-2$, então o último pivô $u_{nn}$ nunca foi testado. Ele é o denominador da primeira divisão da substituição regressiva, por isso é verificado aqui.

**Substituição regressiva**

```python
F_calc = np.zeros(n)
F_calc[n - 1] = b[n - 1] / A[n - 1, n - 1]

for i in range(n - 2, -1, -1):
    soma = b[i]
    for j in range(i + 1, n):
        soma = soma - A[i, j] * F_calc[j]
    F_calc[i] = soma / A[i, i]
```

- Resolve a última incógnita diretamente: $F_n = b_n/u_{nn}$.
- `range(n-2, -1, -1)` percorre $i=n-2,n-3,\dots,0$ (de trás para frente).
- Para cada $i$, parte de `b[i]`, subtrai as contribuições das incógnitas já conhecidas (`A[i,j]*F_calc[j]` para $j>i$) e divide pelo elemento da diagonal.

### 4.4 `rodar_cenario(nome, th1, th2)`

Função de conveniência para executar e avaliar um cenário completo.

```python
A_orig, b_orig, F_ref = criar_sistema_barras(th1, th2)
F_calc = eliminacao_gaussiana_pivoteamento(A_orig, b_orig)
```

Monta o sistema e resolve. O solucionador **não recebe** `F_ref`.

```python
erro = F_calc - F_ref
norma_inf_erro = np.max(np.abs(erro))

residuo = np.dot(A_orig, F_calc) - b_orig
norma_inf_residuo = np.max(np.abs(residuo))
```

Calcula as duas métricas de qualidade:

- **Erro:** $e=\hat F-F_{ref}$ e $\lVert e\rVert_\infty$.
- **Resíduo:** $r=A_{orig}\hat F-b_{orig}$ e $\lVert r\rVert_\infty$. Note o uso de `A_orig` e `b_orig`, pois `A` e `b` internos ao solucionador foram modificados.

```python
print(f"F_calc : {np.round(F_calc, 4)}")
print(f"... {norma_inf_erro:.4e}")
```

`np.round(..., 4)` só afeta a **exibição**; os valores usados nos cálculos de erro e resíduo continuam com a precisão completa. O formato `:.4e` imprime em notação científica com 4 casas.

A função retorna `A_orig, b_orig, F_ref` para reutilização (ou `None, None, None` se o solucionador falhar).

### 4.5 Bloco principal

```python
if __name__ == "__main__":
```

O código abaixo só roda quando o arquivo é executado diretamente (e não quando importado como módulo).

```python
A_A, b_A, F_A = rodar_cenario("A (Afastadas)", 15.0, 65.0)
A_C, b_C, F_C = rodar_cenario("C (Quase Paralelas)", 0.01, 0.02)
```

Executa os cenários A (barras afastadas) e C (quase paralelas).

**Teste de sensibilidade (condicionamento)**

```python
b_C_perturbado = b_C.copy()
b_C_perturbado[1] += 0.01

F_calc_C_pert = eliminacao_gaussiana_pivoteamento(A_C, b_C_perturbado)
delta_F = np.max(np.abs(F_calc_C_pert - F_C))
```

- Copia `b_C` (para não alterar o original) e soma 0,01 N à **segunda** componente.
- Resolve o sistema perturbado **com a mesma matriz**.
- Mede $\lVert\Delta F\rVert_\infty$, a maior mudança nas forças. Compare com $\lVert\Delta b\rVert_\infty = 0{,}01$.

Isso mede a sensibilidade do **problema**: a mesma perturbação feita em qualquer algoritmo exato daria o mesmo $\Delta F$.

### 4.6 Correspondência pseudocódigo → Python

| Pseudocódigo (enunciado) | Python (código) |
|---|---|
| `do k = 1,n-1` | `for k in range(n - 1)` |
| `p = k` | `p = k` |
| `if (abs(A(i,k)) > abs(A(p,k))) p = i` | `if abs(A[i, k]) > abs(A[p, k]): p = i` |
| `if (abs(A(p,k)) <= limite_pivo) parar` | `if abs(A[p, k]) <= limite_pivo: ... return None` |
| troca de linhas `do j=1,n ... end do` | `A[[k, p]] = A[[p, k]]` (vetorizado) |
| `fator = A(i,k)/A(k,k)` | `fator = A[i, k] / A[k, k]` |
| `A(i,j) = A(i,j) - fator*A(k,j)` | `A[i, j] = A[i, j] - fator * A[k, j]` |
| `A(i,k) = 0.0` | `A[i, k] = 0.0` |
| `do i = n-1,1,-1` | `for i in range(n - 2, -1, -1)` |

---

## 5. Resultados obtidos e interpretação

### 5.1 Saída do programa

```
--- Analisando Cenário A (Afastadas) (θ1 = 15.0°, θ2 = 65.0°) ---
F_ref  : [100. -80.  20.  35.  15.]
F_calc : [100. -80.  20.  35.  15.]
Norma infinito do erro (||ΔF||_inf): 7.1054e-15
Norma infinito do resíduo (||r||_inf): 1.4211e-14

--- Analisando Cenário C (Quase Paralelas) (θ1 = 0.01°, θ2 = 0.02°) ---
F_ref  : [100. -80.  20.  35.  15.]
F_calc : [100. -80.  20.  35.  15.]
Norma infinito do erro (||ΔF||_inf): 5.1514e-11
Norma infinito do resíduo (||r||_inf): 3.5527e-15

--- Teste de Sensibilidade (Perturbação +0.01 N na 2ª carga) ---
Variação em b (||Δb||_inf): 0.01
Variação nas forças calculadas no Cenário C (||ΔF||_inf): 2.0000e-02
```

### 5.2 Cargas e passo a passo da eliminação

Vetores $b = A F_{ref}$ gerados:

| Cenário | $b$ |
|---|---|
| A | (−26,6227; 62,7831; 4,7487; −39,7487; 15) |
| C | (19,9895; 20,0000; 4,7487; −39,7487; 15) |

**Etapa k = 1 (índice 0):** a coluna 1 tem apenas dois elementos não nulos, nas linhas 1 e 2.

| | Linha 1: $\lvert\sin\theta_1\rvert$ | Linha 2: $\lvert\cos\theta_1\rvert$ | Linha pivô |
|---|---|---|---|
| Cenário A | 0,259 | 0,966 | **2** (troca 1↔2) |
| Cenário C | 1,7·10⁻⁴ | ≈ 1 | **2** (troca 1↔2) |

Após a troca, o fator é $m=\tan\theta_1$, e a nova linha 2 fica:

- **A:** $[0,\ 0{,}7931,\ 1,\ 0,\ 0]$ com $b_2' = -43{,}445$
- **C:** $[0,\ 1{,}745\cdot10^{-4},\ 1,\ 0,\ 0]$ com $b_2' = 19{,}986$

**Etapas k = 2, 3, 4:** a matriz já está triangular superior depois da etapa 1 (as linhas 3, 4 e 5 já tinham zeros à esquerda da diagonal). Os fatores calculados são todos 0 e nada muda. Os pivôs da diagonal ficam:

$$
U_{\text{diag}} = \bigl(\cos\theta_1,\ \ \sin\theta_2-\tan\theta_1\cos\theta_2,\ \ -1,\ \ -\sin\phi,\ \ 1\bigr)
$$

**Substituição regressiva** (de baixo para cima):

1. $F_5 = b_5/1 = 15$
2. $F_4 = (b_4 + F_5)/(-\sin\phi)$, com $b_4=-39{,}7487$, dá 35
3. $F_3 = (b_3 - \cos\phi\,F_4)/(-1)$, dá 20
4. $F_2 = (b_2' - F_3)/(\sin\theta_2-\tan\theta_1\cos\theta_2)$, dá $-80$
5. $F_1 = (b_1 - \cos\theta_2 F_2)/\cos\theta_1$, dá 100

No cenário C o passo 4 divide por $1{,}745\cdot10^{-4}$, e é aí que o mau condicionamento "amplifica" os erros de arredondamento.

### 5.3 Teste de sensibilidade

| | Cenário A | Cenário C |
|---|---|---|
| κ∞(A) | 11,65 | 45 841 |
| erro ‖ΔF‖∞ (resolução original) | 7,1·10⁻¹⁵ | 5,2·10⁻¹¹ |
| resíduo ‖r‖∞ | 1,4·10⁻¹⁴ | 3,6·10⁻¹⁵ |

**Resultado do teste:** com $\Delta b_2 = +0{,}01$ N, o cenário C teve $\lVert\Delta F\rVert_\infty = 0{,}02$ N, um fator de apenas 2.

**Isso não significa que C é bem condicionado.** A razão é a **direção** da perturbação:

- $\Delta F_1 = -\dfrac{\sin\theta_2}{\det}\Delta b_2$, e em C, $\sin\theta_2$ é minúsculo, o que cancela quase todo o fator $1/\det$.
- Uma perturbação **vertical** em barras quase verticais é uma direção "benigna".
- A direção perigosa é a **horizontal**: perturbar $b_1$ (primeira componente) em 0,01 N daria $\Delta F \approx \pm 0{,}01/1{,}745\cdot10^{-4} \approx \pm 57{,}3$ N no cenário C, uma amplificação de ~5700 vezes. (Em A, a mesma perturbação dá apenas −0,0055 e +0,0126 N.)

O número de condição $\kappa$ é um limite para o **pior caso**. Para ver a sensibilidade de C, teste também `b_C_perturbado[0] += 0.01`.

### 5.4 O pivoteamento "não ajudou" em C?

Em precisão dupla (~16 dígitos) o pivô "pequeno" de C ($1{,}7\cdot10^{-4}$) ainda é enorme comparado ao épsilon de máquina ($\approx 2\cdot10^{-16}$). O fator $\tan 0{,}01° \approx 1{,}7\cdot10^{-4}$ é pequeno, não gigante, e o erro cometido com ou sem pivoteamento é da mesma ordem (~$10^{-11}$), consistente com $\kappa\cdot\varepsilon\cdot\lVert F\rVert\approx 5\cdot10^{4}\cdot10^{-16}\cdot100$.

A diferença aparece em precisão reduzida. Em testes com arredondamento a 4 algarismos significativos em cada operação, o cenário C teve erro ~20 N com resíduo ~5,4 (eliminação sem pivô) contra erro ~23 N com resíduo ~0,004 (com pivoteamento). O pivoteamento reduz o **resíduo** (estabilidade), mas o **erro** continua grande, porque é governado pelo condicionamento.

---

## 6. Limitações e melhorias

Este código resolve o núcleo da atividade (solucionador com pivoteamento, erro, resíduo e sensibilidade), mas há itens a observar e pontos pedidos no enunciado que ele ainda não cobre.

### 6.1 Pontos de atenção no código atual

1. **Comentários `[cite: 2]`.** Ficaram resíduos de marcação de citação nos comentários e docstrings (por exemplo `...fornecidas[cite: 2]`). Remova-os antes de publicar.

2. **Limiar absoluto fixo.** `limite_pivo=1e-12` é absoluto. O enunciado alerta que um limiar fixo não serve para todas as unidades e escalas (se a matriz for multiplicada por $10^{-9}$, todos os pivôs parecerão nulos). Alternativa relativa:

   ```python
   limite = limite_pivo * n * np.max(np.abs(A_orig))
   ```

   e usar `limite` nas duas comparações. Registre o critério adotado no README/relatório.

3. **Tratamento de `None` no teste de sensibilidade.** Se a eliminação falhar, `F_calc_C_pert` será `None` e `F_calc_C_pert - F_C` lançará `TypeError`. Convém testar:

   ```python
   if F_calc_C_pert is not None:
       ...
   ```

4. **Referência do ΔF.** O código calcula `F_calc_C_pert - F_C` (comparando com $F_{ref}$). O mais rigoroso é comparar com a solução calculada **sem** perturbação, para que $\Delta F$ meça só o efeito da perturbação. Aqui a diferença é desprezível (~$10^{-11}$), mas é o procedimento correto. O mesmo teste não foi feito no cenário A.

5. **Variáveis não usadas.** `A_A, b_A, F_A` não são reaproveitadas depois. Não é erro, mas pode ser simplificado.

### 6.2 Itens do enunciado ainda não implementados

| Item da atividade | Situação |
|---|---|
| Imprimir, por etapa, k, linha p, pivô e matriz | Não implementado |
| Registrar as trocas de linhas realizadas | Não implementado |
| Comparar com eliminação **ingênua** (sem pivô) | Não implementado |
| Cenário C em precisão reduzida (4 algarismos significativos) | Não implementado |
| Perturbação em A e em C com comparação ∥Δb∥ × ∥ΔF∥ | Só no cenário C |
| Contadores de FLOPs (mult/div e soma/sub) | Não implementado |
| Cálculo do número de condição | Não implementado |
| Parte 1 (sistema 2×2, determinantes das configurações A, B, C) | Não implementado |

### 6.3 Sugestões de extensão

**Eliminação ingênua** (para comparar): acrescentar um parâmetro `pivotear=True` e pular a busca/troca quando for `False`:

```python
def eliminacao_gaussiana(A_orig, b_orig, pivotear=True, limite_pivo=1e-12):
    ...
    p = k
    if pivotear:
        for i in range(k + 1, n):
            if abs(A[i, k]) > abs(A[p, k]):
                p = i
    ...
```

**Registro das trocas e impressão por etapa:**

```python
trocas = []
...
if p != k:
    A[[k, p]] = A[[p, k]]
    b[[k, p]] = b[[p, k]]
    trocas.append((k + 1, p + 1))     # +1 para numerar a partir de 1
print(f"k={k+1}, p={p+1}, pivô={A[k, k]:.6g}")
print(A)
```

**Número de condição:**

```python
kappa = np.linalg.cond(A_orig, np.inf)
```

**Contadores de FLOPs:** crie `mult_div = 0` e `soma_sub = 0`. Incremente:

- `mult_div += 1` em `fator = ...` (divisão) e em cada `fator * A[k, j]`, `fator * b[k]` e `A[i, j] * F_calc[j]`, mais as divisões da substituição regressiva;
- `soma_sub += 1` em cada subtração de `A[i, j]`, `b[i]` e `soma`.

Não conte comparações, trocas de linha nem operações de controle. Para $n=5$ o esperado é **65** multiplicações/divisões e **50** somas/subtrações (a mesma contagem da versão ingênua).

**Precisão reduzida:** a simulação deve arredondar o resultado de **cada operação**, e não apenas a impressão. Exemplo de função auxiliar:

```python
from math import floor, log10

def arred_sig(x, sig=4):
    if x == 0:
        return 0.0
    return round(x, sig - 1 - int(floor(log10(abs(x)))))
```

aplicada em torno de cada multiplicação, divisão, soma e subtração do algoritmo.

---

## 7. Glossário

| Termo | Significado |
|---|---|
| **Pivô** | Elemento da diagonal $a_{kk}$ usado para eliminar a coluna $k$ |
| **Fator de eliminação (multiplicador)** | $m_{ik}=a_{ik}/a_{kk}$ |
| **Pivoteamento parcial** | Trocar linhas para que o pivô seja o maior valor absoluto da coluna |
| **Matriz triangular superior** | Matriz com zeros abaixo da diagonal |
| **Substituição regressiva** | Resolver um sistema triangular superior de baixo para cima |
| **Singular** | Matriz sem inversa ($\det=0$) |
| **Número de condição κ(A)** | $\lVert A\rVert\lVert A^{-1}\rVert$; mede a sensibilidade do problema |
| **Norma infinito** | Maior valor absoluto das componentes de um vetor |
| **Resíduo** | $r=A\hat F-b$ |
| **Épsilon de máquina** | Menor erro relativo de arredondamento; ≈ $2{,}2\cdot10^{-16}$ em dupla precisão |
| **Backward stable** | Algoritmo que dá a solução exata de um problema levemente perturbado |
| **Treliça** | Estrutura de barras articuladas que só transmitem esforço axial |

---

## 8. Referências

- Burden, R. L.; Faires, J. D. *Numerical Analysis*. Cengage.
- Chapra, S. C.; Canale, R. P. *Métodos Numéricos para Engenharia*. McGraw-Hill.
- Trefethen, L. N.; Bau, D. *Numerical Linear Algebra*. SIAM.
- Golub, G. H.; Van Loan, C. F. *Matrix Computations*. Johns Hopkins University Press.

---

## Autor

*Seu nome* · *Curso / Instituição* · *Ano*
