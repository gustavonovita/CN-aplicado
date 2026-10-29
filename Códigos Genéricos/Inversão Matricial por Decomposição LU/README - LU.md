# Inversão Matricial por Decomposição LU

Implementação em Python do cálculo da matriz inversa por meio da decomposição LU. O código contém duas versões:

1. **Decomposição LU sem pivoteamento** — constrói explicitamente as matrizes triangulares inferior \(L\) e superior \(U\).
2. **Decomposição LU com pivoteamento parcial** — usa trocas de linhas e armazena os fatores \(L\) e \(U\) em uma única matriz, além de um vetor de permutação.

O objetivo deste README é servir como documentação de consulta: apresenta a fundamentação matemática, explica as etapas do algoritmo e descreve como executar e testar as funções.

---

## Sumário

- [1. Objetivo](#1-objetivo)
- [2. Conceitos matemáticos](#2-conceitos-matemáticos)
- [3. Como obter a inversa usando LU](#3-como-obter-a-inversa-usando-lu)
- [4. Versão 1: LU sem pivoteamento](#4-versão-1-lu-sem-pivoteamento)
- [5. Versão 2: LU com pivoteamento parcial](#5-versão-2-lu-com-pivoteamento-parcial)
- [6. Tolerância numérica e limitações](#6-tolerância-numérica-e-limitações)
- [7. Como usar o código](#7-como-usar-o-código)
- [8. Exemplo completo de teste](#8-exemplo-completo-de-teste)
- [9. Como verificar o resultado](#9-como-verificar-o-resultado)
- [10. Complexidade computacional](#10-complexidade-computacional)
- [11. Diferenças entre as implementações](#11-diferenças-entre-as-implementações)
- [12. Possíveis melhorias](#12-possíveis-melhorias)
- [13. Referências](#13-referências)

---

## 1. Objetivo

Dada uma matriz quadrada \(A\), queremos calcular sua inversa \(A^{-1}\), isto é, uma matriz que satisfaça

\[
AA^{-1}=A^{-1}A=I,
\]

em que \(I\) é a matriz identidade.

A matriz inversa só existe se \(A\) for quadrada e não singular, ou seja, se seu determinante for diferente de zero. Em termos computacionais, uma matriz singular ou numericamente muito próxima de singular pode causar divisões por zero, amplificação de erros de arredondamento ou resultados pouco confiáveis.

Em vez de calcular a inversa diretamente por determinantes e cofatores, este programa usa a decomposição LU e resolve uma sequência de sistemas lineares. Essa abordagem é especialmente útil para estudar métodos numéricos de álgebra linear.

## 2. Conceitos matemáticos

### 2.1 Matriz triangular inferior \(L\)

Uma matriz triangular inferior tem zeros acima da diagonal principal. Na versão sem pivoteamento, o algoritmo constrói \(L\) com diagonal unitária:

\[
L=
\begin{bmatrix}
1 & 0 & 0\\
\ell_{21} & 1 & 0\\
\ell_{31} & \ell_{32} & 1
\end{bmatrix}.
\]

Os elementos \(\ell_{ik}\) guardam os fatores usados para eliminar os elementos abaixo dos pivôs durante a eliminação de Gauss.

### 2.2 Matriz triangular superior \(U\)

Uma matriz triangular superior tem zeros abaixo da diagonal principal:

\[
U=
\begin{bmatrix}
u_{11} & u_{12} & u_{13}\\
0 & u_{22} & u_{23}\\
0 & 0 & u_{33}
\end{bmatrix}.
\]

Ela é obtida pela eliminação dos elementos abaixo da diagonal de \(A\). Os elementos da diagonal de \(U\) são os pivôs; para usar substituição retroativa, eles não podem ser zero.

### 2.3 Decomposição LU

Na versão sem pivoteamento, procura-se fatorar a matriz como

\[
A=LU.
\]

Na versão com pivoteamento parcial, podem ocorrer trocas de linhas. O código registra essas trocas no vetor `P` e trabalha com a relação

\[
PA=LU,
\]

em que \(P\) representa uma matriz de permutação. No código, `P` é armazenado como uma lista de índices, e não como uma matriz bidimensional.

### 2.4 Eliminação de Gauss e fator de eliminação

Em uma etapa \(k\), para zerar a entrada \(a_{ik}\), abaixo do pivô \(a_{kk}\), calcula-se o fator

\[
m_{ik}=\frac{a_{ik}}{a_{kk}}.
\]

Depois, atualizam-se as entradas restantes da linha:

\[
a_{ij}\leftarrow a_{ij}-m_{ik}a_{kj}.
\]

Na versão sem pivoteamento, esse fator é salvo em `L[i][k]`. Na versão com pivoteamento, ele é guardado diretamente na posição `LU[i][k]`, enquanto a parte triangular superior fica armazenada na mesma matriz.

### 2.5 Vetores da base e matriz identidade

Se \(e_c\) representa a coluna \(c\) da matriz identidade, então

\[
e_0=\begin{bmatrix}1\\0\\\vdots\\0\end{bmatrix},\quad
e_1=\begin{bmatrix}0\\1\\\vdots\\0\end{bmatrix},\quad \ldots
\]

Usar índices começando em zero é a convenção do Python. Matematicamente, essas colunas costumam ser numeradas de 1 até \(n\).

A identidade é formada por esses vetores:

\[
I=[e_0\ e_1\ \cdots\ e_{n-1}].
\]

## 3. Como obter a inversa usando LU

A definição da inversa é

\[
AA^{-1}=I.
\]

Escrevendo a inversa por colunas,

\[
A^{-1}=[x_0\ x_1\ \cdots\ x_{n-1}],
\]

cada coluna \(x_c\) precisa satisfazer

\[
Ax_c=e_c.
\]

Portanto, em vez de calcular a inversa de uma só vez, o algoritmo resolve \(n\) sistemas lineares, um para cada coluna da identidade.

### 3.1 Sem pivoteamento

Como \(A=LU\),

\[
LUx_c=e_c.
\]

Define-se um vetor intermediário \(d\) tal que

\[
Ld=e_c,\qquad Ux_c=d.
\]

A primeira equação é resolvida por **substituição progressiva**; a segunda, por **substituição retroativa**.

### 3.2 Com pivoteamento parcial

Como \(PA=LU\), multiplicando \(Ax_c=e_c\) por \(P\),

\[
PAx_c=Pe_c
\]

e, portanto,

\[
LUx_c=Pe_c.
\]

Assim, resolvemos

\[
Ld=Pe_c,\qquad Ux_c=d.
\]

O vetor de permutação é importante: após as trocas de linhas, o lado direito também precisa respeitar a mesma permutação. No código, a componente \(i\) de \(Pe_c\) vale 1 quando `P[i] == c`, e vale 0 nos demais casos.

### 3.3 Substituição progressiva

Para um sistema \(Ld=b\), com \(L\) triangular inferior e diagonal unitária,

\[
d_i=b_i-\sum_{j=0}^{i-1}L_{ij}d_j.
\]

O cálculo começa na primeira linha e avança para baixo, pois cada novo valor depende dos anteriores.

Na versão sem pivoteamento, \(b=e_c\), logo o termo inicial é 1 se \(i=c\), e 0 caso contrário. Na versão com pivoteamento, \(b=Pe_c\).

### 3.4 Substituição retroativa

Para um sistema \(Ux_c=d\),

\[
x_{i,c}=\frac{d_i-\sum_{j=i+1}^{n-1}U_{ij}x_{j,c}}{U_{ii}}.
\]

O cálculo começa na última linha e sobe, porque cada valor depende dos valores já calculados abaixo dele.

O resultado \(x_c\) é armazenado como a coluna `c` de `A_inv`. Quando todas as colunas forem calculadas, `A_inv` será a matriz inversa.

## 4. Versão 1: LU sem pivoteamento

Função: `inversao_lu(A, tol=1e-15)`

### 4.1 Preparação

```python
n = len(A)
L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
U = [[float(val) for val in linha] for linha in A]
```

- `n = len(A)` obtém o número de linhas. O algoritmo pressupõe que a entrada seja uma matriz quadrada.
- `L` começa como a identidade: diagonal igual a 1 e demais entradas iguais a 0.
- `U` recebe uma cópia dos valores de `A`, convertidos para `float`. Assim, as operações não alteram diretamente a lista original fornecida pelo usuário.

### 4.2 Fatoração

```python
for k in range(n - 1):
    if abs(U[k][k]) < tol:
        raise ZeroDivisionError(...)
    for i in range(k + 1, n):
        fator = U[i][k] / U[k][k]
        L[i][k] = fator
        U[i][k] = 0.0
        for j in range(k + 1, n):
            U[i][j] -= fator * U[k][j]
```

- `k` identifica a coluna e a linha do pivô atual.
- A verificação de `U[k][k]` evita dividir por um pivô considerado nulo ou muito pequeno.
- `i` percorre as linhas abaixo do pivô.
- `fator` calcula quanto da linha do pivô deve ser subtraído para zerar a entrada abaixo dele.
- `L[i][k] = fator` registra esse multiplicador.
- `U[i][k] = 0.0` registra o zero produzido pela eliminação.
- `j` percorre as colunas à direita do pivô e atualiza a linha.

Após as etapas, `L` e `U` representam a fatoração \(A=LU\), desde que os pivôs necessários sejam válidos.

### 4.3 Verificação do último pivô

```python
if abs(U[n - 1][n - 1]) < tol:
    raise ZeroDivisionError(...)
```

O laço principal vai até `n - 2`, então o último pivô precisa ser verificado separadamente. Se ele for muito pequeno, a matriz é tratada como singular ou numericamente problemática.

### 4.4 Cálculo das colunas da inversa

```python
A_inv = [[0.0] * n for _ in range(n)]
for c in range(n):
    d = [0.0] * n
    for i in range(n):
        soma = 1.0 if i == c else 0.0
        for j in range(i):
            soma -= L[i][j] * d[j]
        d[i] = soma
```

- `A_inv` começa como uma matriz \(n\times n\) de zeros.
- `c` escolhe a coluna da identidade a ser usada no sistema atual.
- `d` é o vetor intermediário da substituição progressiva.
- `soma` começa com a componente correspondente de \(e_c\).
- O laço sobre `j` subtrai as contribuições já conhecidas de \(L\).
- `d[i] = soma` finaliza a componente atual de \(d\).

Em seguida:

```python
for i in range(n - 1, -1, -1):
    soma = d[i]
    for j in range(i + 1, n):
        soma -= U[i][j] * A_inv[j][c]
    A_inv[i][c] = soma / U[i][i]
```

- `range(n - 1, -1, -1)` percorre as linhas de baixo para cima.
- `soma` começa em `d[i]`.
- O laço sobre `j` remove as contribuições das incógnitas já calculadas.
- A divisão por `U[i][i]` isola a incógnita da linha.
- `A_inv[i][c]` guarda o elemento da linha `i` e coluna `c` da inversa.

Ao concluir as \(n\) colunas, `return A_inv` devolve a matriz inversa.

### 4.5 Limitação central

Esta versão **não troca linhas**. Se um pivô for zero ou muito pequeno, a função lança uma exceção mesmo que a matriz original seja inversível. Por exemplo, a matriz

\[
A=\begin{bmatrix}0&1\\1&0\end{bmatrix}
\]

é inversível, mas o primeiro pivô é zero. A versão com pivoteamento parcial é mais adequada para esse caso.

## 5. Versão 2: LU com pivoteamento parcial

Função: `inversao_lu_pivoteamento(A, tol=1e-15)`

Esta versão troca linhas para escolher, na coluna atual, o elemento de maior valor absoluto entre a linha do pivô e as linhas abaixo dela. Essa estratégia é chamada **pivoteamento parcial**.

### 5.1 Preparação e permutação

```python
n = len(A)
LU = [[float(val) for val in linha] for linha in A]
P = list(range(n))
```

- `LU` é uma cópia da matriz original.
- `P` começa como `[0, 1, 2, ..., n-1]`. Cada posição registra qual linha original está naquela posição após as permutações.
- Em vez de criar matrizes `L` e `U` separadas, o programa armazena os multiplicadores de `L` abaixo da diagonal e os elementos de `U` na diagonal e acima dela. É uma **fatoração LU compacta**.

### 5.2 Escolha do pivô e troca de linhas

```python
p = max(range(k, n), key=lambda i: abs(LU[i][k]))
if abs(LU[p][k]) < tol:
    raise ZeroDivisionError(...)
if p != k:
    LU[k], LU[p] = LU[p], LU[k]
    P[k], P[p] = P[p], P[k]
```

- `max(...)` seleciona o índice `p` da linha com maior valor absoluto na coluna `k`, entre as linhas ainda disponíveis.
- `lambda i: abs(LU[i][k])` é uma função curta usada como critério de comparação.
- Se o maior candidato for menor que `tol`, o pivô é considerado inválido.
- Se `p != k`, as linhas `k` e `p` são trocadas tanto em `LU` quanto em `P`. Atualizar ambos é essencial para manter o registro da permutação.

### 5.3 Eliminação compacta

```python
for i in range(k + 1, n):
    LU[i][k] /= LU[k][k]
    for j in range(k + 1, n):
        LU[i][j] -= LU[i][k] * LU[k][j]
```

- A divisão armazena o multiplicador de eliminação em `LU[i][k]`.
- O laço interno atualiza as entradas à direita do pivô.
- Diferentemente da primeira versão, a função não substitui a entrada abaixo do pivô por zero: ela guarda ali o multiplicador de \(L\). A estrutura resultante contém as informações necessárias para reconstruir \(L\) e \(U\).

A diagonal de \(L\) é considerada unitária, embora não esteja explicitamente armazenada.

### 5.4 Resolução dos sistemas

Para cada coluna `c`, a função resolve

\[
Ld=Pe_c,\qquad Ux_c=d.
\]

Na substituição progressiva:

```python
soma = 1.0 if P[i] == c else 0.0
for j in range(i):
    soma -= LU[i][j] * d[j]
d[i] = soma
```

A condição `P[i] == c` representa a componente \(i\) do vetor permutado \(Pe_c\). Os valores abaixo da diagonal em `LU` são os multiplicadores de \(L\), então são usados na soma.

A substituição retroativa é:

```python
soma = d[i]
for j in range(i + 1, n):
    soma -= LU[i][j] * A_inv[j][c]
A_inv[i][c] = soma / LU[i][i]
```

Nesta etapa, as entradas na diagonal e acima dela são os elementos de \(U\). O procedimento é equivalente ao da primeira função.

### 5.5 Por que usar pivoteamento?

O pivoteamento parcial:
- permite continuar quando o pivô atual é zero, mas existe um candidato não nulo abaixo;
- costuma melhorar a estabilidade numérica em relação à eliminação sem trocas;
- não garante, por si só, precisão para toda matriz mal condicionada;
- não transforma uma matriz singular em inversível.

A tolerância absoluta fixa também não é um teste universal de singularidade. Para aplicações numéricas mais exigentes, é importante considerar a escala dos dados e o condicionamento da matriz.

## 6. Tolerância numérica e limitações

As funções recebem `tol=1e-15`. Em Python, `1e-15` significa \(10^{-15}\).

O teste

```python
abs(pivo) < tol
```

compara o valor absoluto do pivô com essa tolerância. Se for menor, o código lança `ZeroDivisionError` para interromper o cálculo antes de uma divisão problemática.

**Atenção:** `tol` não é uma prova matemática de que a matriz é singular. Uma matriz pode ter pivôs pequenos e ainda ser inversível; por outro lado, pivôs acima da tolerância não garantem que a inversa seja numericamente precisa. O valor adequado depende da escala dos elementos e da precisão desejada.

Limitações da implementação atual:

1. Não há validação explícita de que a entrada seja quadrada.
2. Não há validação explícita de que todas as linhas tenham o mesmo comprimento.
3. Não há validação de que os valores sejam números finitos.
4. A tolerância é absoluta, não relativa à escala da matriz.
5. O código usa aritmética de ponto flutuante (`float`), sujeita a arredondamento.
6. A versão sem pivoteamento pode falhar para uma matriz inversível cujo pivô seja zero ou muito pequeno.
7. Para matrizes mal condicionadas, a inversa pode ser muito sensível a pequenas alterações nos dados.

O código também pressupõe uma matriz não vazia: em ambas as funções, a verificação do último pivô usa o índice `n - 1`. Uma matriz vazia causaria erro de índice.

## 7. Como usar o código

### Requisitos

- Python 3.
- Nenhuma biblioteca externa é necessária: as funções usam apenas recursos nativos do Python.

Salve o código em um arquivo, por exemplo `inversao_lu.py`. Se o arquivo tiver outro nome, use esse nome no comando de execução.

### Importar a função em outro script

```python
from inversao_lu import inversao_lu, inversao_lu_pivoteamento

A = [
    [4, 3],
    [6, 3]
]

inversa_1 = inversao_lu(A)
inversa_2 = inversao_lu_pivoteamento(A)

print("Inversa sem pivoteamento:")
for linha in inversa_1:
    print(linha)

print("Inversa com pivoteamento:")
for linha in inversa_2:
    print(linha)
```

O arquivo que contém as funções deve estar na mesma pasta do script que as importa, ou em um diretório disponível no caminho de importação do Python.

### Observação sobre execução direta

O arquivo original define as funções, mas não contém um bloco de exemplo `if __name__ == "__main__":`. Portanto, executá-lo diretamente não imprime um resultado: ele apenas define as funções. Para ver resultados, importe e chame uma delas, ou acrescente um bloco de teste.

## 8. Exemplo completo de teste

O exemplo abaixo pode ser colocado em um arquivo separado, como `teste_inversao.py`:

```python
from inversao_lu import inversao_lu, inversao_lu_pivoteamento

A = [
    [4, 3],
    [6, 3]
]

for nome, funcao in [
    ("LU sem pivoteamento", inversao_lu),
    ("LU com pivoteamento parcial", inversao_lu_pivoteamento),
]:
    print(f"\n{nome}:")
    inversa = funcao(A)
    for linha in inversa:
        print(["{:.6f}".format(valor) for valor in linha])
```

Para essa matriz,

\[
A=\begin{bmatrix}4&3\\6&3\end{bmatrix},
\qquad \det(A)=4\cdot3-3\cdot6=-6\ne0.
\]

Logo, ela é inversível. A inversa esperada é

\[
A^{-1}=
\begin{bmatrix}
-\frac12 & \frac12\\
1 & -\frac23
\end{bmatrix}.
\]

Em forma decimal:

```text
[-0.500000, 0.500000]
[1.000000, -0.666667]
```

### Teste que exige troca de linhas

Para demonstrar a vantagem do pivoteamento, teste também:

```python
B = [
    [0, 1],
    [1, 0]
]

print(inversao_lu_pivoteamento(B))
```

A inversa de `B` é a própria matriz `B`. A versão sem pivoteamento interrompe o cálculo porque o primeiro pivô é zero; a versão com pivoteamento troca as linhas e consegue prosseguir.

## 9. Como verificar o resultado

A verificação mais direta é multiplicar a matriz original pela inversa calculada. O resultado deve ser aproximadamente a identidade:

\[
AA^{-1}\approx I.
\]

Como as operações usam `float`, não se deve exigir igualdade exata em todos os elementos. Pequenas diferenças, como \(10^{-15}\), podem surgir por arredondamento.

Exemplo sem bibliotecas externas:

```python
def multiplicar_matrizes(A, B):
    linhas_A = len(A)
    colunas_A = len(A[0])
    colunas_B = len(B[0])

    resultado = [
        [0.0 for _ in range(colunas_B)]
        for _ in range(linhas_A)
    ]

    for i in range(linhas_A):
        for j in range(colunas_B):
            for k in range(colunas_A):
                resultado[i][j] += A[i][k] * B[k][j]

    return resultado


def imprimir_matriz(M):
    for linha in M:
        print(["{:.10f}".format(valor) for valor in linha])


A = [
    [4, 3],
    [6, 3]
]

A_inv = inversao_lu_pivoteamento(A)
produto = multiplicar_matrizes(A, A_inv)

print("A * A_inv:")
imprimir_matriz(produto)
```

O resultado deve estar muito próximo de

\[
I=\begin{bmatrix}1&0\\0&1\end{bmatrix}.
\]

Uma medida útil do erro é o maior desvio absoluto entre os elementos do produto e os elementos correspondentes da identidade:

\[
r_{\max}=\max_{i,j}|(AA^{-1})_{ij}-I_{ij}|.
\]

Quanto menor esse resíduo, melhor a identidade foi satisfeita numericamente. Um resíduo pequeno é uma verificação útil, mas não substitui uma análise do condicionamento da matriz.

## 10. Complexidade computacional

Para uma matriz \(n\times n\):

- A fatoração LU exige, em geral, \(O(n^3)\) operações.
- Cada substituição progressiva ou retroativa exige \(O(n^2)\) operações.
- Como são resolvidos \(n\) sistemas, as substituições juntas exigem \(O(n^3)\) operações.
- O armazenamento das matrizes e vetores exige \(O(n^2)\) memória.

Assim, o custo total para calcular a inversa é da ordem de \(O(n^3)\) operações e \(O(n^2)\) memória.

Se o objetivo for resolver apenas um sistema \(Ax=b\), normalmente não é necessário calcular a inversa completa. É mais eficiente fatorar \(A\) e resolver o sistema diretamente. Bibliotecas numéricas consolidadas, como NumPy e SciPy, também oferecem rotinas altamente otimizadas e testadas.

## 11. Diferenças entre as implementações

| Característica | `inversao_lu` | `inversao_lu_pivoteamento` |
|---|---|---|
| Fatoração | \(A=LU\) | \(PA=LU\) |
| Troca de linhas | Não | Sim, pivoteamento parcial |
| Armazenamento | `L` e `U` separados | Fatores compactados em `LU` |
| Vetor de permutação | Não | Sim, `P` |
| Tolerância a pivôs nulos | Baixa | Melhor quando há candidato não nulo abaixo |
| Custo assintótico | \(O(n^3)\) | \(O(n^3)\) |

**Recomendação prática:** para experimentação didática, comparar as duas versões ajuda a entender o papel do pivoteamento. Para matrizes gerais, prefira a versão com pivoteamento, sem esquecer que matrizes mal condicionadas exigem cuidados adicionais.

## 12. Possíveis melhorias

Se este projeto for evoluído, considere:

- validar a dimensão da matriz e rejeitar entradas vazias ou não quadradas com mensagens claras;
- validar linhas com comprimentos iguais e valores numéricos finitos;
- usar uma tolerância relativa baseada na escala da matriz;
- acrescentar testes automatizados para matrizes pequenas, matrizes que exigem pivoteamento e matrizes singulares;
- verificar o resíduo \(AA^{-1}-I\);
- comparar os resultados com `numpy.linalg.inv` apenas como referência de teste;
- separar a fatoração LU da etapa de resolução, permitindo reutilizar a fatoração para vários vetores \(b\);
- documentar as funções com *docstrings* e adicionar exemplos de uso;
- em aplicações reais, resolver sistemas lineares diretamente em vez de formar a inversa quando ela não for necessária.

## 13. Referências

Para aprofundar os assuntos envolvidos, consulte materiais sobre:

- Álgebra linear: matriz inversa, matriz identidade e sistemas lineares.
- Métodos numéricos: eliminação de Gauss, decomposição LU, substituição progressiva e retroativa.
- Pivoteamento parcial e estabilidade numérica.
- Documentação oficial do Python: https://docs.python.org/3/
- Documentação oficial do NumPy, útil para comparação numérica: https://numpy.org/doc/

---

## Licença

Nenhuma licença específica é definida por este README. Se você pretende tornar o repositório público e permitir reutilização do código, escolha e adicione uma licença apropriada (por exemplo, MIT), respeitando as condições que desejar estabelecer.
