# Método da Eliminação Gaussiana com Pivoteamento Parcial

Implementação em Python para resolver sistemas lineares $Ax=b$ pelo método da eliminação gaussiana com pivoteamento parcial. O algoritmo transforma a matriz dos coeficientes em uma matriz triangular superior e, em seguida, calcula as incógnitas por substituição retroativa.

Este README documenta a teoria matemática, explica o código linha por linha, apresenta exemplos de uso, descreve as limitações numéricas e mostra como verificar a solução.

---

## Sumário

- [1. Objetivo do programa](#1-objetivo-do-programa)
- [2. Formulação matemática](#2-formulação-matemática)
- [3. Eliminação gaussiana](#3-eliminação-gaussiana)
- [4. Pivoteamento parcial](#4-pivoteamento-parcial)
- [5. Substituição retroativa](#5-substituição-retroativa)
- [6. Explicação do código](#6-explicação-do-código)
- [7. Como executar](#7-como-executar)
- [8. Exemplo completo](#8-exemplo-completo)
- [9. Verificação da solução](#9-verificação-da-solução)
- [10. Tolerância e limitações](#10-tolerância-e-limitações)
- [11. Complexidade computacional](#11-complexidade-computacional)
- [12. Possíveis melhorias](#12-possíveis-melhorias)
- [13. Referências](#13-referências)

---

## 1. Objetivo do programa

O programa resolve um sistema de $n$ equações lineares com $n$ incógnitas:

$$
Ax=b,
$$

em que:

- $A$ é a matriz quadrada dos coeficientes;
- $x$ é o vetor das incógnitas;
- $b$ é o vetor dos termos independentes.

A função implementada é:

```python
eliminacao_gaussiana_pivoteamento_parcial(A, b, tol=1e-15)
```

Ela recebe a matriz $A$, o vetor $b$ e uma tolerância numérica opcional. Se o algoritmo conseguir resolver o sistema, retorna uma lista com os valores das incógnitas, na ordem $x_0,x_1,\ldots,x_{n-1}$.

O método não calcula a matriz inversa. Seu objetivo é encontrar diretamente a solução de um sistema linear.

## 2. Formulação matemática

Um sistema linear pode ser escrito como

$$
\begin{cases}
a_{11}x_1+a_{12}x_2+\cdots+a_{1n}x_n=b_1,\\
a_{21}x_1+a_{22}x_2+\cdots+a_{2n}x_n=b_2,\\
\vdots\\
a_{n1}x_1+a_{n2}x_2+\cdots+a_{nn}x_n=b_n.
\end{cases}
$$

Na notação matricial, isso é $Ax=b$.

Por exemplo,

$$
\begin{cases}
2x+y=5,\\
x-y=1
\end{cases}
$$

corresponde a

$$
A=\begin{bmatrix}2&1\\1&-1\end{bmatrix},\qquad
x=\begin{bmatrix}x\\y\end{bmatrix},\qquad
b=\begin{bmatrix}5\\1\end{bmatrix}.
$$

O método da eliminação gaussiana transforma o sistema em outro equivalente, mais fácil de resolver. As operações elementares usadas nas linhas são:

1. trocar duas linhas;
2. multiplicar uma linha por um número não nulo;
3. somar a uma linha um múltiplo de outra linha.

Essas operações preservam o conjunto de soluções quando realizadas corretamente sobre a matriz dos coeficientes e o vetor dos termos independentes.

## 3. Eliminação gaussiana

A ideia é transformar $A$ em uma matriz triangular superior $U$, zerando os elementos abaixo da diagonal principal.

Para uma matriz $3\times3$, a forma desejada é

$$
U=\begin{bmatrix}
u_{11}&u_{12}&u_{13}\\
0&u_{22}&u_{23}\\
0&0&u_{33}
\end{bmatrix}.
$$

Em cada etapa $k$, o elemento $M_{kk}$ é o pivô. Para eliminar o elemento da linha $i$, coluna $k$, calcula-se

$$
m_{ik}=\frac{M_{ik}}{M_{kk}}.
$$

Em seguida, para cada coluna à direita do pivô, atualiza-se

$$
M_{ij}\leftarrow M_{ij}-m_{ik}M_{kj},
$$

e também se atualiza o termo independente:

$$
v_i\leftarrow v_i-m_{ik}v_k.
$$

Ao terminar todas as etapas, o sistema transformado tem a forma triangular superior $Ux=\tilde b$, onde `M` contém $U$ e `v` contém o vetor de termos independentes transformado.

## 4. Pivoteamento parcial

### 4.1 Por que ele é necessário?

Se um pivô for zero, a divisão necessária para calcular o fator de eliminação não será possível. Mesmo quando o pivô não é exatamente zero, um valor muito pequeno pode amplificar erros de arredondamento.

O pivoteamento parcial escolhe, na coluna atual, o elemento de maior valor absoluto entre a linha do pivô e todas as linhas abaixo dela. A linha escolhida é trocada com a linha do pivô antes da eliminação.

Matematicamente, para a coluna $k$, escolhe-se o índice $p$ que maximiza

$$
|M_{ik}|,\qquad i=k,k+1,\ldots,n-1.
$$

O uso do valor absoluto significa que se compara o módulo dos elementos, não apenas seu valor algébrico.

### 4.2 Exemplo simples

Considere

$$
A=\begin{bmatrix}0&1\\1&0\end{bmatrix},\qquad
b=\begin{bmatrix}2\\3\end{bmatrix}.
$$

O primeiro elemento da diagonal é zero, então não pode ser usado como pivô. O pivoteamento parcial identifica o valor 1 na segunda linha e troca as linhas. O sistema passa a ser

$$
\begin{bmatrix}1&0\\0&1\end{bmatrix}
\begin{bmatrix}x_1\\x_2\end{bmatrix}
=
\begin{bmatrix}3\\2\end{bmatrix},
$$

resultando em $x_1=3$ e $x_2=2$.

### 4.3 Troca do vetor de termos independentes

Quando duas linhas de $A$ são trocadas, as posições correspondentes do vetor $b$ também precisam ser trocadas. Caso contrário, o programa estaria resolvendo um sistema diferente.

No código, as duas trocas são feitas juntas:

```python
M[k], M[p] = M[p], M[k]
v[k], v[p] = v[p], v[k]
```

## 5. Substituição retroativa

Depois da eliminação, o sistema tem forma triangular superior. Para três incógnitas, por exemplo,

$$
\begin{cases}
u_{11}x_1+u_{12}x_2+u_{13}x_3=\tilde b_1,\\
u_{22}x_2+u_{23}x_3=\tilde b_2,\\
u_{33}x_3=\tilde b_3.
\end{cases}
$$

Começa-se pela última equação:

$$
x_{n-1}=\frac{\tilde b_{n-1}}{u_{n-1,n-1}}.
$$

Depois, para $i=n-2,n-3,\ldots,0$,

$$
x_i=\frac{v_i-\displaystyle\sum_{j=i+1}^{n-1}M_{ij}x_j}{M_{ii}}.
$$

O algoritmo calcula primeiro a última incógnita, depois a penúltima e assim por diante até a primeira.

## 6. Explicação do código

### 6.1 Cabeçalho e cópias dos dados

```python
def eliminacao_gaussiana_pivoteamento_parcial(A, b, tol=1e-15):
    n = len(b)
    M = [[float(val) for val in linha] for linha in A]
    v = [float(val) for val in b]
```

- `def` define uma função reutilizável.
- `A` é a matriz dos coeficientes do sistema.
- `b` é o vetor dos termos independentes.
- `tol=1e-15` define o valor padrão da tolerância para avaliar pivôs pequenos.
- `n = len(b)` considera que o número de equações é igual ao número de incógnitas.
- `M` é uma cópia de `A`, com todos os elementos convertidos para `float`.
- `v` é uma cópia de `b`, também convertida para `float`.

Usar cópias evita modificar diretamente as listas recebidas como argumentos. O algoritmo pressupõe que `A` seja quadrada e compatível com o tamanho de `b`; a função original não valida explicitamente essas condições.

### 6.2 Laço principal da eliminação

```python
for k in range(n - 1):
```

`k` identifica a etapa atual e a coluna do pivô. O laço vai até a penúltima linha porque, depois de eliminar as entradas abaixo dela, a última equação já não precisa de uma nova etapa de eliminação.

### 6.3 Busca do pivô

```python
p = max(range(k, n), key=lambda i: abs(M[i][k]))
```

- `range(k, n)` considera a linha atual e todas as linhas abaixo.
- `lambda i: abs(M[i][k])` calcula o módulo do elemento na coluna `k` de cada linha candidata.
- `max` devolve o índice da linha que possui o maior módulo nessa coluna.
- `p` é o índice da linha escolhida como pivô.

Essa busca implementa o pivoteamento parcial.

### 6.4 Verificação do pivô

```python
if abs(M[p][k]) < tol:
    raise ZeroDivisionError(
        f"Matriz singular ou quase singular na coluna {k}."
    )
```

Se o maior candidato da coluna tiver módulo inferior à tolerância, o código interrompe a execução com uma exceção. Isso evita prosseguir com um pivô considerado nulo ou numericamente muito pequeno.

A mensagem inclui a coluna em que o problema foi identificado. Esse teste é um critério numérico prático, não uma prova absoluta de singularidade.

### 6.5 Permutação das linhas

```python
if p != k:
    M[k], M[p] = M[p], M[k]
    v[k], v[p] = v[p], v[k]
```

- A troca só é necessária quando `p` e `k` são diferentes.
- `M[k], M[p] = M[p], M[k]` troca as linhas da matriz.
- `v[k], v[p] = v[p], v[k]` troca os termos independentes correspondentes.

As duas operações mantêm a equivalência do sistema linear.

### 6.6 Cálculo do fator de eliminação

```python
for i in range(k + 1, n):
    m = M[i][k] / M[k][k]
```

O laço percorre todas as linhas abaixo do pivô. A variável `m` guarda o fator necessário para eliminar o elemento `M[i][k]`.

### 6.7 Atualização da matriz e do vetor

```python
for j in range(k + 1, n):
    M[i][j] -= m * M[k][j]
v[i] -= m * v[k]
M[i][k] = 0.0
```

- `j` percorre as colunas à direita do pivô.
- `M[i][j] -= m * M[k][j]` equivale a $M_{ij}=M_{ij}-mM_{kj}$.
- `v[i] -= m * v[k]` aplica a mesma operação ao termo independente.
- `M[i][k] = 0.0` registra explicitamente o zero criado pela eliminação.

Ao final, os elementos abaixo da diagonal principal são zerados progressivamente. A matriz `M` passa a ter estrutura triangular superior.

### 6.8 Verificação do último pivô

```python
if abs(M[n - 1][n - 1]) < tol:
    raise ZeroDivisionError(
        "Matriz singular: último pivô nulo ou muito próximo de zero."
    )
```

O laço de eliminação não verifica o último pivô, pois não há linhas abaixo dele. Como a substituição retroativa divide por esse valor, ele precisa ser examinado antes de calcular a solução.

### 6.9 Inicialização do vetor solução

```python
x = [0.0] * n
```

Cria uma lista com `n` posições inicializadas em zero. Cada posição receberá o valor de uma incógnita.

### 6.10 Cálculo da última incógnita

```python
x[n - 1] = v[n - 1] / M[n - 1][n - 1]
```

A última linha da matriz triangular superior contém apenas uma incógnita. Por isso, ela pode ser isolada diretamente dividindo o termo independente pelo último pivô.

### 6.11 Cálculo das demais incógnitas

```python
for i in range(n - 2, -1, -1):
    soma = v[i]
    for j in range(i + 1, n):
        soma -= M[i][j] * x[j]
    x[i] = soma / M[i][i]
```

- `range(n - 2, -1, -1)` percorre as linhas da penúltima até a primeira, em ordem decrescente.
- `soma = v[i]` inicia o cálculo com o termo independente da linha atual.
- O laço `j` percorre as colunas à direita da diagonal, cujas incógnitas já foram calculadas.
- `soma -= M[i][j] * x[j]` subtrai as contribuições dessas incógnitas.
- `x[i] = soma / M[i][i]` isola a incógnita atual.

### 6.12 Retorno

```python
return x
```

Devolve a lista com a solução calculada. A posição `x[0]` corresponde à primeira incógnita, `x[1]` à segunda e assim por diante.

## 7. Como executar

### Requisitos

- Python 3.
- Nenhuma biblioteca externa é necessária.

Salve o código em um arquivo chamado, por exemplo, `eliminacao_gaussiana.py`. O nome pode ser diferente, desde que seja usado corretamente na importação.

### Exemplo de chamada

```python
from eliminacao_gaussiana import eliminacao_gaussiana_pivoteamento_parcial

A = [
    [2, 1],
    [1, -1]
]

b = [5, 1]

solucao = eliminacao_gaussiana_pivoteamento_parcial(A, b)
print(solucao)
```

A saída será aproximadamente:

```text
[2.0, 1.0]
```

O arquivo original contém a definição da função, mas não um bloco de execução principal. Portanto, executá-lo diretamente apenas define a função; para resolver um sistema, é necessário chamá-la em outro script ou adicionar um exemplo ao final do arquivo.

## 8. Exemplo completo

Considere o sistema:

$$
\begin{cases}
2x+y=5,\\
x-y=1.
\end{cases}
$$

A matriz e o vetor são:

$$
A=\begin{bmatrix}2&1\\1&-1\end{bmatrix},
\qquad b=\begin{bmatrix}5\\1\end{bmatrix}.
$$

Código:

```python
from eliminacao_gaussiana import eliminacao_gaussiana_pivoteamento_parcial

A = [
    [2, 1],
    [1, -1]
]

b = [5, 1]

x = eliminacao_gaussiana_pivoteamento_parcial(A, b)

print("Solução:")
for i, valor in enumerate(x, start=1):
    print(f"x{i} = {valor:.6f}")
```

Saída esperada:

```text
Solução:
x1 = 2.000000
x2 = 1.000000
```

Verificação manual:

$$
2(2)+1=5,\qquad 2-1=1.
$$

As duas equações são satisfeitas.

## 9. Verificação da solução

Uma solução aproximada pode ser avaliada substituindo o vetor calculado no sistema original. O resíduo é definido por

$$
r=Ax-b.
$$

Se a solução for exata, $r$ será o vetor nulo. Em cálculos com ponto flutuante, espera-se que as componentes do resíduo sejam pequenas, dentro da precisão e da escala do problema.

Exemplo usando apenas Python:

```python
from eliminacao_gaussiana import eliminacao_gaussiana_pivoteamento_parcial

def calcular_residuo(A, b, x):
    residuo = []
    for i in range(len(A)):
        valor = sum(A[i][j] * x[j] for j in range(len(x)))
        residuo.append(valor - b[i])
    return residuo

A = [
    [2, 1],
    [1, -1]
]
b = [5, 1]

x = eliminacao_gaussiana_pivoteamento_parcial(A, b)
r = calcular_residuo(A, b, x)

print("Solução:", x)
print("Resíduo:", r)
```

Um resíduo próximo de zero indica que a solução calculada satisfaz aproximadamente o sistema original. Entretanto, um resíduo pequeno não garante sozinho que a solução seja pouco sensível a perturbações: essa sensibilidade também depende do condicionamento da matriz.

## 10. Tolerância e limitações

O valor padrão `tol=1e-15` equivale a $10^{-15}$. O programa considera um pivô suspeito quando

```python
abs(pivo) < tol
```

A tolerância absoluta não se ajusta automaticamente à escala dos elementos da matriz. Por isso, seu significado pode variar conforme os dados de entrada.

Limitações da versão atual:

1. Não valida explicitamente se `A` é quadrada.
2. Não verifica se o número de linhas de `A` coincide com o comprimento de `b`.
3. Não verifica se todas as linhas de `A` têm o mesmo tamanho.
4. Não verifica se os dados são finitos, como `NaN` ou infinito.
5. Não lida separadamente com sistemas sem solução ou com infinitas soluções; pode interromper o cálculo ao encontrar um pivô pequeno.
6. Usa números de ponto flutuante, sujeitos a erros de arredondamento.
7. A matriz pode ser mal condicionada mesmo quando todos os pivôs passam pelo teste de tolerância.
8. A função pressupõe `n > 0`; um vetor `b` vazio causa erro de índice na verificação do último pivô.

O pivoteamento parcial melhora o comportamento numérico em muitos casos, mas não garante precisão para qualquer matriz.

## 11. Complexidade computacional

Para uma matriz quadrada de dimensão $n$:

- A eliminação gaussiana realiza, em geral, $O(n^3)$ operações.
- A substituição retroativa exige $O(n^2)$ operações.
- O armazenamento da matriz e dos vetores exige $O(n^2)$ memória.

O custo total é dominado pela eliminação, portanto é da ordem de $O(n^3)$ operações.

Para matrizes grandes ou aplicações de produção, bibliotecas como NumPy e SciPy oferecem solvers numéricos otimizados e amplamente testados. Para resolver um único sistema $Ax=b$, não é necessário calcular a inversa de $A$.

## 12. Possíveis melhorias

Ideias para evoluir o projeto:

- adicionar validação de dimensões e entradas vazias;
- fornecer mensagens de erro específicas para sistemas incompatíveis;
- usar uma tolerância relativa em função da escala da matriz;
- implementar testes automatizados para sistemas com e sem necessidade de pivoteamento;
- testar matrizes singulares e quase singulares;
- calcular e exibir o resíduo da solução;
- separar a eliminação da substituição retroativa em funções distintas;
- adicionar *docstrings*, exemplos e anotações de tipo;
- comparar resultados com `numpy.linalg.solve` em testes;
- estudar técnicas de escalonamento e condicionamento numérico.

## 13. Referências

Tópicos recomendados para estudo:

- Sistemas lineares e operações elementares sobre linhas.
- Eliminação gaussiana.
- Pivoteamento parcial.
- Substituição retroativa.
- Erros de arredondamento, estabilidade e condicionamento numérico.

Documentação oficial:

- Python: https://docs.python.org/3/
- NumPy: https://numpy.org/doc/
- SciPy: https://docs.scipy.org/doc/scipy/

---

## Licença

Este README não define uma licença para o código. Se o repositório for público e você quiser permitir reutilização, adicione um arquivo de licença e escolha os termos adequados ao seu projeto.
