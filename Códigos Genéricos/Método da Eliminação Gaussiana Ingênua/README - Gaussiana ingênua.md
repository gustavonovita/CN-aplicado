# Método da Eliminação Gaussiana Ingênua

## Descrição

Este projeto implementa o **método da eliminação gaussiana ingênua** para resolver sistemas de equações lineares da forma:

\[
A\mathbf{x}=\mathbf{b}
\]

Em que:

- \(A\) é a matriz dos coeficientes;
- \(\mathbf{x}\) é o vetor das incógnitas;
- \(\mathbf{b}\) é o vetor dos termos independentes.

O método transforma o sistema original em um sistema triangular superior por meio de operações elementares nas linhas da matriz. Em seguida, calcula as incógnitas utilizando a substituição retroativa.

A implementação não realiza pivoteamento: em cada etapa, utiliza diretamente o elemento da diagonal principal como pivô. Por isso, é importante que os pivôs não sejam nulos nem muito próximos de zero.

## Fundamentação matemática

O método é dividido em duas etapas principais:

### 1. Eliminação progressiva

O objetivo é zerar os elementos abaixo da diagonal principal da matriz \(A\), transformando-a em uma matriz triangular superior.

Em cada etapa \(k\), calcula-se o multiplicador:

\[
m_{ik}=\frac{a_{ik}}{a_{kk}}
\]

Em seguida, atualizam-se os elementos da linha:

\[
a_{ij}\leftarrow a_{ij}-m_{ik}a_{kj}
\]

e o termo independente correspondente:

\[
b_i\leftarrow b_i-m_{ik}b_k
\]

O elemento abaixo do pivô é definido como zero.

### 2. Substituição retroativa

Após a eliminação progressiva, o sistema apresenta uma matriz triangular superior. A solução é obtida de baixo para cima, começando pela última equação.

A expressão geral para calcular a incógnita \(x_i\) é:

\[
x_i=
\frac{
b_i-\displaystyle\sum_{j=i+1}^{n-1}a_{ij}x_j
}{
a_{ii}
}
\]

Dessa forma, cada incógnita é calculada utilizando os valores que já foram determinados nas equações anteriores da substituição retroativa.

## Código implementado

```python
def eliminacao_gaussiana_ingenua(A, b, tol=1e-15):
    n = len(b)
    M = [[float(val) for val in linha] for linha in A]
    v = [float(val) for val in b]

    # Eliminação progressiva
    for k in range(n - 1):
        if abs(M[k][k]) < tol:
            raise ZeroDivisionError(f"Pivô nulo ou muito próximo de zero na linha {k}.")
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k + 1, n):
                M[i][j] -= m * M[k][j]
            v[i] -= m * v[k]
            M[i][k] = 0.0

    if abs(M[n - 1][n - 1]) < tol:
        raise ZeroDivisionError("Último pivô nulo ou muito próximo de zero.")

    # Substituição retroativa
    x = [0.0] * n
    x[n - 1] = v[n - 1] / M[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = v[i]
        for j in range(i + 1, n):
            soma -= M[i][j] * x[j]
        x[i] = soma / M[i][i]

    return x
```

## Explicação do código passo a passo

### 1. Definição da função

```python
def eliminacao_gaussiana_ingenua(A, b, tol=1e-15):
```

Define a função responsável por resolver o sistema linear.

Parâmetros:

- `A`: matriz quadrada dos coeficientes;
- `b`: vetor dos termos independentes;
- `tol`: tolerância utilizada para detectar pivôs nulos ou muito próximos de zero. Seu valor padrão é \(10^{-15}\).

A função retorna um vetor contendo os valores calculados para as incógnitas.

### 2. Determinação da dimensão do sistema

```python
n = len(b)
```

Obtém o número de equações a partir do comprimento do vetor `b`.

A implementação pressupõe que a matriz `A` seja quadrada, de dimensão \(n\times n\), e compatível com o vetor de termos independentes.

### 3. Conversão dos dados para números reais

```python
M = [[float(val) for val in linha] for linha in A]
v = [float(val) for val in b]
```

Cria cópias da matriz e do vetor, convertendo seus elementos para o tipo `float`.

As operações são realizadas em `M` e `v`, preservando as estruturas originais recebidas pela função.

### 4. Início da eliminação progressiva

```python
for k in range(n - 1):
```

Percorre as colunas de pivoteamento, da primeira até a penúltima.

Em cada etapa, o elemento diagonal `M[k][k]` é utilizado como pivô para eliminar os elementos abaixo dele.

### 5. Verificação do pivô

```python
if abs(M[k][k]) < tol:
    raise ZeroDivisionError(
        f"Pivô nulo ou muito próximo de zero na linha {k}."
    )
```

Verifica se o módulo do pivô é inferior à tolerância definida.

Se isso acontecer, a função interrompe a execução e lança uma exceção `ZeroDivisionError`.

Essa verificação evita a divisão por um valor considerado numericamente inadequado. Como não há pivoteamento, a função não tenta trocar linhas para encontrar outro pivô.

### 6. Percorrer as linhas abaixo do pivô

```python
for i in range(k + 1, n):
```

Percorre as linhas localizadas abaixo da linha do pivô.

O objetivo é eliminar o elemento da coluna `k` de cada uma dessas linhas.

### 7. Cálculo do multiplicador

```python
m = M[i][k] / M[k][k]
```

Calcula o multiplicador utilizado para zerar o elemento `M[i][k]`.

Matematicamente:

\[
m_{ik}=\frac{a_{ik}}{a_{kk}}
\]

Esse valor determina quanto da linha do pivô deve ser subtraído da linha atual.

### 8. Atualização dos coeficientes

```python
for j in range(k + 1, n):
    M[i][j] -= m * M[k][j]
```

Atualiza os coeficientes da linha atual utilizando a operação elementar:

\[
L_i\leftarrow L_i-mL_k
\]

O laço começa em `k + 1`, pois o elemento da coluna `k` será definido como zero explicitamente.

### 9. Atualização do termo independente

```python
v[i] -= m * v[k]
```

Aplica a mesma operação ao termo independente da equação.

Essa atualização é necessária para manter a equivalência entre o sistema original e o sistema transformado.

### 10. Registro do elemento eliminado

```python
M[i][k] = 0.0
```

Define explicitamente como zero o elemento que foi eliminado.

Isso torna clara a estrutura triangular superior que está sendo construída.

### 11. Verificação do último pivô

```python
if abs(M[n - 1][n - 1]) < tol:
    raise ZeroDivisionError(
        "Último pivô nulo ou muito próximo de zero."
    )
```

Após a eliminação progressiva, verifica se o último elemento da diagonal principal é suficientemente diferente de zero.

Esse elemento será utilizado para calcular a última incógnita durante a substituição retroativa.

### 12. Inicialização do vetor solução

```python
x = [0.0] * n
```

Cria um vetor de tamanho `n` para armazenar as incógnitas calculadas.

Inicialmente, todos os elementos são definidos como zero.

### 13. Cálculo da última incógnita

```python
x[n - 1] = v[n - 1] / M[n - 1][n - 1]
```

Calcula a última incógnita dividindo o termo independente da última equação pelo coeficiente diagonal correspondente.

Como a matriz já é triangular superior, a última equação contém apenas uma incógnita.

### 14. Percorrer as equações em ordem inversa

```python
for i in range(n - 2, -1, -1):
```

Percorre as linhas da penúltima até a primeira.

A ordem inversa permite utilizar as incógnitas já calculadas nas linhas inferiores.

### 15. Inicialização da soma

```python
soma = v[i]
```

Inicializa a variável `soma` com o termo independente da equação atual.

A partir desse valor, serão subtraídas as contribuições das incógnitas já conhecidas.

### 16. Subtração das contribuições conhecidas

```python
for j in range(i + 1, n):
    soma -= M[i][j] * x[j]
```

Subtrai os produtos entre os coeficientes da linha atual e as incógnitas já calculadas.

Matematicamente, essa etapa corresponde a:

\[
b_i-\sum_{j=i+1}^{n-1}a_{ij}x_j
\]

### 17. Cálculo da incógnita atual

```python
x[i] = soma / M[i][i]
```

Divide a soma calculada pelo elemento diagonal da linha atual, determinando o valor da incógnita `x[i]`.

Esse processo é repetido até que todas as incógnitas tenham sido calculadas.

### 18. Retorno da solução

```python
return x
```

Retorna o vetor com os valores encontrados para as incógnitas do sistema linear.

## Exemplo de utilização

Considere o sistema de equações:

\[
\begin{cases}
2x+y-z=8\\
-3x-y+2z=-11\\
-2x+y+2z=-3
\end{cases}
\]

A representação matricial é:

\[
A=
\begin{bmatrix}
2&1&-1\\
-3&-1&2\\
-2&1&2
\end{bmatrix}
\]

\[
b=
\begin{bmatrix}
8\\
-11\\
-3
\end{bmatrix}
\]

O código para resolver esse sistema é:

```python
A = [
    [2, 1, -1],
    [-3, -1, 2],
    [-2, 1, 2]
]

b = [8, -11, -3]

solucao = eliminacao_gaussiana_ingenua(A, b)

print("Solução:", solucao)
```

Resultado esperado:

```text
Solução: [2.0, 3.0, -1.0]
```

Portanto, a solução do sistema é:

\[
\boxed{x=2,\qquad y=3,\qquad z=-1}
\]

## Verificação da solução

Uma forma de verificar a solução é substituir os valores encontrados nas equações originais.

Primeira equação:

\[
2(2)+3-(-1)=8
\]

\[
4+3+1=8
\]

Segunda equação:

\[
-3(2)-3+2(-1)=-11
\]

\[
-6-3-2=-11
\]

Terceira equação:

\[
-2(2)+3+2(-1)=-3
\]

\[
-4+3-2=-3
\]

As três equações são satisfeitas, confirmando a solução do exemplo.

### Verificação numérica com NumPy

Também é possível calcular o resíduo do sistema:

\[
\mathbf{r}=A\mathbf{x}-\mathbf{b}
\]

Um resíduo próximo de zero indica que a solução calculada satisfaz aproximadamente o sistema original.

```python
import numpy as np

A_np = np.array(A, dtype=float)
b_np = np.array(b, dtype=float)
x_np = np.array(solucao, dtype=float)

residuo = A_np @ x_np - b_np

print("Resíduo:", residuo)
print("Norma do resíduo:", np.linalg.norm(residuo))
```

O resíduo é uma ferramenta de verificação numérica. Um resíduo pequeno, isoladamente, não garante que o sistema esteja bem condicionado.

## Cuidados e limitações

- **Ausência de pivoteamento:** o método utiliza diretamente os elementos da diagonal principal como pivôs, sem realizar trocas de linhas. Um sistema que poderia ser resolvido após uma troca de linhas pode gerar uma exceção nesta implementação.
- **Tolerância absoluta:** a comparação com `tol` não considera diretamente a escala dos coeficientes. Em sistemas com valores muito grandes ou muito pequenos, a tolerância pode exigir ajustes.
- **Dimensões compatíveis:** a função pressupõe uma matriz quadrada e um vetor de termos independentes de dimensão correspondente. O código não realiza uma validação completa dessas dimensões.
- **Sistemas singulares ou mal condicionados:** erros de arredondamento podem afetar a solução, especialmente quando os coeficientes tornam o sistema numericamente sensível.
- **Precisão de ponto flutuante:** a utilização de `float` pode introduzir pequenas diferenças numéricas.
- **Validação da solução:** em aplicações importantes, recomenda-se verificar o resíduo utilizando a matriz e o vetor originais.

## Complexidade computacional

Para um sistema linear de dimensão \(n\times n\), os custos aproximados são:

| Etapa | Complexidade |
|---|---|
| Eliminação progressiva | \(O(n^3)\) |
| Substituição retroativa | \(O(n^2)\) |
| Armazenamento da matriz copiada | \(O(n^2)\) |

O custo total é dominado pela eliminação progressiva, resultando em complexidade computacional da ordem de:

\[
\boxed{O(n^3)}
\]

## Resumo

A função `eliminacao_gaussiana_ingenua(A, b, tol=1e-15)` resolve sistemas lineares por meio de duas etapas: eliminação progressiva e substituição retroativa.

A implementação inclui verificações para pivôs nulos ou muito pequenos, mas não realiza pivoteamento. Por isso, seu uso exige atenção às características numéricas do sistema e à validação da solução encontrada.
