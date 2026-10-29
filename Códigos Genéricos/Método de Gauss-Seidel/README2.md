# Método de Gauss-Seidel para Sistemas Lineares

## 1. Objetivo

Este projeto implementa o **método iterativo de Gauss-Seidel**, utilizado para encontrar soluções aproximadas de sistemas de equações lineares.

Um sistema linear pode ser representado matricialmente por:

$$
A\mathbf{x}=\mathbf{b}
$$

Em que:

- $A$ é a matriz dos coeficientes;
- $\mathbf{x}$ é o vetor das incógnitas;
- $\mathbf{b}$ é o vetor dos termos independentes.

O método de Gauss-Seidel começa com uma estimativa inicial para as incógnitas e a atualiza sucessivamente até que a diferença entre duas aproximações consecutivas seja inferior à tolerância estabelecida, ou até que o número máximo de iterações seja atingido.

A função implementada neste projeto é:

```python
def gauss_seidel(A, b, x0=None, tol=1e-6, max_it=100):
```

Ela recebe a matriz dos coeficientes, o vetor dos termos independentes e parâmetros de controle da iteração. Ao final, retorna o vetor solução aproximado e o número de iterações realizadas.

---

## 2. Fundamentação matemática

### 2.1. Sistema de equações lineares

Considere um sistema linear com $n$ equações e $n$ incógnitas:

$$
\begin{cases}
a_{11}x_1+a_{12}x_2+\cdots+a_{1n}x_n=b_1\\
a_{21}x_1+a_{22}x_2+\cdots+a_{2n}x_n=b_2\\
\vdots\\
a_{n1}x_1+a_{n2}x_2+\cdots+a_{nn}x_n=b_n
\end{cases}
$$

O método de Gauss-Seidel reorganiza cada equação para isolar sua incógnita principal.

Para uma equação de índice $i$:

$$
a_{ii}x_i+\sum_{j\ne i}a_{ij}x_j=b_i
$$

Isolando $x_i$:

$$
x_i=\frac{1}{a_{ii}}
\left(b_i-\sum_{j\ne i}a_{ij}x_j\right)
$$

Essa expressão permite calcular uma nova aproximação para cada incógnita.

### 2.2. Fórmula iterativa de Gauss-Seidel

A fórmula geral do método é:

$$
x_i^{(k+1)}
=
\frac{1}{a_{ii}}
\left[
b_i
-\sum_{j=1}^{i-1}a_{ij}x_j^{(k+1)}
-\sum_{j=i+1}^{n}a_{ij}x_j^{(k)}
\right]
$$

Em que:

- $x_i^{(k+1)}$ é a nova aproximação da incógnita $i$;
- $x_j^{(k+1)}$ representa os valores já atualizados na iteração atual;
- $x_j^{(k)}$ representa os valores ainda não atualizados da iteração anterior;
- $k$ identifica a iteração.

A principal característica do método é que **cada valor novo é utilizado imediatamente nos cálculos seguintes da mesma iteração**.

Isso diferencia Gauss-Seidel de métodos que calculam todas as novas incógnitas exclusivamente a partir do vetor anterior, como o método de Jacobi.

### 2.3. Critério de parada

O código calcula a diferença absoluta entre o valor novo e o valor anterior de cada incógnita:

$$
e_i=\left|x_i^{(k+1)}-x_i^{(k)}\right|
$$

Em seguida, determina o maior erro entre todas as incógnitas:

$$
e_{\max}=\max_i e_i
$$

O algoritmo considera que atingiu o critério de parada quando:

$$
e_{\max}<\text{tol}
$$

A tolerância `tol` define a diferença máxima aceita entre aproximações consecutivas.

É importante observar que esse critério mede a variação das incógnitas entre iterações. Ele não calcula diretamente o resíduo do sistema $A\mathbf{x}-\mathbf{b}$.

---

## 3. Parâmetros da função

A função possui cinco parâmetros:

```python
gauss_seidel(A, b, x0=None, tol=1e-6, max_it=100)
```

| Parâmetro | Descrição |
|---|---|
| `A` | Matriz quadrada dos coeficientes do sistema linear. |
| `b` | Lista com os termos independentes das equações. |
| `x0` | Vetor com a estimativa inicial das incógnitas. Se não for informado, o código utiliza zeros. |
| `tol` | Tolerância utilizada no critério de convergência. O padrão é `1e-6`. |
| `max_it` | Número máximo de iterações. O padrão é `100`. |

### Exemplo dos parâmetros

```python
A = [
    [4, -1, 0],
    [-1, 4, -1],
    [0, -1, 3]
]

b = [15, 10, 10]

x0 = [0, 0, 0]

tol = 1e-6
max_it = 100
```

Nesse exemplo, `A` contém os coeficientes do sistema, `b` contém os termos independentes e `x0` fornece a estimativa inicial.

---

## 4. Explicação detalhada do código

### 4.1. Definição da função e inicialização

```python
def gauss_seidel(A, b, x0=None, tol=1e-6, max_it=100):
    n = len(b)
    x = [0.0] * n if x0 is None else [float(val) for val in x0]
```

A primeira linha define a função e seus parâmetros.

O comando:

```python
n = len(b)
```

Determina o número de equações do sistema a partir do tamanho do vetor `b`.

O número `n` será utilizado para controlar os laços de repetição.

A linha seguinte inicializa o vetor das incógnitas:

```python
x = [0.0] * n if x0 is None else [float(val) for val in x0]
```

Essa expressão condicional possui dois casos:

- Se `x0` for `None`, o vetor `x` será inicializado com zeros.
- Se `x0` for fornecido, seus valores serão convertidos para `float` e utilizados como estimativa inicial.

Por exemplo, se `n = 3` e `x0` não for informado, o vetor inicial será:

$$
\mathbf{x}^{(0)}=
\begin{bmatrix}
0\\
0\\
0
\end{bmatrix}
$$

A escolha do vetor inicial pode influenciar a velocidade de convergência do método.

### 4.2. Verificação dos elementos da diagonal

```python
for i in range(n):
    if abs(A[i][i]) < 1e-15:
        raise ZeroDivisionError(
            f"Elemento diagonal nulo ou muito próximo de zero na linha {i}."
        )
```

Esse bloco verifica os elementos da diagonal principal da matriz.

O elemento `A[i][i]` corresponde ao coeficiente da incógnita que será isolada na equação de índice `i`.

Na fórmula de Gauss-Seidel, esse elemento aparece no denominador:

$$
x_i=\frac{\text{expressão}}{a_{ii}}
$$

Portanto, se $a_{ii}=0$, a divisão não pode ser realizada.

O código também considera problemáticos valores cujo módulo seja inferior a $10^{-15}$:

```python
if abs(A[i][i]) < 1e-15:
```

A função `abs()` calcula o valor absoluto do elemento.

Caso o elemento diagonal seja nulo ou muito próximo de zero, o comando:

```python
raise ZeroDivisionError(...)
```

interrompe a execução e apresenta uma mensagem de erro indicando a linha problemática.

**Importante:** essa verificação evita divisões por valores diagonais muito pequenos, mas não garante que o método convergirá. Além disso, a implementação não realiza uma validação completa das dimensões da matriz.

### 4.3. Laço principal das iterações

```python
for it in range(1, max_it + 1):
    maior_erro = 0.0
```

O laço principal executa as iterações do método.

A expressão:

```python
range(1, max_it + 1)
```

faz com que o contador `it` assuma valores de 1 até `max_it`, inclusive.

Por exemplo, se `max_it = 100`, o algoritmo poderá realizar até 100 iterações.

No início de cada iteração:

```python
maior_erro = 0.0
```

O maior erro é inicializado com zero. Durante a atualização das incógnitas, esse valor será comparado com os erros individuais para determinar a maior variação encontrada na iteração.

### 4.4. Laço de atualização das incógnitas

```python
for i in range(n):
    x_ant = x[i]
    soma = float(b[i])
```

Esse laço percorre todas as equações do sistema.

O comando:

```python
x_ant = x[i]
```

Armazena o valor atual da incógnita antes de sua atualização. Esse valor será utilizado posteriormente para calcular o erro.

Em seguida:

```python
soma = float(b[i])
```

Inicializa a variável `soma` com o termo independente da equação atual.

Matematicamente, o cálculo começa com:

$$
S=b_i
$$

Depois, os produtos dos coeficientes pelas demais incógnitas serão subtraídos desse valor.

### 4.5. Subtração dos termos já atualizados

```python
for j in range(i):
    soma -= A[i][j] * x[j]
```

Esse laço percorre os índices anteriores ao índice `i`.

Para cada índice `j`, o código subtrai:

$$
a_{ij}x_j
$$

da variável `soma`.

Como as incógnitas são atualizadas sequencialmente, os valores `x[j]` com $j<i$ já foram calculados na iteração atual.

Assim, o código utiliza as aproximações mais recentes disponíveis, conforme a característica do método de Gauss-Seidel.

Após esse laço, a expressão representa:

$$
S=b_i-\sum_{j=0}^{i-1}a_{ij}x_j
$$

A implementação utiliza índices iniciados em zero, enquanto a notação matemática costuma identificar as incógnitas a partir de 1.

### 4.6. Subtração dos termos ainda não atualizados

```python
for j in range(i + 1, n):
    soma -= A[i][j] * x[j]
```

Esse laço percorre os índices posteriores a `i`.

Para cada termo, subtrai:

$$
a_{ij}x_j
$$

Como as incógnitas de índices maiores ainda não foram atualizadas na iteração atual, seus valores correspondem às aproximações anteriores.

Depois dos dois laços, a variável `soma` contém:

$$
S=b_i
-\sum_{j=0}^{i-1}a_{ij}x_j^{(k+1)}
-\sum_{j=i+1}^{n-1}a_{ij}x_j^{(k)}
$$

Essa expressão corresponde ao numerador da fórmula de Gauss-Seidel.

### 4.7. Atualização da incógnita

```python
x[i] = soma / A[i][i]
```

O código divide a soma pelo elemento diagonal da linha atual.

Matematicamente:

$$
x_i^{(k+1)}=\frac{S}{a_{ii}}
$$

Esse resultado substitui o valor anterior de `x[i]`.

Como o vetor `x` é atualizado diretamente, os valores recém-calculados ficam disponíveis para os cálculos seguintes da mesma iteração.

Essa atualização imediata é a característica central do método de Gauss-Seidel.

### 4.8. Cálculo do erro individual

```python
erro = abs(x[i] - x_ant)
```

Depois de atualizar a incógnita, o código calcula a diferença absoluta entre o valor novo e o valor anterior:

$$
e_i=\left|x_i^{(k+1)}-x_i^{(k)}\right|
$$

A função `abs()` garante que o erro seja não negativo.

Por exemplo, se o valor anterior era 2,5 e o novo valor é 2,8:

$$
e_i=|2,8-2,5|=0,3
$$

Esse erro é utilizado para verificar se as aproximações estão se estabilizando.

### 4.9. Determinação do maior erro

```python
if erro > maior_erro:
    maior_erro = erro
```

O código compara o erro individual com o maior erro registrado até aquele momento.

Se o erro atual for maior, o valor de `maior_erro` é atualizado.

Ao final do laço que percorre todas as incógnitas, essa variável contém:

$$
e_{\max}=\max_i e_i
$$

Ou seja, a maior variação absoluta entre as aproximações consecutivas de todas as incógnitas naquela iteração.

### 4.10. Verificação da convergência

```python
if maior_erro < tol:
    return x, it
```

Depois de atualizar todas as incógnitas, o código compara o maior erro com a tolerância.

O critério é:

$$
e_{\max}<\text{tol}
$$

Se a condição for satisfeita, a função retorna imediatamente:

```python
return x, it
```

O retorno contém dois elementos:

1. `x`: vetor solução aproximado.
2. `it`: número de iterações realizadas.

Por exemplo, o retorno pode ser interpretado como:

```python
solucao, iteracoes = gauss_seidel(A, b)
```

Nesse caso, `solucao` recebe o vetor das incógnitas e `iteracoes` recebe a quantidade de iterações executadas até a parada.

É importante destacar que o critério verifica a diferença entre aproximações consecutivas, e não diretamente a precisão da solução em relação ao sistema original.

### 4.11. Retorno ao atingir o número máximo de iterações

```python
return x, max_it
```

Se o laço terminar sem satisfazer o critério de convergência, a função retorna o vetor calculado após `max_it` iterações.

O retorno possui a mesma estrutura do caso de convergência:

- Primeiro elemento: vetor aproximado das incógnitas.
- Segundo elemento: número de iterações realizadas.

**Atenção:** nesse caso, o código não informa explicitamente que a convergência não foi atingida. O usuário deve verificar a qualidade do resultado antes de considerá-lo uma solução suficientemente precisa.

---

## 5. Exemplo completo de utilização

Considere o seguinte sistema linear:

$$
\begin{cases}
4x_1-x_2=15\\
-x_1+4x_2-x_3=10\\
-x_2+3x_3=10
\end{cases}
$$

A matriz dos coeficientes é:

$$
A=
\begin{bmatrix}
4&-1&0\\
-1&4&-1\\
0&-1&3
\end{bmatrix}
$$

O vetor dos termos independentes é:

$$
\mathbf{b}=
\begin{bmatrix}
15\\
10\\
10
\end{bmatrix}
$$

### 5.1. Código Python

```python
def gauss_seidel(A, b, x0=None, tol=1e-6, max_it=100):
    n = len(b)
    x = [0.0] * n if x0 is None else [float(val) for val in x0]

    for i in range(n):
        if abs(A[i][i]) < 1e-15:
            raise ZeroDivisionError(
                f"Elemento diagonal nulo ou muito próximo de zero na linha {i}."
            )

    for it in range(1, max_it + 1):
        maior_erro = 0.0

        for i in range(n):
            x_ant = x[i]
            soma = float(b[i])

            for j in range(i):
                soma -= A[i][j] * x[j]

            for j in range(i + 1, n):
                soma -= A[i][j] * x[j]

            x[i] = soma / A[i][i]
            erro = abs(x[i] - x_ant)

            if erro > maior_erro:
                maior_erro = erro

        if maior_erro < tol:
            return x, it

    return x, max_it


A = [
    [4, -1, 0],
    [-1, 4, -1],
    [0, -1, 3]
]

b = [15, 10, 10]

x0 = [0, 0, 0]

solucao, iteracoes = gauss_seidel(
    A,
    b,
    x0=x0,
    tol=1e-6,
    max_it=100
)

print("Solução aproximada:", solucao)
print("Número de iterações:", iteracoes)
```

### 5.2. Interpretação

A matriz e o vetor são fornecidos à função, juntamente com a estimativa inicial:

$$
\mathbf{x}^{(0)}=
\begin{bmatrix}
0\\
0\\
0
\end{bmatrix}
$$

O algoritmo calcula sucessivamente os valores de $x_1$, $x_2$ e $x_3$, utilizando os valores atualizados assim que ficam disponíveis.

A solução exata desse sistema é:

$$
x_1=5,\qquad x_2=5,\qquad x_3=5
$$

O resultado numérico deve se aproximar desses valores, respeitando a tolerância e as condições de convergência.

---

## 6. Verificação da solução

Uma forma de verificar a qualidade da solução encontrada é calcular o resíduo do sistema.

O resíduo é definido por:

$$
\mathbf{r}=A\mathbf{x}-\mathbf{b}
$$

Se a solução for exata, espera-se:

$$
\mathbf{r}=\mathbf{0}
$$

Em uma solução numérica aproximada, os componentes do resíduo devem ser suficientemente pequenos para a aplicação considerada.

O código abaixo pode ser utilizado para verificar o resultado:

```python
def calcular_residuo(A, b, x):
    residuo = []

    for i in range(len(b)):
        soma = 0.0

        for j in range(len(x)):
            soma += A[i][j] * x[j]

        residuo.append(soma - b[i])

    return residuo


residuo = calcular_residuo(A, b, solucao)

print("Resíduo:", residuo)
print(
    "Maior resíduo absoluto:",
    max(abs(valor) for valor in residuo)
)
```

O resíduo permite avaliar se o vetor calculado satisfaz aproximadamente as equações originais.

É possível que o critério baseado na diferença entre iterações seja satisfeito sem que o resíduo seja suficientemente pequeno. Por isso, verificar ambos os indicadores torna a análise mais confiável.

---

## 7. Condições de convergência

O método de Gauss-Seidel não converge necessariamente para qualquer matriz.

Uma condição suficiente comum para a convergência é que a matriz seja **estritamente diagonal dominante por linhas**.

Isso significa que, em cada linha $i$:

$$
|a_{ii}|>\sum_{j\ne i}|a_{ij}|
$$

Em palavras, o módulo do elemento da diagonal principal deve ser maior que a soma dos módulos dos demais elementos da mesma linha.

Para a matriz do exemplo:

$$
A=
\begin{bmatrix}
4&-1&0\\
-1&4&-1\\
0&-1&3
\end{bmatrix}
$$

A primeira linha satisfaz:

$$
|4|>|-1|+|0|
$$

$$
4>1
$$

A segunda linha satisfaz:

$$
|4|>|-1|+|-1|
$$

$$
4>2
$$

A terceira linha satisfaz:

$$
|3|>|0|+|-1|
$$

$$
3>1
$$

Portanto, essa matriz é estritamente diagonal dominante por linhas, o que garante a convergência do método para qualquer vetor inicial, sob as condições usuais do problema linear.

A diagonal dominante é uma condição suficiente, mas não necessária. Existem matrizes que não são diagonalmente dominantes e ainda assim permitem a convergência.

O código original não verifica a dominância diagonal. Ele apenas verifica se os elementos da diagonal principal são suficientemente diferentes de zero.

---

## 8. Limitações e cuidados

### 8.1. Elementos diagonais nulos ou muito pequenos

O código interrompe a execução se encontrar um elemento diagonal cujo módulo seja inferior a $10^{-15}$.

Entretanto, elementos diagonais diferentes de zero não garantem a convergência.

### 8.2. Dimensões dos dados

A implementação pressupõe que `A` seja uma matriz quadrada compatível com `b` e que `x0`, quando informado, tenha o tamanho adequado.

O código não verifica explicitamente todas essas condições. Entradas incompatíveis podem provocar erros de índice ou resultados incorretos.

### 8.3. Dependência da matriz

O método pode convergir rapidamente, lentamente ou não convergir, dependendo das propriedades da matriz dos coeficientes.

A escolha da estimativa inicial também pode influenciar o comportamento numérico.

### 8.4. Critério de parada

O critério implementado utiliza a maior diferença absoluta entre duas aproximações consecutivas:

$$
e_{\max}<\text{tol}
$$

Ele não utiliza diretamente o resíduo $A\mathbf{x}-\mathbf{b}$, nem considera uma tolerância relativa.

Para aplicações que exigem maior confiabilidade, é recomendável verificar também o resíduo e, quando apropriado, utilizar critérios de parada relativos.

### 8.5. Ausência de indicador explícito de convergência

Quando a função atinge `max_it`, ela retorna a última aproximação e o número máximo de iterações.

Não há um indicador separado para informar se o critério de convergência foi atingido. Por isso, é importante verificar o resíduo ou adaptar a função para informar explicitamente o estado da execução.

### 8.6. Uso de ponto flutuante

Os cálculos são realizados com números de ponto flutuante. Assim, podem ocorrer pequenas diferenças em relação à solução exata devido ao arredondamento numérico.

---

## 9. Complexidade computacional

Considere um sistema com $n$ equações e $n$ incógnitas.

Em cada iteração, o código percorre as linhas e os coeficientes da matriz. O número de operações é proporcional a $n^2$.

Assim, o custo computacional aproximado por iteração é:

$$
O(n^2)
$$

Se o método realizar $k$ iterações, o custo total é aproximadamente:

$$
O(kn^2)
$$

Essa estimativa considera a implementação apresentada, com a matriz armazenada como uma lista de listas e sem incluir operações adicionais de validação.

O método pode ser interessante para sistemas em que a solução iterativa seja adequada, especialmente quando a matriz apresenta propriedades favoráveis à convergência.

---

## 10. Possíveis melhorias

A implementação pode ser aprimorada com as seguintes funcionalidades:

1. **Validação das dimensões:** verificar se a matriz é quadrada e compatível com os vetores `b` e `x0`.
2. **Verificação da convergência:** informar explicitamente se o critério de parada foi satisfeito ou se o número máximo de iterações foi atingido.
3. **Cálculo do resíduo:** avaliar $A\mathbf{x}-\mathbf{b}$ ao final da execução.
4. **Tolerância relativa:** considerar a escala das incógnitas e dos dados do sistema.
5. **Aviso sobre convergência:** verificar condições suficientes, como a dominância diagonal estrita, e informar ao usuário quando elas não forem satisfeitas.
6. **Registro das iterações:** armazenar os valores aproximados em cada iteração para estudar a evolução do método.
7. **Tratamento de entradas inválidas:** verificar valores não finitos, coeficientes inadequados e parâmetros de entrada.

Essas melhorias são sugestões para versões futuras e não estão implementadas no código original.

---

## 11. Resumo do algoritmo

O funcionamento do método pode ser resumido nas seguintes etapas:

1. Receber a matriz $A$ e o vetor $\mathbf{b}$.
2. Definir a estimativa inicial das incógnitas.
3. Verificar os elementos da diagonal principal.
4. Iniciar uma nova iteração.
5. Atualizar cada incógnita usando os valores mais recentes disponíveis.
6. Calcular a diferença absoluta entre o valor novo e o anterior de cada incógnita.
7. Determinar a maior diferença da iteração.
8. Comparar essa diferença com a tolerância.
9. Se a tolerância for satisfeita, retornar a solução aproximada e o número de iterações.
10. Caso contrário, continuar até atingir o limite máximo de iterações.
11. Se o limite for atingido, retornar a última aproximação disponível.

O método de Gauss-Seidel é uma ferramenta importante para a resolução iterativa de sistemas lineares. Sua principal característica é utilizar imediatamente os valores recém-calculados, o que pode favorecer a convergência em comparação com outros métodos iterativos.

Para obter resultados confiáveis, é importante analisar as propriedades da matriz, escolher parâmetros adequados e verificar a solução por meio do resíduo do sistema.
