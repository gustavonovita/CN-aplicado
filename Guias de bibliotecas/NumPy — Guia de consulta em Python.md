# NumPy — Guia de consulta em Python

Este arquivo reúne as principais funções do **NumPy**, com explicações e exemplos práticos para consulta no dia a dia.

O NumPy é uma biblioteca para computação numérica em Python. Seu principal objeto é o `ndarray`, um array multidimensional eficiente, usado em cálculos matemáticos, estatística, álgebra linear, simulações e análise de dados.

## Instalação

```bash
python -m pip install numpy matplotlib
```

Se estiver usando um ambiente virtual:

```bash
python -m venv .venv
```

No Windows:

```bash
.venv\Scripts\activate
```

No Linux ou macOS:

```bash
source .venv/bin/activate
```


## Importação

```python
import numpy as np
```

O apelido `np` é a convenção mais utilizada para o NumPy.

# NumPy: arrays

O objeto principal do NumPy é o `ndarray`, um array capaz de armazenar dados numéricos em uma ou mais dimensões.

## Lista Python versus array NumPy

```python
import numpy as np

lista = [1, 2, 3, 4]
array = np.array([1, 2, 3, 4])

print(lista)
print(array)
```

A principal vantagem do array NumPy é permitir operações vetorizadas, sem a necessidade de escrever um `for` para cada elemento:

```python
array = np.array([1, 2, 3, 4])

print(array * 2)
# [2 4 6 8]

print(array + 10)
# [11 12 13 14]
```

Com uma lista Python, `lista * 2` repete a lista duas vezes, em vez de multiplicar os valores.

---

## Criação de arrays

### `np.array()`

Converte uma lista, tupla ou sequência em um array NumPy.

```python
import numpy as np

vetor = np.array([10, 20, 30])
matriz = np.array([[1, 2, 3], [4, 5, 6]])

print(vetor)
print(matriz)
```

### `np.zeros()`

Cria um array preenchido com zeros.

```python
zeros = np.zeros(5)
print(zeros)
# [0. 0. 0. 0. 0.]

matriz_zeros = np.zeros((2, 3))
print(matriz_zeros)
# [[0. 0. 0.]
#  [0. 0. 0.]]
```

### `np.ones()`

Cria um array preenchido com um.

```python
uns = np.ones(4)
matriz_uns = np.ones((2, 2))
```

### `np.full()`

Cria um array preenchido com um valor específico.

```python
matriz = np.full((2, 3), 7)
print(matriz)
# [[7 7 7]
#  [7 7 7]]
```

### `np.empty()`

Cria um array sem inicializar os valores. O conteúdo inicial pode conter valores que já estavam na memória.

```python
array = np.empty((2, 2))
```

Use essa função somente quando os valores forem preenchidos posteriormente.

### `np.arange()`

Cria valores igualmente espaçados em um intervalo. O limite final não é incluído.

```python
print(np.arange(5))
# [0 1 2 3 4]

print(np.arange(2, 10, 2))
# [2 4 6 8]
```

Sintaxe:

```python
np.arange(início, fim, passo)
```

### `np.linspace()`

Cria uma quantidade definida de valores igualmente espaçados entre dois limites, incluindo o limite final por padrão.

```python
valores = np.linspace(0, 10, 5)
print(valores)
# [ 0.   2.5  5.   7.5 10. ]
```

Sintaxe:

```python
np.linspace(início, fim, quantidade)
```

É muito útil para criar o eixo `x` de gráficos.

### `np.eye()` e `np.identity()`

Criam uma matriz identidade.

```python
identidade_1 = np.eye(3)
identidade_2 = np.identity(3)

print(identidade_1)
# [[1. 0. 0.]
#  [0. 1. 0.]
#  [0. 0. 1.]]
```

### Criar arrays a partir de valores aleatórios

```python
rng = np.random.default_rng(seed=42)
array = rng.random(5)
print(array)
```

A seção [Números aleatórios](#números-aleatórios) explica esse assunto com mais detalhes.

---

## Atributos importantes de um array

```python
import numpy as np

matriz = np.array([[1, 2, 3], [4, 5, 6]])

print(matriz.ndim)   # Número de dimensões: 2
print(matriz.shape)  # Formato: (2, 3)
print(matriz.size)   # Quantidade total de elementos: 6
print(matriz.dtype)  # Tipo dos dados, por exemplo: int64
print(matriz.itemsize)  # Quantidade de bytes por elemento
```

### Atributos principais

| Atributo | Descrição |
|---|---|
| `ndim` | Número de dimensões do array |
| `shape` | Tamanho de cada dimensão |
| `size` | Número total de elementos |
| `dtype` | Tipo dos elementos |
| `itemsize` | Quantidade de bytes ocupados por cada elemento |

### Definir o tipo de dado

```python
inteiros = np.array([1, 2, 3], dtype=np.int32)
decimais = np.array([1, 2, 3], dtype=float)
booleanos = np.array([True, False, True], dtype=bool)

print(inteiros.dtype)
print(decimais.dtype)
print(booleanos.dtype)
```

Tipos comuns:

- `np.int32`, `np.int64`: inteiros.
- `np.float32`, `np.float64`: números decimais.
- `np.bool_`: valores booleanos.
- `np.str_`: textos.

### Converter o tipo de um array com `astype()`

```python
valores = np.array([1.2, 2.7, 3.9])
inteiros = valores.astype(int)

print(inteiros)
# [1 2 3]
```

A conversão para inteiro elimina a parte decimal; ela não arredonda automaticamente.

---

## Indexação e fatiamento

### Arrays de uma dimensão

```python
valores = np.array([10, 20, 30, 40, 50])

print(valores[0])    # 10
print(valores[-1])   # 50
print(valores[1:4])  # [20 30 40]
print(valores[:3])   # [10 20 30]
print(valores[::2])  # [10 30 50]
```

A indexação começa em zero, assim como em listas Python.

### Arrays de duas dimensões

```python
matriz = np.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
])

print(matriz[0, 1])  # Linha 0, coluna 1: 2
print(matriz[2, 0])  # Linha 2, coluna 0: 7
print(matriz[1, :])  # Segunda linha: [4 5 6]
print(matriz[:, 1])  # Segunda coluna: [2 5 8]
print(matriz[:2, :2])
# [[1 2]
#  [4 5]]
```

### Indexação booleana

```python
valores = np.array([5, 12, 3, 20, 8])

maiores_que_10 = valores > 10
print(maiores_que_10)
# [False  True False  True False]

print(valores[maiores_que_10])
# [12 20]
```

Forma mais curta:

```python
print(valores[valores > 10])
```

### `np.where()`

Retorna índices ou escolhe valores de acordo com uma condição.

```python
valores = np.array([1, 2, 3, 4, 5])

indices = np.where(valores % 2 == 0)
print(indices)

resultado = np.where(valores >= 3, "aprovado", "reprovado")
print(resultado)
# ['reprovado' 'reprovado' 'aprovado' 'aprovado' 'aprovado']
```

---

## Alteração de arrays

### Alterar elementos

```python
valores = np.array([10, 20, 30, 40])
valores[1] = 99
print(valores)
# [10 99 30 40]
```

### Alterar vários elementos com uma condição

```python
valores = np.array([5, 15, 8, 22, 3])
valores[valores < 10] = 0
print(valores)
# [ 0 15  0 22  0]
```

### `reshape()`

Altera o formato do array sem mudar a quantidade total de elementos.

```python
valores = np.arange(12)
matriz = valores.reshape(3, 4)

print(matriz)
# [[ 0  1  2  3]
#  [ 4  5  6  7]
#  [ 8  9 10 11]]
```

A quantidade de elementos deve ser compatível:

```python
# 12 elementos podem virar 3x4, 2x6, 4x3 etc.
matriz = valores.reshape(2, 6)
```

É possível deixar uma dimensão para o NumPy calcular:

```python
matriz = valores.reshape(3, -1)
```

### `flatten()` e `ravel()`

Transformam um array multidimensional em uma dimensão.

```python
matriz = np.array([[1, 2], [3, 4]])

copia_1d = matriz.flatten()
visao_1d = matriz.ravel()

print(copia_1d)
print(visao_1d)
```

- `flatten()` retorna uma cópia.
- `ravel()` normalmente retorna uma visão quando possível, sendo mais econômico em memória.

### `transpose()` e `.T`

Trocam linhas por colunas.

```python
matriz = np.array([[1, 2, 3], [4, 5, 6]])

print(matriz.T)
# [[1 4]
#  [2 5]
#  [3 6]]
```

### `resize()`

Altera o tamanho do array original.

```python
valores = np.array([1, 2, 3, 4])
valores.resize((2, 2))
print(valores)
```

Diferentemente de `reshape()`, `resize()` pode modificar o próprio array e pode repetir ou remover dados dependendo do novo tamanho.

---

## Combinação e divisão de arrays

### `concatenate()`

Une arrays ao longo de um eixo.

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

resultado = np.concatenate((a, b))
print(resultado)
# [1 2 3 4 5 6]
```

Com matrizes:

```python
m1 = np.array([[1, 2], [3, 4]])
m2 = np.array([[5, 6], [7, 8]])

por_linhas = np.concatenate((m1, m2), axis=0)
por_colunas = np.concatenate((m1, m2), axis=1)
```

### `vstack()` e `hstack()`

```python
a = np.array([1, 2])
b = np.array([3, 4])

print(np.vstack((a, b)))
# [[1 2]
#  [3 4]]

print(np.hstack((a, b)))
# [1 2 3 4]
```

### `split()`

Divide um array em partes.

```python
valores = np.array([1, 2, 3, 4, 5, 6])
partes = np.split(valores, 3)

print(partes)
# [array([1, 2]), array([3, 4]), array([5, 6])]
```

---

## Operações matemáticas

As operações são feitas elemento por elemento quando aplicadas entre arrays do mesmo formato ou entre um array e um escalar.

```python
import numpy as np

a = np.array([10, 20, 30])
b = np.array([2, 4, 5])

print(a + b)  # [12 24 35]
print(a - b)  # [ 8 16 25]
print(a * b)  # [ 20  80 150]
print(a / b)  # [5. 5. 6.]
print(a ** 2) # [100 400 900]
print(a % b)  # [0 0 0]
```

### Funções matemáticas comuns

```python
valores = np.array([1, 4, 9, 16])

print(np.sqrt(valores))       # Raiz quadrada
print(np.square(valores))     # Quadrado
print(np.abs([-3, -1, 2]))    # Valor absoluto
print(np.exp([0, 1, 2]))      # Exponencial de base e
print(np.log([1, np.e]))      # Logaritmo natural
print(np.log10([1, 10, 100])) # Logaritmo na base 10
```

### Arredondamento

```python
valores = np.array([1.234, 2.567, 3.891])

print(np.round(valores, 2)) # [1.23 2.57 3.89]
print(np.floor(valores))    # Arredonda para baixo
print(np.ceil(valores))     # Arredonda para cima
print(np.trunc(valores))    # Remove a parte decimal
```

### Trigonometria

As funções trigonométricas recebem ângulos em radianos.

```python
angulo = np.pi / 2

print(np.sin(angulo))
print(np.cos(angulo))
print(np.tan(np.pi / 4))
```

Conversões:

```python
graus = np.array([0, 90, 180])
radianos = np.deg2rad(graus)

print(np.rad2deg(radianos))
```

### Broadcasting

Broadcasting é o mecanismo que permite operar arrays de formatos compatíveis sem repetir manualmente os valores.

```python
notas = np.array([
    [8.0, 7.0, 9.0],
    [6.0, 8.0, 7.0]
])

# Adiciona 1 ponto a cada nota
novas_notas = notas + 1
print(novas_notas)
```

Outro exemplo: adicionar valores diferentes a cada coluna.

```python
notas = np.array([[7, 8, 9], [6, 7, 8]])
bonus = np.array([1, 0, 0.5])

resultado = notas + bonus
print(resultado)
```

---

## Funções estatísticas

Considere o array:

```python
notas = np.array([6, 8, 7, 10, 9])
```

### Soma, produto, mínimo e máximo

```python
print(np.sum(notas))  # Soma
print(np.prod(notas)) # Produto
print(np.min(notas))  # Menor valor
print(np.max(notas))  # Maior valor
```

Também é possível usar métodos do próprio array:

```python
print(notas.sum())
print(notas.min())
print(notas.max())
```

### Média, mediana e percentis

```python
print(np.mean(notas))   # Média aritmética
print(np.median(notas)) # Mediana
print(np.percentile(notas, 25)) # Primeiro quartil
print(np.percentile(notas, 75)) # Terceiro quartil
```

### Desvio padrão e variância

```python
print(np.std(notas)) # Desvio padrão
print(np.var(notas)) # Variância
```

Por padrão, o NumPy calcula a variância populacional. Para a variância amostral, use `ddof=1`:

```python
variancia_amostral = np.var(notas, ddof=1)
desvio_amostral = np.std(notas, ddof=1)
```

### Estatísticas por eixo

```python
dados = np.array([
    [10, 20, 30],
    [40, 50, 60]
])

print(np.sum(dados))
print(np.sum(dados, axis=0)) # Soma por coluna: [50 70 90]
print(np.sum(dados, axis=1)) # Soma por linha: [60 150]
print(np.mean(dados, axis=0))
```

- `axis=0`: percorre as linhas e calcula por coluna.
- `axis=1`: percorre as colunas e calcula por linha.

### `argmin()` e `argmax()`

Retornam o índice do menor e do maior valor.

```python
valores = np.array([12, 5, 20, 8])

print(np.argmin(valores)) # 1
print(np.argmax(valores)) # 2
```

### `cumsum()` e `cumprod()`

Calculam soma e produto acumulados.

```python
valores = np.array([1, 2, 3, 4])

print(np.cumsum(valores))  # [1 3 6 10]
print(np.cumprod(valores)) # [1 2 6 24]
```

---

## Álgebra linear

### Produto elemento a elemento e produto matricial

```python
A = np.array([[1, 2], [3, 4]])
B = np.array([[5, 6], [7, 8]])

produto_elemento = A * B
produto_matricial = A @ B

print(produto_elemento)
print(produto_matricial)
```

Também é possível usar:

```python
produto_matricial = np.matmul(A, B)
```

### `dot()`

Calcula produto escalar, produto de matrizes ou produtos compatíveis entre arrays.

```python
vetor_a = np.array([1, 2, 3])
vetor_b = np.array([4, 5, 6])

print(np.dot(vetor_a, vetor_b))
# 32: 1*4 + 2*5 + 3*6
```

### Determinante e inversa

```python
A = np.array([[4, 7], [2, 6]])

print(np.linalg.det(A))
print(np.linalg.inv(A))
```

Uma matriz precisa ser quadrada e invertível para possuir inversa.

### Resolver sistemas lineares

Para resolver `Ax = b`:

```python
A = np.array([[2, 1], [1, 3]])
b = np.array([8, 13])

x = np.linalg.solve(A, b)
print(x)
```

### Norma de um vetor

```python
vetor = np.array([3, 4])
print(np.linalg.norm(vetor))
# 5.0
```

---

## Números aleatórios

A recomendação atual é criar um gerador com `np.random.default_rng()`.

```python
import numpy as np

rng = np.random.default_rng(seed=42)
```

Usar uma `seed` torna os resultados reproduzíveis.

### `random()`

Gera números decimais entre 0 e 1.

```python
print(rng.random())
print(rng.random(5))
print(rng.random((2, 3)))
```

### `integers()`

Gera números inteiros aleatórios. O limite superior é exclusivo.

```python
print(rng.integers(1, 7))
print(rng.integers(1, 7, size=10))
```

O exemplo acima simula lançamentos de um dado.

### `normal()`

Gera valores de uma distribuição normal.

```python
amostra = rng.normal(loc=100, scale=15, size=1000)
```

- `loc`: média.
- `scale`: desvio padrão.
- `size`: quantidade ou formato da amostra.

### `choice()`

Escolhe elementos aleatoriamente.

```python
opcoes = np.array(["A", "B", "C"])
escolhas = rng.choice(opcoes, size=5)
print(escolhas)
```

Também pode escolher sem repetição:

```python
amostra = rng.choice(np.arange(10), size=4, replace=False)
```

### Embaralhar com `shuffle()`

```python
valores = np.array([1, 2, 3, 4, 5])
rng.shuffle(valores)
print(valores)
```

O array é alterado no próprio lugar.

---

## Valores especiais e condições

### `np.nan`

Representa um valor ausente ou não numérico.

```python
valores = np.array([10, 20, np.nan, 40])

print(np.isnan(valores))
# [False False  True False]
```

Funções comuns para ignorar valores `NaN`:

```python
print(np.nanmean(valores))
print(np.nanmedian(valores))
print(np.nansum(valores))
print(np.nanstd(valores))
```

### `np.inf`

Representa infinito positivo ou negativo.

```python
valores = np.array([1, np.inf, -np.inf, 5])

print(np.isinf(valores))
print(np.isfinite(valores))
```

### `all()` e `any()`

```python
valores = np.array([True, True, False])

print(np.all(valores)) # False: todos precisam ser True
print(np.any(valores)) # True: pelo menos um é True
```

### `clip()`

Limita valores a um intervalo.

```python
valores = np.array([-10, 2, 8, 15])
limitados = np.clip(valores, 0, 10)

print(limitados)
# [ 0  2  8 10]
```

### `unique()`

Retorna valores únicos.

```python
valores = np.array([1, 2, 2, 3, 3, 3])
print(np.unique(valores))
# [1 2 3]
```

Com contagem:

```python
unicos, contagens = np.unique(valores, return_counts=True)
print(unicos)
print(contagens)
```

---

## Erros comuns relacionados ao NumPy

### Confundir `*` com `@`

```python
A * B  # Multiplicação elemento por elemento
A @ B  # Multiplicação matricial
```

### Tentar fazer `reshape()` com quantidade incompatível

```python
valores = np.arange(10)

# Erro: 10 elementos não podem formar uma matriz 3x3
# valores.reshape(3, 3)
```

### Esquecer que o limite final de `arange()` não é incluído

```python
np.arange(1, 5)
# [1 2 3 4]
```

Para controlar a quantidade de pontos e incluir o final, use `linspace()`:

```python
np.linspace(1, 5, 5)
# [1. 2. 3. 4. 5.]
```

### Divisão inteira ou conversão inesperada

```python
valores = np.array([1.9, 2.9])
print(valores.astype(int))
# [1 2]
```

A conversão não arredonda; ela remove a parte decimal.

### Problemas com `NaN`

```python
valores = np.array([10, 20, np.nan])
print(np.mean(valores))
# nan

print(np.nanmean(valores))
# 15.0
```

## Tabela de consulta rápida

### NumPy

| Função | Uso |
|---|---|
| `np.array()` | Cria um array |
| `np.zeros()` | Array preenchido com zeros |
| `np.ones()` | Array preenchido com uns |
| `np.full()` | Array preenchido com um valor |
| `np.arange()` | Valores espaçados por passo |
| `np.linspace()` | Quantidade definida de valores espaçados |
| `np.reshape()` | Altera o formato |
| `.T` | Transposta |
| `np.concatenate()` | Concatena arrays |
| `np.vstack()` | Empilha verticalmente |
| `np.hstack()` | Empilha horizontalmente |
| `np.split()` | Divide um array |
| `np.where()` | Seleção condicional |
| `np.unique()` | Valores únicos |
| `np.sum()` | Soma |
| `np.mean()` | Média |
| `np.median()` | Mediana |
| `np.std()` | Desvio padrão |
| `np.var()` | Variância |
| `np.min()` / `np.max()` | Menor / maior valor |
| `np.argmin()` / `np.argmax()` | Índice do menor / maior valor |
| `np.sqrt()` | Raiz quadrada |
| `np.round()` | Arredondamento |
| `np.floor()` / `np.ceil()` | Arredonda para baixo / cima |
| `np.abs()` | Valor absoluto |
| `np.sin()` / `np.cos()` | Funções trigonométricas |
| `np.isnan()` | Identifica `NaN` |
| `np.isfinite()` | Identifica valores finitos |
| `np.clip()` | Limita valores a um intervalo |
| `np.dot()` | Produto escalar ou matricial |
| `np.linalg.solve()` | Resolve sistema linear |
| `np.linalg.inv()` | Calcula matriz inversa |
| `np.random.default_rng()` | Cria gerador aleatório |

## Boas práticas

1. Use nomes de variáveis descritivos, como `notas`, `meses` e `vendas`.
2. Use `np.random.default_rng(seed)` quando precisar reproduzir resultados.
3. Evite loops quando uma operação vetorizada do NumPy resolver o problema de forma mais simples e rápida.
4. Verifique `shape` e `dtype` quando um cálculo produzir resultados inesperados.
5. Trate valores ausentes com funções como `np.nanmean()`, `np.nanmedian()` e `np.nansum()`.

## Referências oficiais

- [Documentação do NumPy](https://numpy.org/doc/)
- [Referência da API do NumPy](https://numpy.org/doc/stable/reference/)

