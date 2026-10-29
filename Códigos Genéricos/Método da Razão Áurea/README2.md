# Método da Razão Áurea

Implementação em Python do método da razão áurea para encontrar aproximadamente um mínimo ou um máximo de uma função real em um intervalo fechado $[a,b]$. O algoritmo reduz progressivamente o intervalo de busca, reutilizando avaliações da função para diminuir o número de cálculos necessários.

Este README explica a teoria matemática, a lógica do algoritmo e cada parte do código, além de apresentar exemplos de uso, critérios de parada e limitações.

---

## Sumário

- [1. Objetivo](#1-objetivo)
- [2. O que é otimização unidimensional](#2-o-que-é-otimização-unidimensional)
- [3. Fundamento matemático da razão áurea](#3-fundamento-matemático-da-razão-áurea)
- [4. Como o algoritmo funciona](#4-como-o-algoritmo-funciona)
- [5. Explicação do código linha por linha](#5-explicação-do-código-linha-por-linha)
- [6. Como utilizar a função](#6-como-utilizar-a-função)
- [7. Exemplos completos](#7-exemplos-completos)
- [8. Critérios de parada e interpretação do retorno](#8-critérios-de-parada-e-interpretação-do-retorno)
- [9. Limitações e cuidados](#9-limitações-e-cuidados)
- [10. Complexidade computacional](#10-complexidade-computacional)
- [11. Possíveis melhorias](#11-possíveis-melhorias)
- [12. Referências](#12-referências)

---

## 1. Objetivo

A função implementada no código é:

```python
razao_aurea(f, a, b, tol=1e-6, max_it=100, minimizar=True)
```

Ela procura um ponto de mínimo ou máximo de uma função $f(x)$ dentro do intervalo $[a,b]$.

Parâmetros:

- `f`: função matemática a ser otimizada;
- `a`: limite inferior do intervalo;
- `b`: limite superior do intervalo;
- `tol`: tolerância usada no critério de parada;
- `max_it`: número máximo de iterações;
- `minimizar`: indica se a busca é por mínimo (`True`) ou máximo (`False`).

A função retorna uma tupla com dois elementos:

```python
(ponto_aproximado, numero_de_iteracoes)
```

O primeiro valor é o ponto médio do intervalo final, usado como aproximação do ótimo. O segundo indica quantas iterações foram realizadas.

**Importante:** o método é uma técnica de otimização numérica. Ele não calcula necessariamente uma raiz de $f(x)$; procura um extremo da função.

## 2. O que é otimização unidimensional

Otimização unidimensional consiste em encontrar o valor de $x$ que minimiza ou maximiza uma função de uma variável, dentro de um intervalo.

Para minimização, o objetivo é encontrar

$$
x^*=\arg\min_{x\in[a,b]} f(x).
$$

Para maximização, o objetivo é encontrar

$$
x^*=\arg\max_{x\in[a,b]} f(x).
$$

O valor $x^*$ é o ponto ótimo. O valor $f(x^*)$ é o valor da função nesse ponto.

Exemplos de aplicações:
- minimizar o erro de um modelo;
- encontrar o menor custo de uma operação;
- maximizar o rendimento de um processo;
- ajustar um parâmetro de um modelo matemático.

O método da razão áurea é um método de busca por intervalo: em cada iteração, compara valores da função em dois pontos internos e descarta uma parte do intervalo que não precisa mais ser considerada, sob a hipótese de unimodalidade.

## 3. Fundamento matemático da razão áurea

### 3.1 Proporção áurea

A razão áurea é associada ao número

$$
\varphi=\frac{1+\sqrt{5}}{2}\approx1{,}618034.
$$

O código utiliza a constante complementar

$$
\rho=\frac{\sqrt{5}-1}{2}\approx0{,}618034.
$$

As duas constantes se relacionam por $\rho=1/\varphi$, e também vale $1-\rho\approx0{,}381966$.

A constante `phi` no código corresponde a $\rho$, e não ao número $\varphi\approx1{,}618$. Essa escolha permite calcular os pontos internos do intervalo com fórmulas compactas.

### 3.2 Escolha dos pontos internos

Dado um intervalo atual $[a,b]$, o algoritmo define:

$$
x_1=b-\rho(b-a),
$$

$$
x_2=a+\rho(b-a).
$$

Esses pontos ficam dentro do intervalo e dividem-no em segmentos relacionados pela razão áurea.

Em geral, $a<x_1<x_2<b$, desde que $a<b$.

O algoritmo avalia $f(x_1)$ e $f(x_2)$. Comparando os resultados, decide qual lado do intervalo pode ser descartado.

### 3.3 Hipótese de unimodalidade

A garantia clássica do método depende de a função ser **unimodal** no intervalo considerado: ela deve apresentar um único mínimo (no caso de minimização) ou um único máximo (no caso de maximização) no intervalo.

Se houver vários extremos locais, a comparação de apenas dois pontos pode descartar a região que contém o melhor extremo global. Portanto, escolher um intervalo adequado é uma parte importante do método.

### 3.4 Por que a razão áurea é útil?

O posicionamento dos pontos permite que, após descartar um dos lados, um dos pontos já avaliados seja reutilizado na iteração seguinte. Assim, em condições normais, o algoritmo precisa avaliar a função em apenas um novo ponto por iteração, além das duas avaliações iniciais.

Isso é especialmente vantajoso quando calcular $f(x)$ é caro.

## 4. Como o algoritmo funciona

Em cada iteração:

1. Existem dois pontos internos, $x_1$ e $x_2$, e seus valores $f(x_1)$ e $f(x_2)$.
2. O algoritmo compara os valores de acordo com o objetivo: mínimo ou máximo.
3. Descarta uma das extremidades do intervalo.
4. Reutiliza o ponto e o valor da função que permanecem válidos.
5. Calcula um novo ponto interno e avalia a função nele.
6. Verifica se a largura do intervalo ficou menor que a tolerância.

O processo termina quando

$$
|b-a|<\text{tol}
$$

ou quando o número de iterações chega a `max_it`.

O ponto retornado é

$$
x_{\text{aprox}}=\frac{a+b}{2},
$$

isto é, o ponto médio do intervalo restante.

### 4.1 Regra para minimização

Se $f(x_1)<f(x_2)$, a região à direita de $x_2$ pode ser descartada para uma função unimodal em minimização. O novo intervalo é $[a,x_2]$.

Caso contrário, o novo intervalo é $[x_1,b]$.

### 4.2 Regra para maximização

Para maximização, a comparação é invertida. Se $f(x_1)>f(x_2)$, a região à direita de $x_2$ pode ser descartada. Caso contrário, o novo intervalo é $[x_1,b]$.

Em caso de igualdade, o código entra no ramo `else`; isso mantém a região $[x_1,b]$. O comportamento é válido como regra de redução para o caso unimodal, embora empates e funções com trechos planos possam reduzir a informação disponível na comparação.

## 5. Explicação do código linha por linha

### 5.1 Definição da função

```python
def razao_aurea(f, a, b, tol=1e-6, max_it=100, minimizar=True):
```

- `def` cria uma função reutilizável.
- `f` é a função objetivo, passada como argumento.
- `a` e `b` delimitam o intervalo de busca.
- `tol=1e-6` é a tolerância padrão, equivalente a $10^{-6}$.
- `max_it=100` limita o algoritmo a 100 iterações por padrão.
- `minimizar=True` seleciona minimização por padrão. Para maximização, deve-se passar `minimizar=False`.

O código pressupõe um intervalo válido, normalmente com `a < b`, e uma função `f` que possa ser avaliada nos pontos internos.

### 5.2 Constante da razão áurea

```python
phi = (5.0 ** 0.5 - 1.0) / 2.0
```

- `5.0 ** 0.5` calcula a raiz quadrada de 5.
- Subtrai-se 1 e divide-se por 2.
- O resultado é aproximadamente $0{,}618034$.

Embora a variável se chame `phi`, ela representa a constante complementar $(\sqrt5-1)/2$, não $(1+\sqrt5)/2$.

### 5.3 Cálculo dos pontos internos

```python
x1 = b - phi * (b - a)
x2 = a + phi * (b - a)
fx1, fx2 = f(x1), f(x2)
```

- `b - a` calcula a largura do intervalo.
- `x1` e `x2` são os pontos internos da razão áurea.
- `f(x1)` e `f(x2)` avaliam a função nesses pontos.
- A atribuição múltipla guarda os dois resultados em `fx1` e `fx2`.

As avaliações iniciais são feitas uma única vez. O reaproveitamento posterior evita recalcular um valor que já é conhecido.

### 5.4 Laço de iterações

```python
for it in range(1, max_it + 1):
```

O laço começa em 1 e vai até `max_it`, inclusive. A variável `it` registra a iteração atual.

### 5.5 Escolha do lado a descartar

```python
if (fx1 < fx2) if minimizar else (fx1 > fx2):
```

Esta é uma expressão condicional em Python. Ela equivale conceitualmente a:

```python
if fx1 < fx2:       # quando minimizar é True
    ...
else:               # condição alternativa usada quando minimizar é False
    ...
```

Mais precisamente, a condição é `fx1 < fx2` quando `minimizar` é verdadeiro e `fx1 > fx2` quando `minimizar` é falso. Assim, o mesmo código pode executar minimização ou maximização.

### 5.6 Primeiro ramo: atualizar o limite superior

```python
b = x2
x2, fx2 = x1, fx1
x1 = b - phi * (b - a)
fx1 = f(x1)
```

Quando a comparação entra neste ramo:

1. `b = x2` reduz o limite superior.
2. `x2, fx2 = x1, fx1` reutiliza o ponto antigo `x1` e seu valor como o novo ponto `x2`.
3. Calcula-se a nova posição de `x1`.
4. Avalia-se `f` apenas nesse novo ponto.

A atribuição múltipla de Python troca os valores de forma simultânea, usando os valores anteriores à atribuição.

### 5.7 Segundo ramo: atualizar o limite inferior

```python
else:
    a = x1
    x1, fx1 = x2, fx2
    x2 = a + phi * (b - a)
    fx2 = f(x2)
```

Quando a comparação não entra no primeiro ramo:

1. `a = x1` aumenta o limite inferior.
2. O ponto e o valor antigos de `x2` passam a ser o novo `x1` e `fx1`.
3. Calcula-se o novo `x2`.
4. Avalia-se a função no novo ponto.

Novamente, um dos pontos anteriores é reaproveitado e somente uma nova avaliação é necessária.

### 5.8 Critério de parada

```python
if abs(b - a) < tol:
    return (a + b) / 2.0, it
```

- `b - a` é a largura atual do intervalo.
- `abs(...)` garante que a largura seja considerada em valor absoluto.
- Se ela for menor que `tol`, a função retorna o ponto médio e o número da iteração em que parou.

O teste é estrito (`<`), não `<=`.

### 5.9 Limite máximo de iterações

```python
return (a + b) / 2.0, max_it
```

Se nenhuma iteração satisfizer o critério de tolerância, a função retorna o ponto médio do intervalo final e `max_it`.

**Importante:** o retorno não informa explicitamente se a convergência ocorreu ou se o limite de iterações foi atingido. Para saber isso, o usuário precisa comparar o intervalo final ou adaptar a função para retornar também um indicador de convergência.

## 6. Como utilizar a função

### Requisitos

- Python 3.
- Nenhuma biblioteca externa é necessária para executar o algoritmo.

Salve a função em um arquivo chamado, por exemplo, `razao_aurea.py`. Se usar outro nome, ajuste a importação.

### Exemplo de chamada

```python
from razao_aurea import razao_aurea

def f(x):
    return (x - 2)**2 + 1

ponto, iteracoes = razao_aurea(f, 0, 5)

print("Ponto de mínimo:", ponto)
print("Número de iterações:", iteracoes)
print("Valor da função:", f(ponto))
```

A função $f(x)=(x-2)^2+1$ tem mínimo em $x=2$, com $f(2)=1$. O método deve retornar uma aproximação de 2, dependendo da tolerância e do número de iterações.

## 7. Exemplos completos

### 7.1 Minimizar uma função

Considere

$$
f(x)=(x-2)^2+1.
$$

A função é uma parábola voltada para cima e seu mínimo ocorre em $x=2$.

```python
from razao_aurea import razao_aurea

def f(x):
    return (x - 2)**2 + 1

x_min, it = razao_aurea(
    f,
    a=0,
    b=5,
    tol=1e-6,
    max_it=100
)

print(f"x mínimo aproximado = {x_min:.8f}")
print(f"f(x) = {f(x_min):.8f}")
print(f"Iterações = {it}")
```

### 7.2 Maximizar uma função

Considere

$$
g(x)=-(x-3)^2+10.
$$

Esta parábola é voltada para baixo e seu máximo ocorre em $x=3$, com $g(3)=10$.

```python
from razao_aurea import razao_aurea

def g(x):
    return -(x - 3)**2 + 10

x_max, it = razao_aurea(
    g,
    a=0,
    b=6,
    tol=1e-6,
    max_it=100,
    minimizar=False
)

print(f"x máximo aproximado = {x_max:.8f}")
print(f"g(x) = {g(x_max):.8f}")
print(f"Iterações = {it}")
```

### 7.3 Interpretando o retorno

O retorno é uma tupla de dois elementos. Por exemplo:

```python
ponto, iteracoes = razao_aurea(f, 0, 5)
```

- `ponto` contém o ponto médio do intervalo final.
- `iteracoes` contém a iteração em que o algoritmo parou, ou `max_it` caso tenha alcançado o limite sem satisfazer a tolerância.

O algoritmo retorna a posição aproximada do extremo, não o valor ótimo da função diretamente. Para obter o valor da função no ponto encontrado, use `f(ponto)`.

## 8. Critérios de parada e interpretação do retorno

O algoritmo para em uma de duas situações:

1. **Convergência por intervalo:** quando $|b-a|<\text{tol}$.
2. **Limite de iterações:** quando o laço chega a `max_it`.

A tolerância controla a precisão desejada na localização do ponto, medida pela largura do intervalo, não diretamente pelo erro no valor $f(x)$.

Se o intervalo final tiver largura $w$, o ponto médio fica a no máximo $w/2$ de qualquer ponto dentro desse intervalo. A interpretação como erro em relação ao ótimo depende de o ótimo estar mantido dentro do intervalo e das hipóteses do método.

A tolerância deve ser escolhida considerando a escala do problema. Uma tolerância muito pequena pode exigir mais iterações e, em cálculos de ponto flutuante, não garante precisão arbitrária.

## 9. Limitações e cuidados

1. **Unimodalidade:** o método pressupõe um único extremo relevante no intervalo. Funções com vários mínimos ou máximos podem levar a um resultado que não é o ótimo global.
2. **Intervalo válido:** o código pressupõe normalmente `a < b`. Não há validação explícita dessa condição.
3. **Função avaliável:** `f` deve aceitar números reais e retornar valores comparáveis. Valores `NaN`, infinitos ou exceções dentro de `f` não são tratados pelo algoritmo.
4. **Máximo de iterações:** se `max_it` for insuficiente, a função retorna uma aproximação sem indicar explicitamente que a tolerância não foi atingida.
5. **Validação de parâmetros:** o código não verifica se `tol` é positiva nem se `max_it` é um inteiro positivo.
6. **Ótimo na fronteira:** a busca por razão áurea é usualmente aplicada a um intervalo que contém o extremo interior. Se o melhor valor estiver na fronteira, o método pode aproximar essa fronteira, mas convém verificar os extremos $a$ e $b$ também.
7. **Empates e regiões planas:** se os valores avaliados forem iguais, a regra do `else` descarta um dos lados; em funções com trechos planos, a localização do extremo pode ser menos informativa.

## 10. Complexidade computacional

Cada iteração realiza uma quantidade constante de operações aritméticas, comparações e, normalmente, uma nova avaliação da função. O custo real depende bastante do custo de calcular $f(x)$.

A cada iteração, o comprimento do intervalo é reduzido aproximadamente pelo fator

$$
\rho=\frac{\sqrt5-1}{2}\approx0{,}618034.
$$

Após $k$ reduções, a largura é aproximadamente

$$
L_k=L_0\rho^k,
$$

em que $L_0=b-a$ é a largura inicial.

Para que a largura fique menor que `tol`, a estimativa idealizada é

$$
k>\frac{\ln(\text{tol}/L_0)}{\ln(\rho)}.
$$

Como $0<\rho<1$, o denominador é negativo; a desigualdade deve ser interpretada com cuidado, mas a fórmula permite estimar quantas reduções são necessárias. Na prática, `max_it` deve ser suficiente para a tolerância e a largura inicial escolhidas.

## 11. Possíveis melhorias

- validar `a < b`, `tol > 0` e `max_it` antes de iniciar;
- retornar também um indicador explícito de convergência;
- retornar o intervalo final e o valor da função no ponto estimado;
- registrar o histórico de intervalos, pontos e avaliações para análise ou gráficos;
- adicionar testes automatizados para minimização, maximização e funções com extremos na fronteira;
- permitir contabilizar o número de avaliações da função;
- documentar a função com *docstring* e anotações de tipo;
- comparar os resultados com uma solução analítica quando disponível.

## 12. Referências

Tópicos recomendados para estudo:

- Otimização unidimensional.
- Método da busca pela razão áurea.
- Unimodalidade de funções.
- Critérios de parada e erro numérico.

Documentação oficial do Python:

- https://docs.python.org/3/

---

## Licença

Este README não define uma licença para o código. Caso o repositório seja público e você queira permitir reutilização, adicione um arquivo de licença e escolha os termos adequados ao projeto.
