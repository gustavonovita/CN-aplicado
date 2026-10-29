# Método de Newton-Raphson

## Descrição

Este projeto implementa o **Método de Newton-Raphson**, um método numérico iterativo utilizado para encontrar raízes de funções não lineares.

O método parte de uma única aproximação inicial, $x_0$, e utiliza a **derivada** da função para construir a reta tangente ao gráfico no ponto atual. O ponto em que essa tangente cruza o eixo $x$ é a nova aproximação da raiz.

Quando a aproximação inicial é adequada, a convergência é muito rápida (quadrática). Em contrapartida, o método exige o cálculo da derivada e não garante convergência para qualquer ponto inicial.

## Fundamentação matemática

Considere uma função diferenciável:

$$
f(x)=0
$$

O objetivo é encontrar um valor $x$ para o qual a função seja nula ou suficientemente próxima de zero.

### 1. Reta tangente

Em torno de um ponto $x_n$, a função pode ser aproximada pela reta tangente:

$$
y=f(x_n)+f'(x_n)\,(x-x_n)
$$

### 2. Fórmula de iteração

Impondo $y=0$ e isolando $x$, obtém-se a nova aproximação:

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
$$

Essa é a expressão utilizada no código.

### 3. Interpretação geométrica

A cada iteração, traça-se a tangente à curva em $(x_n,f(x_n))$. A nova aproximação $x_{n+1}$ é o ponto em que essa tangente intercepta o eixo $x$.

### 4. Critérios de parada

O processo termina quando pelo menos uma das condições abaixo é satisfeita:

$$
|x_{n+1}-x_n|<\text{tol}
$$

$$
|f(x_{n+1})|<\text{tol}
$$

### 5. Ordem de convergência

Para raízes simples e aproximação inicial suficientemente próxima da raiz, a convergência é quadrática:

$$
|e_{n+1}|\approx C\,|e_n|^2
$$

em que $e_n=x_n-x^*$ é o erro na iteração $n$. Na prática, o número de dígitos corretos tende a dobrar a cada iteração.

## Código implementado

```python
def newton_raphson(f, df, x0, tol=1e-6, max_it=100):
    x_atual = x0
    fx = f(x_atual)
    if abs(fx) < tol:
        return x_atual, 0

    for it in range(1, max_it + 1):
        dfx = df(x_atual)
        if abs(dfx) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: derivada f'(x) muito próxima de zero.")

        x_novo = x_atual - fx / dfx
        fx_novo = f(x_novo)
        erro = abs(x_novo - x_atual)

        if erro < tol or abs(fx_novo) < tol:
            return x_novo, it

        x_atual, fx = x_novo, fx_novo

    return x_atual, max_it
```

## Explicação do código passo a passo

### 1. Definição da função

```python
def newton_raphson(f, df, x0, tol=1e-6, max_it=100):
```

Define a função que implementa o Método de Newton-Raphson.

Parâmetros:

- `f`: função matemática cuja raiz se deseja encontrar;
- `df`: derivada de `f`, fornecida pelo usuário;
- `x0`: aproximação inicial;
- `tol`: tolerância utilizada nos critérios de parada, com valor padrão de $10^{-6}$;
- `max_it`: número máximo de iterações, com valor padrão de 100.

A função retorna uma tupla contendo a aproximação da raiz e o número de iterações realizadas.

### 2. Avaliação da função no ponto inicial

```python
x_atual = x0
fx = f(x_atual)
```

Define o ponto atual e calcula o valor da função nele. Esse valor é guardado em `fx` e reaproveitado na iteração seguinte.

### 3. Verificação da aproximação inicial

```python
if abs(fx) < tol:
    return x_atual, 0
```

Se o ponto inicial já satisfaz a tolerância, ele é retornado como raiz aproximada, sem realizar iterações.

### 4. Laço iterativo

```python
for it in range(1, max_it + 1):
```

Inicia o processo iterativo, que pode executar até `max_it` vezes. A variável `it` registra o número da iteração atual.

### 5. Cálculo da derivada

```python
dfx = df(x_atual)
```

Calcula o valor da derivada no ponto atual, que corresponde à inclinação da reta tangente.

### 6. Verificação da derivada

```python
if abs(dfx) < 1e-15:
    raise ZeroDivisionError("Divisão por zero: derivada f'(x) muito próxima de zero.")
```

A fórmula divide por $f'(x_n)$. Se a derivada for nula ou muito pequena, a tangente é praticamente horizontal e não cruza o eixo $x$ de forma útil. Nesse caso, a função lança uma exceção `ZeroDivisionError`.

### 7. Cálculo da nova aproximação

```python
x_novo = x_atual - fx / dfx
```

Aplica a fórmula de Newton-Raphson:

$$
x_{\text{novo}}=x_{\text{atual}}-\frac{f(x_{\text{atual}})}{f'(x_{\text{atual}})}
$$

### 8. Avaliação da função e cálculo do erro

```python
fx_novo = f(x_novo)
erro = abs(x_novo - x_atual)
```

Calcula o valor da função na nova aproximação e o erro, definido como a distância entre aproximações sucessivas:

$$
\text{erro}=|x_{\text{novo}}-x_{\text{atual}}|
$$

### 9. Critério de convergência

```python
if erro < tol or abs(fx_novo) < tol:
    return x_novo, it
```

A função retorna `x_novo` e o número de iterações se pelo menos uma das condições for satisfeita:

- a distância entre aproximações sucessivas é menor que a tolerância;
- o valor absoluto da função na nova aproximação é menor que a tolerância.

### 10. Atualização do ponto

```python
x_atual, fx = x_novo, fx_novo
```

Se o método ainda não convergiu, a nova aproximação passa a ser o ponto atual, junto com o valor da função já calculado.

### 11. Retorno após atingir o limite de iterações

```python
return x_atual, max_it
```

Se o número máximo de iterações for atingido sem satisfazer os critérios de parada, a função retorna a última aproximação e `max_it`.

**Importante:** esse retorno não indica explicitamente se houve convergência.

## Exemplo de utilização

Considere a função e sua derivada:

$$
f(x)=x^3-x-2,\qquad f'(x)=3x^2-1
$$

Utilizando a aproximação inicial:

$$
x_0=1.5
$$

O código de utilização é:

```python
def f(x):
    return x**3 - x - 2

def df(x):
    return 3 * x**2 - 1

raiz, iteracoes = newton_raphson(f, df, 1.5, tol=1e-6, max_it=100)

print("Raiz aproximada:", raiz)
print("Número de iterações:", iteracoes)
print("f(raiz):", f(raiz))
```

Resultado obtido na execução:

```text
Raiz aproximada: 1.5213798059647863
Número de iterações: 2
f(raiz): 5.893874259754739e-07
```

A raiz real dessa função está aproximadamente em:

$$
x\approx 1.52138
$$

## Cuidados e limitações

- **Necessidade da derivada:** o usuário precisa fornecer `df` corretamente. Uma derivada errada leva a resultados incorretos sem aviso.
- **Derivada nula ou muito pequena:** se $f'(x_n)\approx 0$, a função lança `ZeroDivisionError`. Pontos de máximo ou mínimo local são especialmente problemáticos.
- **Dependência do ponto inicial:** uma aproximação inicial distante da raiz pode levar à divergência, a oscilações ou à convergência para outra raiz.
- **Ausência de intervalo:** diferentemente da bisseção e da falsa posição, não há garantia de que a raiz permaneça em um intervalo conhecido.
- **Raízes múltiplas:** em raízes de multiplicidade maior que 1, a convergência deixa de ser quadrática e torna-se mais lenta.
- **Critérios de parada alternativos:** basta uma das condições ser satisfeita. Um deslocamento pequeno não garante, isoladamente, que o valor da função também seja pequeno.
- **Limite de iterações:** se o limite for atingido, a função retorna a última aproximação sem indicar explicitamente se houve convergência.

## Complexidade computacional

Em cada iteração, o método realiza um número constante de operações aritméticas, uma avaliação de $f$ na nova aproximação e uma avaliação de $f'$ no ponto atual.

Se o custo de avaliar $f(x)$ e $f'(x)$ for considerado constante, o custo computacional é aproximadamente:

$$
O(N)
$$

em que $N$ representa o número de iterações realizadas. O armazenamento adicional é constante:

$$
O(1)
$$

## Resumo

A função `newton_raphson(f, df, x0, tol=1e-6, max_it=100)` implementa o Método de Newton-Raphson para encontrar raízes aproximadas de funções diferenciáveis.

O método utiliza a derivada para construir a reta tangente e converge rapidamente quando o ponto inicial é adequado. Por outro lado, exige a derivada, não garante convergência e pode falhar quando $f'(x)$ se aproxima de zero. É importante escolher bem o ponto inicial e verificar o valor da função na raiz encontrada.
