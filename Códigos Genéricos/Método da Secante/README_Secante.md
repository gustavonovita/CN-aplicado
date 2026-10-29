# Método da Secante

## Descrição

Este projeto implementa o **Método da Secante**, um método numérico iterativo utilizado para encontrar raízes de funções não lineares.

O método utiliza duas aproximações iniciais, $x_0$ e $x_1$, e substitui a derivada do método de Newton-Raphson por uma aproximação obtida a partir da reta secante que passa pelos dois pontos. Assim, não é necessário conhecer a derivada analítica da função.

Diferentemente da bisseção e da falsa posição, o método da secante **não mantém um intervalo** com mudança de sinal. Isso o torna mais rápido em geral, mas sem garantia de convergência.

## Fundamentação matemática

Considere uma função não linear:

$$
f(x)=0
$$

O objetivo é encontrar um valor $x$ para o qual a função seja nula ou suficientemente próxima de zero.

### 1. Origem do método

O método de Newton-Raphson utiliza a derivada:

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
$$

No método da secante, a derivada é aproximada pela inclinação da reta que passa pelos dois últimos pontos:

$$
f'(x_n)\approx\frac{f(x_n)-f(x_{n-1})}{x_n-x_{n-1}}
$$

### 2. Fórmula de iteração

Substituindo a aproximação da derivada, obtém-se:

$$
x_{n+1}=x_n-f(x_n)\,\frac{x_n-x_{n-1}}{f(x_n)-f(x_{n-1})}
$$

Essa é a expressão utilizada no código, com $x_{n-1}=x_0$ e $x_n=x_1$.

### 3. Interpretação geométrica

A nova aproximação $x_{n+1}$ é o ponto em que a reta que passa por $(x_{n-1},f(x_{n-1}))$ e $(x_n,f(x_n))$ cruza o eixo $x$.

### 4. Atualização dos pontos

Após calcular o novo ponto, descarta-se o mais antigo:

$$
x_{n-1}\leftarrow x_n,\qquad x_n\leftarrow x_{n+1}
$$

### 5. Critérios de parada

O processo termina quando pelo menos uma das condições abaixo é satisfeita:

$$
|x_{n+1}-x_n|<\text{tol}
$$

$$
|f(x_{n+1})|<\text{tol}
$$

### 6. Ordem de convergência

Para raízes simples e boas aproximações iniciais, a ordem de convergência do método da secante é:

$$
p=\frac{1+\sqrt{5}}{2}\approx 1.618
$$

Ou seja, é mais rápida que a bisseção (linear), porém um pouco mais lenta que Newton-Raphson (quadrática), com a vantagem de não exigir a derivada.

## Código implementado

```python
def secante(f, x0, x1, tol=1e-6, max_it=100):
    fx0, fx1 = f(x0), f(x1)
    if abs(fx0) < tol:
        return x0, 0
    if abs(fx1) < tol:
        return x1, 0

    for it in range(1, max_it + 1):
        if abs(fx1 - fx0) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: f(x1) e f(x0) são muito próximos.")

        x_novo = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        fx_novo = f(x_novo)
        erro = abs(x_novo - x1)

        if erro < tol or abs(fx_novo) < tol:
            return x_novo, it

        x0, fx0 = x1, fx1
        x1, fx1 = x_novo, fx_novo

    return x1, max_it
```

## Explicação do código passo a passo

### 1. Definição da função

```python
def secante(f, x0, x1, tol=1e-6, max_it=100):
```

Define a função que implementa o Método da Secante.

Parâmetros:

- `f`: função matemática cuja raiz se deseja encontrar;
- `x0`, `x1`: duas aproximações iniciais;
- `tol`: tolerância utilizada nos critérios de parada, com valor padrão de $10^{-6}$;
- `max_it`: número máximo de iterações, com valor padrão de 100.

A função retorna uma tupla contendo a aproximação da raiz e o número de iterações realizadas.

### 2. Avaliação da função nos pontos iniciais

```python
fx0, fx1 = f(x0), f(x1)
```

Calcula os valores da função nos dois pontos iniciais. Esses valores são reaproveitados nas iterações seguintes, de modo que cada iteração exige apenas uma nova avaliação de $f$.

### 3. Verificação das aproximações iniciais

```python
if abs(fx0) < tol:
    return x0, 0
if abs(fx1) < tol:
    return x1, 0
```

Se algum dos pontos iniciais já satisfaz a tolerância, ele é retornado como raiz aproximada, sem realizar iterações.

### 4. Laço iterativo

```python
for it in range(1, max_it + 1):
```

Inicia o processo iterativo, que pode executar até `max_it` vezes. A variável `it` registra o número da iteração atual.

### 5. Verificação do denominador

```python
if abs(fx1 - fx0) < 1e-15:
    raise ZeroDivisionError("Divisão por zero: f(x1) e f(x0) são muito próximos.")
```

A fórmula divide por $f(x_1)-f(x_0)$. Se esse valor for nulo ou muito pequeno, a reta secante é praticamente horizontal e a função interrompe a execução com uma exceção `ZeroDivisionError`.

### 6. Cálculo da nova aproximação

```python
x_novo = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
```

Aplica a fórmula da secante:

$$
x_{\text{novo}}=x_1-f(x_1)\,\frac{x_1-x_0}{f(x_1)-f(x_0)}
$$

### 7. Avaliação da função e cálculo do erro

```python
fx_novo = f(x_novo)
erro = abs(x_novo - x1)
```

Calcula o valor da função na nova aproximação e o erro, definido como a distância entre aproximações sucessivas:

$$
\text{erro}=|x_{\text{novo}}-x_1|
$$

### 8. Critério de convergência

```python
if erro < tol or abs(fx_novo) < tol:
    return x_novo, it
```

A função retorna `x_novo` e o número de iterações se pelo menos uma das condições for satisfeita:

- a distância entre aproximações sucessivas é menor que a tolerância;
- o valor absoluto da função na nova aproximação é menor que a tolerância.

### 9. Atualização dos pontos

```python
x0, fx0 = x1, fx1
x1, fx1 = x_novo, fx_novo
```

Se o método ainda não convergiu, o ponto mais antigo é descartado e a nova aproximação passa a ser o ponto atual.

### 10. Retorno após atingir o limite de iterações

```python
return x1, max_it
```

Se o número máximo de iterações for atingido sem satisfazer os critérios de parada, a função retorna a última aproximação e `max_it`.

**Importante:** esse retorno não indica explicitamente se houve convergência.

## Exemplo de utilização

Considere a função:

$$
f(x)=x^3-x-2
$$

Utilizando as aproximações iniciais:

$$
x_0=1,\qquad x_1=2
$$

O código de utilização é:

```python
def f(x):
    return x**3 - x - 2

raiz, iteracoes = secante(f, 1.0, 2.0, tol=1e-6, max_it=100)

print("Raiz aproximada:", raiz)
print("Número de iterações:", iteracoes)
print("f(raiz):", f(raiz))
```

Resultado obtido na execução:

```text
Raiz aproximada: 1.5213797079848717
Número de iterações: 6
f(raiz): 7.015477976324291e-09
```

A raiz real dessa função está aproximadamente em:

$$
x\approx 1.52138
$$

## Cuidados e limitações

- **Ausência de garantia de convergência:** como não há intervalo com mudança de sinal, o método pode divergir ou convergir para uma raiz diferente da desejada.
- **Escolha dos pontos iniciais:** os dois pontos devem estar razoavelmente próximos da raiz. Pontos mal escolhidos podem levar a divergência.
- **Pontos com mesmo valor de função:** se $f(x_0)\approx f(x_1)$, a reta secante é quase horizontal e a função lança `ZeroDivisionError`.
- **Raízes múltiplas:** em raízes de multiplicidade maior que 1, a convergência é mais lenta.
- **Critérios de parada alternativos:** basta uma das condições ser satisfeita. Um deslocamento pequeno não garante, isoladamente, que o valor da função também seja pequeno.
- **Limite de iterações:** se o limite for atingido, a função retorna a última aproximação sem indicar explicitamente se houve convergência.

## Complexidade computacional

Em cada iteração, o método realiza um número constante de operações aritméticas e **uma única** avaliação da função, pois os valores anteriores são reaproveitados.

Se o custo de avaliar $f(x)$ for considerado constante, o custo computacional é aproximadamente:

$$
O(N)
$$

em que $N$ representa o número de iterações realizadas. O armazenamento adicional é constante:

$$
O(1)
$$

## Resumo

A função `secante(f, x0, x1, tol=1e-6, max_it=100)` implementa o Método da Secante para encontrar raízes aproximadas de funções não lineares.

O método dispensa o cálculo da derivada e converge de forma superlinear para raízes simples, mas não garante convergência. Para obter resultados confiáveis, é importante escolher bem os pontos iniciais, verificar o valor da função na raiz encontrada e considerar a possibilidade de não convergência.
