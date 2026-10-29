# Método de Müller

## Descrição

Este projeto implementa o **Método de Müller**, um método numérico iterativo utilizado para encontrar raízes de funções não lineares.

O método utiliza três aproximações iniciais, $x_0$, $x_1$ e $x_2$, para construir um polinômio interpolador de segundo grau. Em seguida, calcula uma das raízes desse polinômio para obter uma nova aproximação, $x_3$, da raiz da função original.

O procedimento é repetido até que o erro entre aproximações sucessivas ou o valor absoluto da função seja inferior à tolerância especificada.

Diferentemente de métodos que utilizam apenas aproximações lineares, como o método da secante, o Método de Müller utiliza uma interpolação quadrática, o que permite representar a função localmente por uma parábola.

## Fundamentação matemática

Considere uma função não linear:

$$
f(x)=0
$$

O objetivo é encontrar um valor $x$ para o qual a função seja nula ou suficientemente próxima de zero.

### 1. Aproximações iniciais

O método começa com três pontos:

$$
x_0,\quad x_1,\quad x_2
$$

Calculam-se os valores da função nesses pontos:

$$
f(x_0),\quad f(x_1),\quad f(x_2)
$$

Esses valores são utilizados para construir um polinômio quadrático interpolador.

### 2. Cálculo das diferenças entre os pontos

Calculam-se as distâncias entre as aproximações:

$$
h_0=x_1-x_0
$$

$$
h_1=x_2-x_1
$$

Essas diferenças são utilizadas no cálculo das diferenças divididas.

### 3. Cálculo das diferenças divididas

As diferenças divididas de primeira ordem são:

$$
d_0=\frac{f(x_1)-f(x_0)}{h_0}
$$

$$
d_1=\frac{f(x_2)-f(x_1)}{h_1}
$$

Elas representam inclinações aproximadas da função entre os pontos.

### 4. Construção do polinômio quadrático

O coeficiente quadrático é calculado por:

$$
a=\frac{d_1-d_0}{h_1+h_0}
$$

O coeficiente linear, na forma deslocada utilizada pelo código, é:

$$
b=ah_1+d_1
$$

O termo constante é:

$$
c=f(x_2)
$$

Assim, o polinômio interpolador em relação ao deslocamento $dx=x-x_2$ é:

$$
p(dx)=a(dx)^2+b(dx)+c
$$

Para encontrar a próxima aproximação, procura-se uma raiz desse polinômio.

### 5. Resolução da equação quadrática

A equação a ser resolvida é:

$$
a(dx)^2+b(dx)+c=0
$$

Pela fórmula quadrática:

$$
dx=\frac{-b\pm\sqrt{b^2-4ac}}{2a}
$$

O código utiliza uma forma equivalente, racionalizada:

$$
dx=\frac{-2c}{b\pm\sqrt{b^2-4ac}}
$$

Essa forma pode reduzir problemas de cancelamento numérico. O código escolhe o denominador de maior valor absoluto entre as duas alternativas.

### 6. Atualização da aproximação

Depois de calcular $dx$, obtém-se:

$$
x_3=x_2+dx
$$

Calcula-se então:

$$
f(x_3)
$$

A nova aproximação é aceita como solução quando:

$$
|dx|<\text{tol}
$$

ou

$$
|f(x_3)|<\text{tol}
$$

Se nenhum desses critérios for satisfeito, as aproximações são atualizadas e o processo continua.

## Código implementado

```python
def muller(f, x0, x1, x2, tol=1e-6, max_it=100):
    fx0, fx1, fx2 = f(x0), f(x1), f(x2)
    if abs(fx0) < tol:
        return x0, 0
    if abs(fx1) < tol:
        return x1, 0
    if abs(fx2) < tol:
        return x2, 0

    for it in range(1, max_it + 1):
        h0 = x1 - x0
        h1 = x2 - x1
        if abs(h0) < 1e-15 or abs(h1) < 1e-15 or abs(h1 + h0) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: pontos muito próximos.")

        d0 = (fx1 - fx0) / h0
        d1 = (fx2 - fx1) / h1

        a = (d1 - d0) / (h1 + h0)
        b = a * h1 + d1
        c = fx2

        raiz_disc = (b**2 - 4 * a * c) ** 0.5
        den_pos = b + raiz_disc
        den_neg = b - raiz_disc
        denominador = den_pos if abs(den_pos) >= abs(den_neg) else den_neg

        if abs(denominador) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: denominador nulo na fórmula quadrática.")

        dx = -2 * c / denominador
        x3 = x2 + dx
        fx3 = f(x3)
        erro = abs(dx)

        if erro < tol or abs(fx3) < tol:
            return x3, it

        x0, fx0 = x1, fx1
        x1, fx1 = x2, fx2
        x2, fx2 = x3, fx3

    return x2, max_it
```

## Explicação do código passo a passo

### 1. Definição da função

```python
def muller(f, x0, x1, x2, tol=1e-6, max_it=100):
```

Define a função que implementa o Método de Müller.

Parâmetros:

- `f`: função matemática cuja raiz se deseja encontrar;
- `x0`, `x1`, `x2`: três aproximações iniciais;
- `tol`: tolerância utilizada para verificar a convergência, com valor padrão de $10^{-6}$;
- `max_it`: número máximo de iterações, com valor padrão de 100.

A função retorna uma tupla contendo a aproximação da raiz e o número de iterações realizadas.

### 2. Avaliação da função nos pontos iniciais

```python
fx0, fx1, fx2 = f(x0), f(x1), f(x2)
```

Calcula os valores da função nos três pontos iniciais.

Esses valores são necessários para determinar as diferenças divididas e construir o polinômio quadrático interpolador.

### 3. Verificação da primeira aproximação

```python
if abs(fx0) < tol:
    return x0, 0
```

Verifica se o primeiro ponto já satisfaz a tolerância. Caso satisfaça, retorna `x0` como raiz aproximada, sem realizar iterações.

### 4. Verificação da segunda aproximação

```python
if abs(fx1) < tol:
    return x1, 0
```

Realiza a mesma verificação para `x1`.

### 5. Verificação da terceira aproximação

```python
if abs(fx2) < tol:
    return x2, 0
```

Verifica se `x2` já é uma solução suficientemente precisa.

### 6. Laço iterativo

```python
for it in range(1, max_it + 1):
```

Inicia o processo iterativo, que pode executar até `max_it` vezes.

A variável `it` registra o número da iteração atual.

### 7. Cálculo dos intervalos

```python
h0 = x1 - x0
h1 = x2 - x1
```

Calcula as diferenças entre os pontos consecutivos.

Esses valores são utilizados no cálculo das diferenças divididas e do polinômio interpolador.

### 8. Verificação de pontos muito próximos

```python
if abs(h0) < 1e-15 or abs(h1) < 1e-15 or abs(h1 + h0) < 1e-15:
    raise ZeroDivisionError("Divisão por zero: pontos muito próximos.")
```

Verifica se alguma diferença utilizada nos denominadores é nula ou extremamente pequena.

Se isso ocorrer, a função interrompe a execução e lança uma exceção para evitar uma divisão numericamente inadequada.

### 9. Cálculo das diferenças divididas

```python
d0 = (fx1 - fx0) / h0
d1 = (fx2 - fx1) / h1
```

Calcula as diferenças divididas de primeira ordem.

Elas representam inclinações aproximadas da função nos intervalos entre os pontos.

### 10. Cálculo do coeficiente quadrático

```python
a = (d1 - d0) / (h1 + h0)
```

Calcula o coeficiente quadrático do polinômio interpolador.

Esse coeficiente determina a curvatura da parábola construída a partir dos três pontos.

### 11. Cálculo do coeficiente linear

```python
b = a * h1 + d1
```

Calcula o coeficiente linear na representação do polinômio em torno de `x2`.

### 12. Cálculo do termo constante

```python
c = fx2
```

Define o termo constante como o valor da função no ponto atual `x2`.

### 13. Cálculo da raiz do discriminante

```python
raiz_disc = (b**2 - 4 * a * c) ** 0.5
```

Calcula a raiz quadrada do discriminante da equação quadrática:

$$
\Delta=b^2-4ac
$$

O valor do discriminante determina a natureza das raízes do polinômio.

**Atenção:** na implementação apresentada, o cálculo é feito com a potência `0.5`. Se o discriminante for negativo, o resultado pode se tornar um número complexo. A função não inclui tratamento explícito para esse caso.

### 14. Cálculo dos denominadores candidatos

```python
den_pos = b + raiz_disc
den_neg = b - raiz_disc
```

Calcula as duas possibilidades para o denominador da fórmula quadrática racionalizada.

### 15. Escolha do denominador

```python
denominador = den_pos if abs(den_pos) >= abs(den_neg) else den_neg
```

Seleciona o denominador com maior valor absoluto.

Essa escolha busca reduzir o cancelamento numérico que pode ocorrer quando dois números de magnitudes semelhantes são subtraídos.

### 16. Verificação do denominador

```python
if abs(denominador) < 1e-15:
    raise ZeroDivisionError(
        "Divisão por zero: denominador nulo na fórmula quadrática."
    )
```

Verifica se o denominador escolhido é nulo ou muito pequeno.

Se for, a execução é interrompida para evitar uma divisão numericamente inadequada.

### 17. Cálculo do deslocamento

```python
dx = -2 * c / denominador
```

Calcula o deslocamento necessário para obter a raiz do polinômio quadrático interpolador.

A expressão corresponde à forma racionalizada da fórmula quadrática.

### 18. Atualização da aproximação

```python
x3 = x2 + dx
```

Calcula a nova aproximação da raiz somando o deslocamento ao ponto atual.

### 19. Avaliação da função na nova aproximação

```python
fx3 = f(x3)
```

Calcula o valor da função original na nova aproximação.

Esse valor será utilizado para verificar se a solução encontrada atende ao critério de convergência.

### 20. Cálculo do erro

```python
erro = abs(dx)
```

Define o erro como o módulo do deslocamento entre a aproximação atual e a nova aproximação:

$$
\text{erro}=|x_3-x_2|=|dx|
$$

### 21. Critério de convergência

```python
if erro < tol or abs(fx3) < tol:
    return x3, it
```

Verifica se pelo menos uma das condições de convergência foi satisfeita:

- O deslocamento é menor que a tolerância;
- O valor absoluto da função na nova aproximação é menor que a tolerância.

Se alguma condição for satisfeita, a função retorna a raiz aproximada e o número de iterações realizadas.

### 22. Atualização dos pontos

```python
x0, fx0 = x1, fx1
x1, fx1 = x2, fx2
x2, fx2 = x3, fx3
```

Se o método ainda não convergiu, atualiza os pontos para a próxima iteração.

O ponto mais antigo é descartado, e a nova aproximação passa a ser o ponto atual.

### 23. Retorno após atingir o limite de iterações

```python
return x2, max_it
```

Se o método atingir o número máximo de iterações sem satisfazer os critérios de convergência, retorna a última aproximação e `max_it`.

**Importante:** esse retorno não garante que a raiz tenha convergido. O programa que chama a função deve considerar a possibilidade de não convergência.

## Exemplo de utilização

Considere a função:

$$
f(x)=x^3-x-2
$$

Deseja-se encontrar uma raiz real utilizando três aproximações iniciais:

$$
x_0=1,\qquad x_1=1.5,\qquad x_2=2
$$

O código de utilização é:

```python
def f(x):
    return x**3 - x - 2

raiz, iteracoes = muller(
    f,
    x0=1.0,
    x1=1.5,
    x2=2.0,
    tol=1e-6,
    max_it=100
)

print("Raiz aproximada:", raiz)
print("Número de iterações:", iteracoes)
print("f(raiz):", f(raiz))
```

A raiz real dessa função está aproximadamente em:

$$
x\approx 1.52138
$$

A quantidade exata de iterações depende dos critérios de parada e das aproximações calculadas durante a execução.

## Cuidados e limitações

- **Escolha dos pontos iniciais:** os três pontos influenciam a parábola interpoladora e o comportamento do método.
- **Ausência de garantia de convergência:** o método pode não convergir para determinadas funções ou escolhas de pontos iniciais.
- **Divisão por zero:** o código verifica diferenças e denominadores muito pequenos, mas essas verificações não eliminam todos os possíveis problemas numéricos.
- **Discriminante negativo:** a implementação não trata explicitamente o caso em que $b^2-4ac<0$. Como o código utiliza a potência `0.5`, esse caso pode produzir valores complexos.
- **Raízes complexas:** o Método de Müller pode ser adaptado para trabalhar com raízes complexas, mas essa implementação não foi preparada para tratar todos os casos de forma explícita.
- **Critérios de parada alternativos:** como o código utiliza `erro < tol or abs(fx3) < tol`, basta uma das condições ser satisfeita. Isso significa que um deslocamento pequeno não garante, isoladamente, que o valor da função também seja pequeno.
- **Limite de iterações:** se o método não convergir, a função retorna a última aproximação sem indicar explicitamente se houve convergência.

## Complexidade computacional

Em cada iteração, o método realiza um número constante de operações aritméticas e calcula a função uma vez na nova aproximação.

Se o custo de avaliar $f(x)$ for considerado constante, o custo computacional é aproximadamente:

$$
O(N)
$$

em que $N$ representa o número de iterações realizadas.

Se a avaliação da função tiver custo computacional significativo, esse custo deverá ser considerado separadamente.

O armazenamento adicional utilizado pela implementação é constante, ou seja:

$$
O(1)
$$

## Resumo

A função `muller(f, x0, x1, x2, tol=1e-6, max_it=100)` implementa o Método de Müller para encontrar raízes aproximadas de funções não lineares.

O método utiliza três pontos iniciais, constrói um polinômio quadrático interpolador e calcula uma nova aproximação a cada iteração. A implementação inclui verificações para diferenças e denominadores muito pequenos, além de critérios de parada baseados no deslocamento e no valor da função.

Para obter resultados confiáveis, é importante escolher adequadamente os pontos iniciais, verificar o valor da função na raiz encontrada e considerar a possibilidade de não convergência ou de ocorrência de valores complexos.
