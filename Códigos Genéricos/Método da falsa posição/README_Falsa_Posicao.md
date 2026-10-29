# Método da Falsa Posição

## Descrição

Este projeto implementa o **Método da Falsa Posição** (*regula falsi*), um método numérico iterativo utilizado para encontrar raízes de funções não lineares em um intervalo no qual a função muda de sinal.

Assim como o método da bisseção, ele mantém um intervalo que contém a raiz. A diferença está na escolha do ponto de divisão: em vez do ponto médio, utiliza-se o ponto em que a **reta secante** que passa pelos extremos do intervalo cruza o eixo $x$. Essa escolha aproveita os valores da função e, em geral, acelera a convergência.

## Fundamentação matemática

Considere uma função contínua $f(x)$ e um intervalo $[a,b]$ tal que:

$$
f(a)\cdot f(b)<0
$$

Então existe pelo menos um valor $x^*\in(a,b)$ tal que $f(x^*)=0$.

### 1. Reta que passa pelos extremos

A reta que passa pelos pontos $(a,f(a))$ e $(b,f(b))$ é:

$$
y=f(a)+\frac{f(b)-f(a)}{b-a}\,(x-a)
$$

### 2. Ponto de interseção com o eixo $x$

Impondo $y=0$ e isolando $x$, obtém-se a nova aproximação:

$$
x_r=\frac{a\,f(b)-b\,f(a)}{f(b)-f(a)}
$$

Essa é a expressão utilizada no código.

### 3. Escolha do novo intervalo

Calcula-se $f(x_r)$ e compara-se seu sinal com o de $f(a)$:

- se $f(a)\cdot f(x_r)<0$, a raiz está em $[a,x_r]$ e faz-se $b\leftarrow x_r$;
- caso contrário, a raiz está em $[x_r,b]$ e faz-se $a\leftarrow x_r$.

### 4. Critérios de parada

O processo termina quando pelo menos uma das condições abaixo é satisfeita:

$$
|f(x_r)|<\text{tol}
$$

$$
|x_r-x_{r,\text{ant}}|<\text{tol}\quad(\text{a partir da segunda iteração})
$$

Diferentemente da bisseção, o comprimento do intervalo $b-a$ não é usado como critério de parada, pois, na falsa posição, um dos extremos pode permanecer praticamente fixo durante muitas iterações.

## Código implementado

```python
def falsa_posicao(f, a, b, tol=1e-6, max_it=100):
    fa, fb = f(a), f(b)
    if abs(fa) < tol:
        return a, 0
    if abs(fb) < tol:
        return b, 0
    if fa * fb > 0:
        raise ValueError(f"O intervalo [{a}, {b}] não apresenta mudança de sinal.")

    xr_ant = a
    for it in range(1, max_it + 1):
        if abs(fb - fa) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: f(b) e f(a) são muito próximos.")

        xr = (a * fb - b * fa) / (fb - fa)
        fxr = f(xr)
        erro = abs(xr - xr_ant) if it > 1 else abs(b - a)

        if abs(fxr) < tol or (it > 1 and erro < tol):
            return xr, it

        if fa * fxr < 0:
            b, fb = xr, fxr
        else:
            a, fa = xr, fxr

        xr_ant = xr

    return xr, max_it
```

## Explicação do código passo a passo

### 1. Definição da função

```python
def falsa_posicao(f, a, b, tol=1e-6, max_it=100):
```

Define a função que implementa o Método da Falsa Posição.

Parâmetros:

- `f`: função matemática cuja raiz se deseja encontrar;
- `a`, `b`: extremos do intervalo inicial;
- `tol`: tolerância utilizada nos critérios de parada, com valor padrão de $10^{-6}$;
- `max_it`: número máximo de iterações, com valor padrão de 100.

A função retorna uma tupla contendo a aproximação da raiz e o número de iterações realizadas.

### 2. Avaliação da função nos extremos

```python
fa, fb = f(a), f(b)
```

Calcula os valores da função nos extremos do intervalo. Esses valores são usados tanto na fórmula da nova aproximação quanto na verificação de sinais.

### 3. Verificação dos extremos

```python
if abs(fa) < tol:
    return a, 0
if abs(fb) < tol:
    return b, 0
```

Se algum dos extremos já satisfaz a tolerância, ele é retornado como raiz aproximada, sem realizar iterações.

### 4. Verificação da mudança de sinal

```python
if fa * fb > 0:
    raise ValueError(f"O intervalo [{a}, {b}] não apresenta mudança de sinal.")
```

Se $f(a)$ e $f(b)$ têm o mesmo sinal, não há garantia de raiz no intervalo e a função lança uma exceção `ValueError`.

### 5. Inicialização da aproximação anterior

```python
xr_ant = a
```

Cria a variável que guarda a aproximação da iteração anterior, usada no cálculo do erro entre aproximações sucessivas.

### 6. Laço iterativo

```python
for it in range(1, max_it + 1):
```

Inicia o processo iterativo, que pode executar até `max_it` vezes.

### 7. Verificação do denominador

```python
if abs(fb - fa) < 1e-15:
    raise ZeroDivisionError("Divisão por zero: f(b) e f(a) são muito próximos.")
```

A fórmula de $x_r$ divide por $f(b)-f(a)$. Se esse valor for nulo ou muito pequeno, a função interrompe a execução para evitar uma divisão numericamente inadequada.

### 8. Cálculo da nova aproximação

```python
xr = (a * fb - b * fa) / (fb - fa)
fxr = f(xr)
```

Calcula o ponto em que a reta que passa pelos extremos cruza o eixo $x$ e o valor da função nesse ponto.

### 9. Cálculo do erro

```python
erro = abs(xr - xr_ant) if it > 1 else abs(b - a)
```

Na primeira iteração, o erro é o comprimento do intervalo inicial. A partir da segunda, é a distância entre aproximações sucessivas:

$$
\text{erro}=|x_r-x_{r,\text{ant}}|
$$

### 10. Critério de convergência

```python
if abs(fxr) < tol or (it > 1 and erro < tol):
    return xr, it
```

A função retorna `xr` e o número de iterações se:

- o valor absoluto da função na nova aproximação é menor que a tolerância; ou
- a partir da segunda iteração, a distância entre aproximações sucessivas é menor que a tolerância.

### 11. Atualização do intervalo

```python
if fa * fxr < 0:
    b, fb = xr, fxr
else:
    a, fa = xr, fxr
```

Mantém o subintervalo em que há mudança de sinal, substituindo o extremo correspondente por `xr`.

### 12. Atualização da aproximação anterior

```python
xr_ant = xr
```

Guarda a aproximação atual para o cálculo do erro na próxima iteração.

### 13. Retorno após atingir o limite de iterações

```python
return xr, max_it
```

Se o número máximo de iterações for atingido sem satisfazer os critérios de parada, a função retorna a última aproximação e `max_it`.

**Importante:** esse retorno não indica explicitamente se houve convergência.

## Exemplo de utilização

Considere a função:

$$
f(x)=x^3-x-2
$$

Como $f(1)=-2$ e $f(2)=4$, há mudança de sinal em $[1,2]$:

```python
def f(x):
    return x**3 - x - 2

raiz, iteracoes = falsa_posicao(f, 1.0, 2.0, tol=1e-6, max_it=100)

print("Raiz aproximada:", raiz)
print("Número de iterações:", iteracoes)
print("f(raiz):", f(raiz))
```

Resultado obtido na execução:

```text
Raiz aproximada: 1.5213794617901568
Número de iterações: 12
f(raiz): -1.4563135963818752e-06
```

A raiz real dessa função está aproximadamente em:

$$
x\approx 1.52138
$$

## Cuidados e limitações

- **Mudança de sinal obrigatória:** o intervalo inicial precisa ter $f(a)\cdot f(b)<0$. Caso contrário, a função lança `ValueError`.
- **Extremo que permanece fixo:** em funções com forte curvatura, um dos extremos pode ficar praticamente parado, e a convergência pode ser lenta. Nesses casos, a bisseção pode ser mais eficiente.
- **Intervalo não encolhe necessariamente:** como o intervalo pode não reduzir de tamanho de forma significativa, o código não usa $b-a$ como critério de parada.
- **Divisão por zero:** o código verifica se $f(b)-f(a)$ é muito pequeno, mas isso não elimina todos os problemas numéricos possíveis.
- **Continuidade:** o método pressupõe que a função seja contínua no intervalo.
- **Critérios de parada alternativos:** basta uma das condições ser satisfeita. Um valor pequeno de $|f(x_r)|$ não garante, isoladamente, que $x_r$ esteja próximo da raiz.
- **Limite de iterações:** se o limite for atingido, a função retorna a última aproximação sem indicar explicitamente se houve convergência.

## Complexidade computacional

Em cada iteração, o método realiza um número constante de operações aritméticas e uma avaliação da função.

Se o custo de avaliar $f(x)$ for considerado constante, o custo computacional é aproximadamente:

$$
O(N)
$$

em que $N$ representa o número de iterações realizadas. O armazenamento adicional é constante:

$$
O(1)
$$

## Resumo

A função `falsa_posicao(f, a, b, tol=1e-6, max_it=100)` implementa o Método da Falsa Posição para encontrar raízes aproximadas de funções contínuas em um intervalo com mudança de sinal.

O método combina a segurança do intervalo da bisseção com a interpolação linear da secante. Em muitos casos converge mais rápido que a bisseção, mas pode ser lento quando um dos extremos permanece fixo. É importante escolher bem o intervalo e verificar o valor da função na raiz encontrada.
