# Método da Bisseção

## Descrição

Este projeto implementa o **Método da Bisseção**, um método numérico iterativo utilizado para encontrar raízes de funções não lineares em um intervalo no qual a função muda de sinal.

A cada iteração, o intervalo é dividido ao meio e é mantida a metade em que ocorre a mudança de sinal. Dessa forma, o intervalo que contém a raiz é reduzido progressivamente até que a aproximação atinja a tolerância especificada.

O método é simples e robusto: se a função for contínua e o intervalo inicial tiver mudança de sinal, a convergência é garantida. Em contrapartida, sua convergência é lenta quando comparada a métodos como Newton-Raphson ou secante.

## Fundamentação matemática

Considere uma função contínua $f(x)$ e um intervalo $[a,b]$ tal que:

$$
f(a)\cdot f(b)<0
$$

Pelo **Teorema do Valor Intermediário**, existe pelo menos um valor $x^*\in(a,b)$ tal que:

$$
f(x^*)=0
$$

### 1. Ponto médio do intervalo

Em cada iteração, calcula-se o ponto médio:

$$
x_r=\frac{a+b}{2}
$$

### 2. Escolha do novo intervalo

Calcula-se $f(x_r)$ e compara-se seu sinal com o de $f(a)$:

- se $f(a)\cdot f(x_r)<0$, a raiz está em $[a,x_r]$ e faz-se $b\leftarrow x_r$;
- caso contrário, a raiz está em $[x_r,b]$ e faz-se $a\leftarrow x_r$.

### 3. Redução do intervalo

A cada iteração, o comprimento do intervalo é reduzido pela metade. Após $k$ iterações, o comprimento é:

$$
b_k-a_k=\frac{b_0-a_0}{2^k}
$$

Como a raiz pertence ao intervalo, o erro da aproximação $x_r$ satisfaz:

$$
|x_r-x^*|\le\frac{b_k-a_k}{2}
$$

### 4. Número de iterações necessário

Para garantir que o erro seja menor que uma tolerância $\text{tol}$, é suficiente que:

$$
k\ge\log_2\!\left(\frac{b_0-a_0}{\text{tol}}\right)
$$

### 5. Critérios de parada

O processo termina quando pelo menos uma das condições abaixo é satisfeita:

$$
|f(x_r)|<\text{tol}
$$

$$
\frac{b-a}{2}<\text{tol}
$$

$$
|x_r-x_{r,\text{ant}}|<\text{tol}\quad(\text{a partir da segunda iteração})
$$

## Código implementado

```python
def bissecao(f, a, b, tol=1e-6, max_it=100):
    fa, fb = f(a), f(b)
    if abs(fa) < tol:
        return a, 0
    if abs(fb) < tol:
        return b, 0
    if fa * fb > 0:
        raise ValueError(f"O intervalo [{a}, {b}] não apresenta mudança de sinal.")

    xr_ant = a
    for it in range(1, max_it + 1):
        xr = (a + b) / 2.0
        fxr = f(xr)
        erro = abs(xr - xr_ant) if it > 1 else (b - a) / 2.0

        if abs(fxr) < tol or (b - a) / 2.0 < tol or (it > 1 and erro < tol):
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
def bissecao(f, a, b, tol=1e-6, max_it=100):
```

Define a função que implementa o Método da Bisseção.

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

Calcula os valores da função nos extremos do intervalo.

Esses valores são usados para verificar a mudança de sinal e, depois, para decidir qual metade do intervalo deve ser mantida.

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

Se $f(a)$ e $f(b)$ têm o mesmo sinal, o produto é positivo e não há garantia de raiz no intervalo. Nesse caso, a função lança uma exceção `ValueError`.

### 5. Inicialização da aproximação anterior

```python
xr_ant = a
```

Cria a variável que guarda a aproximação da iteração anterior, usada no cálculo do erro entre aproximações sucessivas.

### 6. Laço iterativo

```python
for it in range(1, max_it + 1):
```

Inicia o processo iterativo, que pode executar até `max_it` vezes. A variável `it` registra o número da iteração atual.

### 7. Cálculo do ponto médio

```python
xr = (a + b) / 2.0
fxr = f(xr)
```

Calcula o ponto médio do intervalo atual e o valor da função nesse ponto.

### 8. Cálculo do erro

```python
erro = abs(xr - xr_ant) if it > 1 else (b - a) / 2.0
```

Na primeira iteração, ainda não há aproximação anterior; por isso, o erro é estimado pela metade do comprimento do intervalo:

$$
\text{erro}=\frac{b-a}{2}
$$

A partir da segunda iteração, o erro é a distância entre aproximações sucessivas:

$$
\text{erro}=|x_r-x_{r,\text{ant}}|
$$

### 9. Critério de convergência

```python
if abs(fxr) < tol or (b - a) / 2.0 < tol or (it > 1 and erro < tol):
    return xr, it
```

A função retorna `xr` e o número de iterações se pelo menos uma das condições for satisfeita:

- o valor absoluto da função no ponto médio é menor que a tolerância;
- a metade do comprimento do intervalo é menor que a tolerância;
- a partir da segunda iteração, a distância entre aproximações sucessivas é menor que a tolerância.

### 10. Atualização do intervalo

```python
if fa * fxr < 0:
    b, fb = xr, fxr
else:
    a, fa = xr, fxr
```

Se $f(a)$ e $f(x_r)$ têm sinais opostos, a raiz está na metade esquerda e o extremo direito passa a ser `xr`. Caso contrário, a raiz está na metade direita e o extremo esquerdo passa a ser `xr`.

### 11. Atualização da aproximação anterior

```python
xr_ant = xr
```

Guarda a aproximação atual para o cálculo do erro na próxima iteração.

### 12. Retorno após atingir o limite de iterações

```python
return xr, max_it
```

Se o número máximo de iterações for atingido sem satisfazer os critérios de parada, a função retorna a última aproximação e `max_it`.

**Importante:** esse retorno não indica explicitamente se houve convergência. O programa que chama a função deve considerar essa possibilidade.

## Exemplo de utilização

Considere a função:

$$
f(x)=x^3-x-2
$$

Como $f(1)=-2$ e $f(2)=4$, há mudança de sinal em $[1,2]$:

```python
def f(x):
    return x**3 - x - 2

raiz, iteracoes = bissecao(f, 1.0, 2.0, tol=1e-6, max_it=100)

print("Raiz aproximada:", raiz)
print("Número de iterações:", iteracoes)
print("f(raiz):", f(raiz))
```

Resultado obtido na execução:

```text
Raiz aproximada: 1.5213804244995117
Número de iterações: 20
f(raiz): 4.265829404825894e-06
```

A raiz real dessa função está aproximadamente em:

$$
x\approx 1.52138
$$

## Cuidados e limitações

- **Mudança de sinal obrigatória:** o intervalo inicial precisa ter $f(a)\cdot f(b)<0$. Caso contrário, a função lança `ValueError`.
- **Continuidade:** o método pressupõe que a função seja contínua no intervalo. Descontinuidades, como uma assíntota vertical, podem produzir mudança de sinal sem que exista raiz.
- **Raízes de multiplicidade par:** em raízes em que a função toca o eixo sem mudar de sinal, como $f(x)=x^2$, o método não consegue isolar a raiz por mudança de sinal.
- **Várias raízes no intervalo:** se houver mais de uma raiz, o método converge para uma delas, e não necessariamente para a desejada.
- **Convergência lenta:** o erro é reduzido pela metade a cada iteração (convergência linear), exigindo mais iterações que métodos mais rápidos.
- **Critérios de parada alternativos:** basta uma das condições ser satisfeita. Um valor pequeno de $|f(x_r)|$ não garante, isoladamente, que $x_r$ esteja próximo da raiz, especialmente quando a função é muito achatada.
- **Limite de iterações:** se o limite for atingido, a função retorna a última aproximação sem indicar explicitamente se houve convergência.

## Complexidade computacional

Em cada iteração, o método realiza um número constante de operações aritméticas e uma avaliação da função.

O número de iterações necessário para atingir a tolerância é da ordem de:

$$
O\!\left(\log_2\frac{b-a}{\text{tol}}\right)
$$

Se o custo de avaliar $f(x)$ for considerado constante, esse também é o custo total. O armazenamento adicional é constante:

$$
O(1)
$$

## Resumo

A função `bissecao(f, a, b, tol=1e-6, max_it=100)` implementa o Método da Bisseção para encontrar raízes aproximadas de funções contínuas em um intervalo com mudança de sinal.

O método divide o intervalo ao meio a cada iteração e mantém a metade que contém a raiz. É simples e tem convergência garantida sob as hipóteses corretas, mas é mais lento que outros métodos. Para obter resultados confiáveis, é importante escolher um intervalo adequado e verificar o valor da função na raiz encontrada.
