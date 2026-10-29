# Método da Interpolação Quadrática para Otimização de Funções

## 1. Objetivo

Este projeto implementa o **método da interpolação quadrática**, uma técnica numérica utilizada para encontrar pontos de mínimo ou máximo de uma função.

A ideia central do método é aproximar uma função por um polinômio quadrático construído a partir de três pontos conhecidos. Em seguida, calcula-se o vértice desse polinômio, que fornece uma estimativa da posição do extremo da função original.

A função implementada no código é:

```python
def interpolacao_quadratica(
    f, x0, x1, x2, tol=1e-6, max_it=100, minimizar=True
):
```

Ela recebe uma função matemática, três pontos iniciais e parâmetros de controle. Ao final, retorna a estimativa da posição do extremo e o número de iterações realizadas.

O método pode ser utilizado para:

- Encontrar mínimos de funções.
- Encontrar máximos de funções.
- Estimar numericamente a posição de um extremo.
- Resolver problemas de otimização unidimensional.

---

## 2. Fundamentação matemática

### 2.1. Interpolação quadrática

Considere uma função \(f(x)\) conhecida em três pontos distintos:

\[
(x_0,f(x_0)),\qquad
(x_1,f(x_1)),\qquad
(x_2,f(x_2))
\]

O método constrói um polinômio quadrático da forma:

\[
p(x)=ax^2+bx+c
\]

Esse polinômio passa pelos três pontos fornecidos, ou seja:

\[
p(x_0)=f(x_0)
\]

\[
p(x_1)=f(x_1)
\]

\[
p(x_2)=f(x_2)
\]

Os coeficientes \(a\), \(b\) e \(c\) são determinados de modo que o polinômio represente a função nesses pontos.

A aproximação quadrática permite estimar o comportamento da função entre os pontos conhecidos.

### 2.2. Determinação do vértice

Para encontrar o extremo do polinômio quadrático, calcula-se sua derivada:

\[
p'(x)=2ax+b
\]

No vértice, a derivada é igual a zero:

\[
2ax+b=0
\]

Isolando \(x\):

\[
x=-\frac{b}{2a}
\]

Portanto, a posição do vértice é:

\[
\boxed{x_v=-\frac{b}{2a}}
\]

Se \(a>0\), o polinômio apresenta um mínimo.

Se \(a<0\), o polinômio apresenta um máximo.

O método da interpolação quadrática utiliza essa propriedade para estimar o extremo da função original.

**Importante:** o vértice do polinômio interpolador é uma estimativa do extremo da função. Ele não garante, por si só, que o extremo global da função original tenha sido encontrado.

### 2.3. Fórmula utilizada no código

O código calcula diretamente a posição do vértice do polinômio interpolador, sem calcular explicitamente os coeficientes \(a\), \(b\) e \(c\).

A fórmula utilizada é:

\[
x_3=
\frac{
f(x_0)(x_1^2-x_2^2)
+f(x_1)(x_2^2-x_0^2)
+f(x_2)(x_0^2-x_1^2)
}{
2\left[
f(x_0)(x_1-x_2)
+f(x_1)(x_2-x_0)
+f(x_2)(x_0-x_1)
\right]
}
\]

O resultado \(x_3\) representa a estimativa da posição do vértice do polinômio quadrático.

Depois de calcular \(x_3\), o algoritmo avalia a função nesse ponto:

\[
f(x_3)
\]

Essa avaliação é utilizada para decidir como atualizar os três pontos de referência para a próxima iteração.

---

## 3. Parâmetros da função

A função possui sete parâmetros:

```python
interpolacao_quadratica(
    f, x0, x1, x2, tol=1e-6, max_it=100, minimizar=True
)
```

| Parâmetro | Descrição |
|---|---|
| `f` | Função matemática que será otimizada. |
| `x0` | Primeiro ponto inicial. |
| `x1` | Segundo ponto inicial. |
| `x2` | Terceiro ponto inicial. |
| `tol` | Tolerância usada no critério de parada. O padrão é `1e-6`. |
| `max_it` | Número máximo de iterações. O padrão é `100`. |
| `minimizar` | Define o objetivo: `True` para minimização e `False` para maximização. |

### 3.1. Parâmetro `f`

Representa a função que será avaliada pelo algoritmo.

Por exemplo:

```python
def f(x):
    return (x - 2)**2 + 1
```

Nesse caso, a função é:

\[
f(x)=(x-2)^2+1
\]

### 3.2. Parâmetros `x0`, `x1` e `x2`

São os três pontos utilizados para construir a aproximação quadrática inicial.

É recomendável escolher pontos distintos e em uma região relevante para o problema de otimização.

A escolha dos pontos influencia a qualidade da aproximação e o comportamento da convergência.

### 3.3. Parâmetro `tol`

Define a tolerância utilizada para verificar a diferença entre duas estimativas consecutivas do extremo.

O valor padrão é:

\[
\text{tol}=10^{-6}
\]

Quanto menor a tolerância, menor será a variação permitida entre as estimativas antes de interromper o algoritmo. Isso não significa necessariamente que o resultado terá essa mesma precisão em relação ao extremo verdadeiro.

### 3.4. Parâmetro `max_it`

Define o número máximo de iterações permitidas.

O valor padrão é:

```python
max_it = 100
```

Esse limite impede que o algoritmo continue indefinidamente caso o critério de parada não seja satisfeito.

### 3.5. Parâmetro `minimizar`

Controla o tipo de otimização:

- `minimizar=True`: procura um mínimo.
- `minimizar=False`: procura um máximo.

O valor padrão é `True`, portanto a função realiza uma tentativa de minimização quando o argumento não é informado.

---

## 4. Explicação detalhada do código

### 4.1. Definição da função e avaliação inicial

```python
def interpolacao_quadratica(
    f, x0, x1, x2, tol=1e-6, max_it=100, minimizar=True
):
    fx0, fx1, fx2 = f(x0), f(x1), f(x2)
    x3_ant = x1
    x3 = x1
```

A primeira linha define a função e seus parâmetros.

Em seguida:

```python
fx0, fx1, fx2 = f(x0), f(x1), f(x2)
```

O código avalia a função nos três pontos iniciais:

\[
f(x_0),\qquad f(x_1),\qquad f(x_2)
\]

Os resultados são armazenados nas variáveis `fx0`, `fx1` e `fx2`.

Esses valores são utilizados para construir a aproximação quadrática.

Depois:

```python
x3_ant = x1
x3 = x1
```

As variáveis `x3_ant` e `x3` são inicializadas com o valor de `x1`.

- `x3`: armazenará a estimativa atual do extremo.
- `x3_ant`: armazenará a estimativa anterior para calcular o erro entre iterações.

A inicialização de `x3_ant` também permite que o código compare estimativas consecutivas nas iterações seguintes.

### 4.2. Laço principal das iterações

```python
for it in range(1, max_it + 1):
```

Esse laço executa as iterações do método.

A expressão:

```python
range(1, max_it + 1)
```

faz com que o contador `it` assuma valores de 1 até `max_it`, inclusive.

Por exemplo, se:

```python
max_it = 100
```

o algoritmo poderá realizar até 100 iterações.

Em cada iteração, o código calcula uma nova estimativa do extremo e atualiza os pontos utilizados na interpolação.

### 4.3. Cálculo do numerador

```python
num = (
    fx0 * (x1**2 - x2**2)
    + fx1 * (x2**2 - x0**2)
    + fx2 * (x0**2 - x1**2)
)
```

O numerador é calculado utilizando as coordenadas dos três pontos e os valores da função nesses pontos.

Matematicamente:

\[
N=
f(x_0)(x_1^2-x_2^2)
+f(x_1)(x_2^2-x_0^2)
+f(x_2)(x_0^2-x_1^2)
\]

O operador `**2` calcula o quadrado de um valor.

Por exemplo:

```python
x1**2
```

representa:

\[
x_1^2
\]

O numerador faz parte da expressão utilizada para determinar o vértice do polinômio interpolador.

### 4.4. Cálculo do denominador

```python
den = 2.0 * (
    fx0 * (x1 - x2)
    + fx1 * (x2 - x0)
    + fx2 * (x0 - x1)
)
```

O denominador é dado por:

\[
D=2\left[
f(x_0)(x_1-x_2)
+f(x_1)(x_2-x_0)
+f(x_2)(x_0-x_1)
\right]
\]

O fator `2.0` corresponde ao fator 2 presente na fórmula do vértice.

O numerador e o denominador são calculados separadamente para facilitar a leitura e a organização do código.

### 4.5. Verificação de divisão por zero

```python
if abs(den) < 1e-15:
    raise ZeroDivisionError(
        "Divisão por zero: os três pontos são colineares."
    )
```

Antes de calcular a nova estimativa, o código verifica se o denominador é muito próximo de zero.

A condição:

```python
abs(den) < 1e-15
```

significa que o módulo do denominador é inferior a \(10^{-15}\).

Se isso ocorrer, a função interrompe a execução e apresenta uma exceção do tipo `ZeroDivisionError`.

Essa verificação evita a divisão por um denominador considerado numericamente muito pequeno.

A mensagem informa que os três pontos são colineares. Mais precisamente, isso corresponde à situação em que os pontos da função não permitem determinar um vértice quadrático finito por meio dessa fórmula, ou em que o cálculo se torna numericamente degenerado.

### 4.6. Cálculo da nova estimativa

```python
x3 = num / den
```

Se o denominador for válido, o código calcula a nova estimativa:

\[
\boxed{x_3=\frac{N}{D}}
\]

O valor `x3` representa a posição estimada do vértice do polinômio quadrático interpolador.

Esse ponto é utilizado como candidato a extremo da função original.

### 4.7. Avaliação da função no novo ponto

```python
fx3 = f(x3)
```

Depois de calcular `x3`, o código avalia a função nesse ponto:

\[
f(x_3)
\]

O resultado é armazenado em `fx3`.

Essa avaliação é necessária para comparar o candidato com os pontos de referência e decidir como atualizar o intervalo de interpolação.

### 4.8. Cálculo do erro

```python
erro = abs(x3 - x3_ant)
```

O erro é calculado pela diferença absoluta entre a estimativa atual e a estimativa anterior:

\[
e=|x_3-x_{3,\text{ant}}|
\]

Em que:

- \(x_3\) é a estimativa atual.
- \(x_{3,\text{ant}}\) é a estimativa anterior.

Por exemplo, se a estimativa anterior era 2,1 e a atual é 2,0001:

\[
e=|2,0001-2,1|=0,0999
\]

Esse valor indica quanto a estimativa mudou entre duas iterações consecutivas.

### 4.9. Critério de parada

```python
if it > 1 and erro < tol:
    return x3, it
```

O código verifica se:

1. A iteração atual é posterior à primeira.
2. O erro é menor que a tolerância.

Matematicamente:

\[
e<\text{tol}
\]

A condição `it > 1` impede que o algoritmo encerre já na primeira iteração com base na comparação entre `x3` e seu valor inicial.

Se as duas condições forem satisfeitas, a função retorna:

```python
return x3, it
```

O retorno contém:

- `x3`: estimativa da posição do extremo.
- `it`: número de iterações realizadas.

**Importante:** esse critério considera apenas a mudança na posição estimada. Ele não verifica diretamente se a derivada da função é próxima de zero, nem garante que o extremo seja global.

### 4.10. Comparação do valor da função

```python
melhor = (fx3 < fx1) if minimizar else (fx3 > fx1)
```

Essa linha compara o valor da função no novo ponto com o valor no ponto central `x1`.

Se `minimizar=True`, a condição é:

\[
f(x_3)<f(x_1)
\]

Se `minimizar=False`, a condição é:

\[
f(x_3)>f(x_1)
\]

A variável `melhor` recebe `True` quando o novo ponto apresenta um valor de função mais favorável ao objetivo escolhido.

Assim, o algoritmo utiliza a mesma estrutura para tentar minimizar ou maximizar a função.

### 4.11. Atualização dos pontos quando `x3 > x1`

```python
if x3 > x1:
    if melhor:
        x0, fx0 = x1, fx1
        x1, fx1 = x3, fx3
    else:
        x2, fx2 = x3, fx3
```

Esse bloco verifica se o novo ponto está à direita do ponto central:

\[
x_3>x_1
\]

Se isso acontecer, o algoritmo avalia se o novo ponto é melhor para o objetivo escolhido.

**Caso 1: o novo ponto é melhor**

```python
x0, fx0 = x1, fx1
x1, fx1 = x3, fx3
```

O ponto central anterior passa para a posição de `x0`, e o novo ponto passa a ocupar a posição de `x1`.

Os valores correspondentes da função também são atualizados.

**Caso 2: o novo ponto não é melhor**

```python
x2, fx2 = x3, fx3
```

O novo ponto passa a ocupar a posição de `x2`.

A ideia é manter os pontos de referência atualizados de acordo com o resultado obtido na interpolação.

### 4.12. Atualização dos pontos quando `x3 <= x1`

```python
else:
    if melhor:
        x2, fx2 = x1, fx1
        x1, fx1 = x3, fx3
    else:
        x0, fx0 = x3, fx3
```

Esse bloco é executado quando o novo ponto está à esquerda ou coincide com o ponto central:

\[
x_3\leq x_1
\]

Novamente, o código verifica se o valor da função no novo ponto é melhor.

**Caso 1: o novo ponto é melhor**

```python
x2, fx2 = x1, fx1
x1, fx1 = x3, fx3
```

O ponto central anterior passa para a posição de `x2`, enquanto o novo ponto passa a ocupar a posição de `x1`.

**Caso 2: o novo ponto não é melhor**

```python
x0, fx0 = x3, fx3
```

O novo ponto passa a ocupar a posição de `x0`.

Essas atualizações modificam os três pontos utilizados na próxima interpolação.

A implementação pressupõe que os pontos escolhidos e suas atualizações sejam adequados ao problema. Ela não valida explicitamente a ordenação inicial dos pontos nem garante que os três pontos continuem formando uma região que contenha o extremo desejado.

### 4.13. Armazenamento da estimativa anterior

```python
x3_ant = x3
```

Ao final de cada iteração, o código atualiza `x3_ant` com a estimativa atual.

Assim, na próxima iteração, será possível calcular a diferença entre a nova estimativa e a anterior.

Essa atualização é essencial para o critério de parada baseado na variação de `x3`.

### 4.14. Retorno ao atingir o limite de iterações

```python
return x3, max_it
```

Se o algoritmo atingir o final do laço sem satisfazer o critério de parada, a função retorna a última estimativa calculada e o limite máximo de iterações.

O retorno possui dois elementos:

1. A última estimativa do extremo.
2. O número máximo de iterações.

Isso não significa necessariamente que o método convergiu. O código original não fornece um indicador separado de convergência.

---

## 5. Exemplo completo de utilização

Considere a função:

\[
f(x)=(x-2)^2+1
\]

Desenvolvendo a expressão:

\[
f(x)=x^2-4x+5
\]

Como o coeficiente de \(x^2\) é positivo, a parábola possui um mínimo.

O mínimo ocorre em:

\[
x=2
\]

Substituindo na função:

\[
f(2)=(2-2)^2+1=1
\]

Portanto, o mínimo exato é:

\[
\boxed{x_{\min}=2}
\]

\[
\boxed{f(x_{\min})=1}
\]

### 5.1. Código Python

```python
def interpolacao_quadratica(
    f, x0, x1, x2, tol=1e-6, max_it=100, minimizar=True
):
    fx0, fx1, fx2 = f(x0), f(x1), f(x2)
    x3_ant = x1
    x3 = x1

    for it in range(1, max_it + 1):
        num = (
            fx0 * (x1**2 - x2**2)
            + fx1 * (x2**2 - x0**2)
            + fx2 * (x0**2 - x1**2)
        )

        den = 2.0 * (
            fx0 * (x1 - x2)
            + fx1 * (x2 - x0)
            + fx2 * (x0 - x1)
        )

        if abs(den) < 1e-15:
            raise ZeroDivisionError(
                "Divisão por zero: os três pontos são colineares."
            )

        x3 = num / den
        fx3 = f(x3)
        erro = abs(x3 - x3_ant)

        if it > 1 and erro < tol:
            return x3, it

        melhor = (fx3 < fx1) if minimizar else (fx3 > fx1)

        if x3 > x1:
            if melhor:
                x0, fx0 = x1, fx1
                x1, fx1 = x3, fx3
            else:
                x2, fx2 = x3, fx3
        else:
            if melhor:
                x2, fx2 = x1, fx1
                x1, fx1 = x3, fx3
            else:
                x0, fx0 = x3, fx3

        x3_ant = x3

    return x3, max_it


def f(x):
    return (x - 2)**2 + 1


x_otimo, iteracoes = interpolacao_quadratica(
    f,
    x0=0,
    x1=1,
    x2=3,
    tol=1e-6,
    max_it=100,
    minimizar=True
)

print("Ponto ótimo:", x_otimo)
print("Valor da função:", f(x_otimo))
print("Número de iterações:", iteracoes)
```

### 5.2. Interpretação

O algoritmo começa com três pontos:

\[
x_0=0,\qquad x_1=1,\qquad x_2=3
\]

A função é avaliada nesses pontos e o polinômio interpolador é construído.

Como a função utilizada no exemplo já é quadrática, a interpolação com três pontos distintos recupera exatamente essa parábola em aritmética exata.

A posição do vértice calculada será \(x=2\), desconsiderando pequenos erros de arredondamento.

A implementação realiza pelo menos duas iterações antes de retornar por convergência, pois o critério de parada exige `it > 1`.

Esse exemplo é útil para verificar o funcionamento básico do código, mas uma função quadrática não representa todos os desafios que podem aparecer quando o método é aplicado a funções não quadráticas.

---

## 6. Exemplo de maximização

O mesmo algoritmo pode ser utilizado para tentar encontrar um máximo.

Considere:

\[
f(x)=-(x-2)^2+5
\]

Essa função é uma parábola voltada para baixo, portanto apresenta um máximo em:

\[
x_{\max}=2
\]

O valor máximo é:

\[
f(2)=5
\]

Para executar o algoritmo com o objetivo de maximização, basta utilizar `minimizar=False`:

```python
def f(x):
    return -(x - 2)**2 + 5


x_otimo, iteracoes = interpolacao_quadratica(
    f,
    x0=0,
    x1=1,
    x2=3,
    tol=1e-6,
    max_it=100,
    minimizar=False
)

print("Ponto ótimo:", x_otimo)
print("Valor da função:", f(x_otimo))
print("Número de iterações:", iteracoes)
```

O método tentará encontrar a posição do máximo utilizando a mesma fórmula de interpolação quadrática.

Como na minimização, é necessário verificar o resultado e não presumir que o método sempre convergirá.

---

## 7. Verificação dos resultados

Depois de obter a estimativa do extremo, é possível avaliar a função nesse ponto:

```python
valor_otimo = f(x_otimo)

print("Ponto estimado:", x_otimo)
print("Valor da função nesse ponto:", valor_otimo)
```

Para funções diferenciáveis, uma verificação adicional consiste em avaliar a derivada próxima do ponto encontrado.

Em um extremo interior suave, espera-se que:

\[
f'(x_{\text{ótimo}})\approx 0
\]

Entretanto, essa condição não é suficiente para identificar, isoladamente, se o ponto é um mínimo ou um máximo.

A segunda derivada pode ajudar na classificação local:

- Se \(f''(x_{\text{ótimo}})>0\), há um mínimo local em condições usuais.
- Se \(f''(x_{\text{ótimo}})<0\), há um máximo local em condições usuais.

Essas verificações são informações adicionais e não fazem parte do código original.

Também é importante lembrar que um ponto estacionário pode não ser um extremo e que o método não garante encontrar um extremo global.

---

## 8. Limitações e cuidados

### 8.1. Escolha dos pontos iniciais

A qualidade da aproximação depende dos três pontos iniciais.

É recomendável escolher pontos distintos e localizados em uma região em que se espera encontrar o extremo.

A implementação não verifica explicitamente se os pontos iniciais estão ordenados.

### 8.2. Denominador próximo de zero

Se o denominador da fórmula for muito pequeno, o cálculo do vértice pode se tornar numericamente instável.

O código interrompe a execução quando:

\[
|D|<10^{-15}
\]

Essa proteção evita uma divisão por um valor considerado muito pequeno, mas não tenta selecionar novos pontos automaticamente.

### 8.3. Convergência não garantida

O método não garante convergência para qualquer função ou qualquer conjunto de pontos iniciais.

A sequência de estimativas pode convergir, oscilar ou produzir resultados inadequados para o problema.

### 8.4. Critério de parada limitado

O código utiliza:

\[
|x_3-x_{3,\text{ant}}|<\text{tol}
\]

Esse critério verifica a mudança na posição estimada, mas não verifica diretamente a derivada da função, a qualidade do polinômio interpolador ou a diferença entre o valor encontrado e o extremo verdadeiro.

### 8.5. Ausência de indicador explícito de convergência

Se o limite de iterações for atingido, a função retorna a última estimativa e `max_it`.

O retorno não distingue explicitamente uma convergência bem-sucedida de uma interrupção por limite de iterações.

### 8.6. Mínimos e máximos locais

O algoritmo utiliza comparações locais entre os valores da função para atualizar os pontos.

Isso não garante que o extremo encontrado seja global. Uma função pode apresentar vários mínimos e máximos em diferentes regiões.

### 8.7. Validação da função

A implementação pressupõe que a função fornecida possa ser avaliada nos pontos gerados.

Funções com descontinuidades, domínios restritos ou valores não finitos podem causar problemas que não são tratados explicitamente pelo código.

### 8.8. Avaliação de funções não quadráticas

Para uma função que não seja quadrática, o polinômio interpolador representa apenas uma aproximação local baseada nos três pontos.

Por isso, o vértice do polinômio pode não coincidir com o extremo real da função.

---

## 9. Complexidade computacional

Em cada iteração, o código realiza um número constante de operações aritméticas e algumas avaliações da função.

Desconsiderando o custo interno da função \(f\), o custo aritmético de cada iteração é:

\[
O(1)
\]

Se o algoritmo realizar \(k\) iterações, o custo aritmético total é:

\[
O(k)
\]

Entretanto, o custo real também depende do tempo necessário para calcular \(f(x)\). Se a avaliação da função for computacionalmente cara, ela poderá dominar o tempo total de execução.

O número de iterações depende da tolerância, dos pontos iniciais, da função e do comportamento da sequência de aproximações.

---

## 10. Possíveis melhorias na implementação

A implementação pode ser aprimorada com as seguintes funcionalidades:

1. **Validação dos pontos iniciais:** verificar se os pontos são distintos e se estão em uma ordem apropriada.
2. **Indicador de convergência:** informar se o critério de parada foi satisfeito ou se o limite de iterações foi atingido.
3. **Verificação da qualidade do extremo:** utilizar critérios adicionais, como a derivada ou a variação dos valores da função.
4. **Tratamento de denominadores pequenos:** permitir uma estratégia de recuperação em vez de interromper imediatamente a execução.
5. **Controle do domínio:** garantir que os pontos calculados estejam dentro do domínio permitido para a função.
6. **Registro das iterações:** armazenar os pontos, valores da função e erros para acompanhar a evolução do algoritmo.
7. **Critério de parada adicional:** combinar a diferença entre estimativas com a variação do valor da função ou outros indicadores adequados ao problema.

Essas melhorias são sugestões para versões futuras e não fazem parte da implementação original.

---

## 11. Resumo do algoritmo

O método implementado pode ser resumido nas seguintes etapas:

1. Receber a função e três pontos iniciais.
2. Avaliar a função nos três pontos.
3. Calcular o numerador e o denominador da fórmula de interpolação quadrática.
4. Verificar se o denominador está suficientemente afastado de zero.
5. Calcular a posição estimada do vértice do polinômio interpolador.
6. Avaliar a função no novo ponto.
7. Calcular a diferença entre a estimativa atual e a anterior.
8. Verificar se o critério de parada foi satisfeito.
9. Comparar o valor da função no novo ponto com o valor no ponto central.
10. Atualizar os três pontos de referência.
11. Repetir o procedimento até satisfazer o critério de parada ou atingir o limite de iterações.
12. Retornar a estimativa final e o número de iterações.

O método da interpolação quadrática é uma ferramenta de otimização unidimensional que utiliza uma aproximação parabólica para estimar extremos de funções. Sua aplicação exige uma escolha adequada dos pontos iniciais e a verificação dos resultados obtidos.

Para obter resultados confiáveis, é importante considerar as limitações do critério de parada, a possibilidade de instabilidade numérica e a diferença entre encontrar um extremo local e encontrar um extremo global.