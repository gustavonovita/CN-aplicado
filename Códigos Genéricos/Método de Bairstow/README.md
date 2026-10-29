# Método de Bairstow para Encontrar Raízes de Polinômios

## 1. Objetivo

Este projeto implementa o **Método de Bairstow**, um método numérico utilizado para encontrar raízes reais e complexas de polinômios.

O método consiste em determinar fatores quadráticos do polinômio por meio de aproximações sucessivas dos parâmetros desses fatores. Após encontrar um fator quadrático, o algoritmo calcula suas duas raízes, reduz o grau do polinômio e repete o procedimento até que reste um polinômio de grau 2 ou 1.

A função implementada no código é:

```python
def bairstow(a, r0=0.1, s0=0.1, tol=1e-6, max_it=100):
```

A função recebe os coeficientes do polinômio, realiza as iterações numéricas e retorna uma lista com as raízes calculadas.

## 2. Representação do polinômio

O código recebe os coeficientes do polinômio em **ordem crescente das potências de x**.

Um polinômio de grau \(n\) pode ser escrito como:

\[
P(x)=a_0+a_1x+a_2x^2+\cdots+a_nx^n
\]

A lista de coeficientes deve seguir esta ordem:

```python
a = [a0, a1, a2, ..., an]
```

### Exemplo

Considere o polinômio:

\[
P(x)=x^3-6x^2+11x-6
\]

Seus coeficientes são:

- \(a_0=-6\)
- \(a_1=11\)
- \(a_2=-6\)
- \(a_3=1\)

Portanto, a representação correta é:

```python
coeficientes = [-6, 11, -6, 1]
```

**Atenção:** não forneça os coeficientes na ordem decrescente das potências. A função espera primeiro o termo constante e, por último, o coeficiente da maior potência.

---

## 3. Fundamentação matemática

### 3.1. Fator quadrático

O método de Bairstow procura um fator quadrático da forma:

\[
x^2-rx-s
\]

Os parâmetros \(r\) e \(s\) são determinados numericamente.

Quando esse fator é encontrado, suas raízes podem ser calculadas resolvendo:

\[
x^2-rx-s=0
\]

Aplicando a fórmula de Bhaskara:

\[
x=\frac{r\pm\sqrt{r^2+4s}}{2}
\]

O discriminante utilizado pelo código é:

\[
\Delta=r^2+4s
\]

Sua interpretação é:

- Se \(\Delta>0\), existem duas raízes reais distintas.
- Se \(\Delta=0\), existem duas raízes reais iguais.
- Se \(\Delta<0\), existem duas raízes complexas conjugadas.

### 3.2. Divisão sintética

Para determinar os parâmetros \(r\) e \(s\), o algoritmo utiliza uma recorrência semelhante à divisão sintética de polinômios.

O vetor `b` armazena os coeficientes calculados durante essa operação.

Os valores `b[0]` e `b[1]` representam os resíduos associados à tentativa de divisão pelo fator quadrático. O objetivo é ajustar os parâmetros até que esses resíduos sejam suficientemente pequenos.

Idealmente:

\[
b_0\approx 0
\]

\[
b_1\approx 0
\]

Quando isso ocorre, o fator quadrático é aproximadamente um divisor do polinômio.

### 3.3. Correção dos parâmetros

O método também calcula um vetor auxiliar `c`, utilizado para determinar as correções dos parâmetros.

As correções são armazenadas em:

- `dr`: correção de \(r\).
- `ds`: correção de \(s\).

Os parâmetros são atualizados por:

\[
r_{\text{novo}}=r+\Delta r
\]

\[
s_{\text{novo}}=s+\Delta s
\]

O processo é repetido até que o critério de convergência seja satisfeito ou o número máximo de iterações seja atingido.

---

## 4. Parâmetros da função

A função possui cinco parâmetros:

```python
bairstow(a, r0=0.1, s0=0.1, tol=1e-6, max_it=100)
```

| Parâmetro | Descrição |
|---|---|
| `a` | Lista de coeficientes do polinômio, em ordem crescente das potências de \(x\). |
| `r0` | Estimativa inicial do parâmetro \(r\). O valor padrão é `0.1`. |
| `s0` | Estimativa inicial do parâmetro \(s\). O valor padrão é `0.1`. |
| `tol` | Tolerância utilizada no critério de convergência. O valor padrão é `1e-6`. |
| `max_it` | Número máximo de iterações realizadas para encontrar cada fator quadrático. O valor padrão é `100`. |

Os valores iniciais de `r0` e `s0` podem influenciar a convergência. Dependendo do polinômio, pode ser necessário experimentar outras estimativas iniciais.

A tolerância define quão pequenas devem ser as correções dos parâmetros ou os resíduos para que o algoritmo interrompa as iterações.

---

## 5. Explicação detalhada do código

### 5.1. Definição da função e inicialização

```python
def bairstow(a, r0=0.1, s0=0.1, tol=1e-6, max_it=100):
    # a = [a0, a1, ..., an] onde a[i] multiplica x^i
    coefs = [float(val) for val in a]
    raizes = []
    r, s = float(r0), float(s0)
```

**Definição da função:**

A primeira linha define a função `bairstow` e seus parâmetros.

**Conversão dos coeficientes:**

```python
coefs = [float(val) for val in a]
```

Converte todos os coeficientes para números de ponto flutuante. Isso permite realizar as operações numéricas com valores decimais.

**Lista de raízes:**

```python
raizes = []
```

Cria uma lista vazia para armazenar as raízes encontradas ao longo da execução.

**Estimativas iniciais:**

```python
r, s = float(r0), float(s0)
```

Inicializa os parâmetros do fator quadrático com os valores fornecidos pelo usuário.

### 5.2. Laço para encontrar fatores quadráticos

```python
while len(coefs) - 1 >= 3:
    n = len(coefs) - 1
    b = [0.0] * (n + 1)
    c = [0.0] * (n + 1)
```

O laço `while` continua enquanto o grau do polinômio restante for igual ou superior a 3.

O grau é calculado por:

\[
n=\operatorname{len}(\text{coefs})-1
\]

O vetor `b` é utilizado para calcular os coeficientes da divisão sintética, enquanto o vetor `c` é utilizado no cálculo das correções de \(r\) e \(s\).

Os vetores são inicializados com zeros e possuem tamanho `n + 1`.

A cada fator quadrático encontrado, o grau do polinômio restante diminui em duas unidades.

### 5.3. Cálculo do vetor `b`

```python
for _ in range(max_it):
    b[n] = coefs[n]
    b[n - 1] = coefs[n - 1] + r * b[n]

    for i in range(n - 2, -1, -1):
        b[i] = coefs[i] + r * b[i + 1] + s * b[i + 2]
```

Esse trecho realiza os cálculos associados à divisão sintética.

Primeiramente:

```python
b[n] = coefs[n]
```

O último coeficiente de `b` recebe o coeficiente de maior grau do polinômio.

Em seguida:

```python
b[n - 1] = coefs[n - 1] + r * b[n]
```

O penúltimo coeficiente é calculado utilizando o parâmetro \(r\).

Por fim, o laço:

```python
for i in range(n - 2, -1, -1):
    b[i] = coefs[i] + r * b[i + 1] + s * b[i + 2]
```

Calcula os demais coeficientes de `b`, percorrendo os índices em ordem decrescente.

A recorrência utilizada é:

\[
b_i=a_i+rb_{i+1}+sb_{i+2}
\]

Os elementos `b[0]` e `b[1]` são utilizados como resíduos para avaliar a qualidade da aproximação do fator quadrático.

### 5.4. Cálculo do vetor auxiliar `c`

```python
c[n] = b[n]
c[n - 1] = b[n - 1] + r * c[n]

for i in range(n - 2, 0, -1):
    c[i] = b[i] + r * c[i + 1] + s * c[i + 2]
```

O vetor `c` é calculado a partir dos elementos do vetor `b`.

Primeiro, são definidos os dois últimos elementos relevantes:

```python
c[n] = b[n]
c[n - 1] = b[n - 1] + r * c[n]
```

Depois, os demais valores são calculados pela recorrência:

\[
c_i=b_i+rc_{i+1}+sc_{i+2}
\]

Esse vetor é utilizado para calcular as correções dos parâmetros \(r\) e \(s\).

### 5.5. Determinante e prevenção de divisão por zero

```python
det = c[2] * c[2] - c[1] * c[3]

if abs(det) < 1e-15:
    r += 0.1
    s += 0.1
    continue
```

O código calcula um determinante:

\[
D=c_2^2-c_1c_3
\]

Esse valor é utilizado como denominador nas expressões de `dr` e `ds`.

Se o módulo do determinante for inferior a \(10^{-15}\), o algoritmo considera que ele está muito próximo de zero.

Nesse caso:

```python
r += 0.1
s += 0.1
```

Os parâmetros são incrementados em `0.1`.

O comando:

```python
continue
```

interrompe a iteração atual e retorna ao início do laço, tentando realizar os cálculos novamente com os novos valores.

Essa estratégia evita a divisão por um determinante considerado muito pequeno, mas não garante a convergência. Além disso, o incremento fixo pode não ser adequado para todos os polinômios.

### 5.6. Cálculo das correções `dr` e `ds`

```python
dr = (-b[1] * c[2] + b[0] * c[3]) / det
ds = (-b[0] * c[2] + b[1] * c[1]) / det
```

Essas expressões calculam as correções utilizadas para atualizar os parâmetros.

Matematicamente:

\[
\Delta r=\frac{-b_1c_2+b_0c_3}{D}
\]

\[
\Delta s=\frac{-b_0c_2+b_1c_1}{D}
\]

As correções dependem dos resíduos `b[0]` e `b[1]`, dos coeficientes auxiliares `c` e do determinante `det`.

### 5.7. Atualização dos parâmetros

```python
r += dr
s += ds
```

Equivale a:

\[
r_{\text{novo}}=r+\Delta r
\]

\[
s_{\text{novo}}=s+\Delta s
\]

Os parâmetros atualizados representam uma nova estimativa do fator quadrático.

O processo é repetido para tentar diminuir os resíduos e aproximar o fator quadrático de um divisor do polinômio.

### 5.8. Critério de convergência

```python
if max(abs(dr), abs(ds)) < tol or max(abs(b[0]), abs(b[1])) < tol:
    break
```

O código verifica se uma das duas condições foi satisfeita.

**Primeira condição:**

```python
max(abs(dr), abs(ds)) < tol
```

Verifica se as correções de \(r\) e \(s\) são menores que a tolerância.

Isso indica que os parâmetros estão sofrendo alterações pequenas.

**Segunda condição:**

```python
max(abs(b[0]), abs(b[1])) < tol
```

Verifica se os dois resíduos são menores que a tolerância.

**Operador `or`:**

O operador lógico `or` significa “ou”. Portanto, basta que uma das condições seja verdadeira para que o laço seja interrompido.

O comando `break` encerra o laço de iterações.

Se nenhuma condição for satisfeita, o laço termina quando atingir `max_it`. O código não registra explicitamente se houve convergência ou se o limite de iterações foi atingido.

### 5.9. Cálculo do discriminante

```python
disc = r * r + 4.0 * s
```

Após o término das iterações, o código calcula o discriminante do fator quadrático:

\[
\Delta=r^2+4s
\]

O valor de `disc` determina se as raízes calculadas serão reais ou complexas.

### 5.10. Cálculo das raízes reais

```python
if disc >= 0.0:
    raizes.append((r + disc**0.5) / 2.0)
    raizes.append((r - disc**0.5) / 2.0)
```

Quando o discriminante é não negativo, o código calcula duas raízes reais:

\[
x_1=\frac{r+\sqrt{\Delta}}{2}
\]

\[
x_2=\frac{r-\sqrt{\Delta}}{2}
\]

O método `append()` adiciona cada raiz à lista `raizes`.

O operador `**0.5` calcula a raiz quadrada de um número não negativo.

### 5.11. Cálculo das raízes complexas

```python
else:
    imag = (-disc) ** 0.5 / 2.0
    raizes.append(complex(r / 2.0, imag))
    raizes.append(complex(r / 2.0, -imag))
```

Quando o discriminante é negativo, as raízes são complexas.

A parte real é:

\[
\operatorname{Re}(x)=\frac{r}{2}
\]

A magnitude da parte imaginária é:

\[
\operatorname{Im}(x)=\frac{\sqrt{-\Delta}}{2}
\]

O código cria as duas raízes conjugadas:

\[
x_1=\frac{r}{2}+i\frac{\sqrt{-\Delta}}{2}
\]

\[
x_2=\frac{r}{2}-i\frac{\sqrt{-\Delta}}{2}
\]

A função `complex()` representa esses valores como números complexos em Python.

### 5.12. Redução do grau do polinômio

```python
coefs = b[2:]
```

Após encontrar um fator quadrático, o código utiliza os elementos de `b` a partir do índice 2 para formar os coeficientes do polinômio restante.

Essa operação reduz o grau do polinômio em duas unidades e permite que o algoritmo procure outro fator quadrático.

O processo se repete enquanto o grau do polinômio restante for pelo menos 3.

**Observação sobre a implementação:** o vetor `b` foi calculado antes da atualização final de `r` e `s` naquela iteração. Se o laço terminar logo após atualizar os parâmetros, as raízes podem ser calculadas com os parâmetros atualizados, enquanto `b[2:]` pode corresponder à divisão calculada com os parâmetros anteriores. Isso pode introduzir inconsistência numérica. Uma implementação mais robusta recalcularia o vetor `b` com os parâmetros finais antes de reduzir o grau.

### 5.13. Determinação do grau restante

```python
grau_restante = len(coefs) - 1
```

Depois que o laço principal termina, o código determina o grau do polinômio restante.

Existem três possibilidades principais:

- Grau 2: calcula as duas raízes pela fórmula quadrática.
- Grau 1: calcula diretamente uma raiz.
- Grau 0: não adiciona novas raízes.

### 5.14. Resolução direta do polinômio de segundo grau

```python
if grau_restante == 2:
    r_fin = -coefs[1] / coefs[2]
    s_fin = -coefs[0] / coefs[2]
    disc = r_fin * r_fin + 4.0 * s_fin
```

Quando o polinômio restante é de segundo grau, não é necessário continuar com as iterações de Bairstow.

O polinômio pode ser escrito como:

\[
P(x)=a_0+a_1x+a_2x^2
\]

Dividindo os coeficientes por \(a_2\), obtém-se a forma normalizada:

\[
x^2+\frac{a_1}{a_2}x+\frac{a_0}{a_2}=0
\]

Comparando com \(x^2-rx-s=0\), o código define:

\[
r_{\text{fin}}=-\frac{a_1}{a_2}
\]

\[
s_{\text{fin}}=-\frac{a_0}{a_2}
\]

Em seguida, calcula:

\[
\Delta=r_{\text{fin}}^2+4s_{\text{fin}}
\]

Se o discriminante for não negativo, adiciona as duas raízes reais à lista. Caso seja negativo, adiciona as duas raízes complexas conjugadas.

### 5.15. Resolução direta do polinômio de primeiro grau

```python
elif grau_restante == 1:
    raizes.append(-coefs[0] / coefs[1])
```

Se restar um polinômio de primeiro grau:

\[
P(x)=a_0+a_1x
\]

A raiz é determinada por:

\[
a_0+a_1x=0
\]

Isolando \(x\):

\[
a_1x=-a_0
\]

\[
x=-\frac{a_0}{a_1}
\]

O código calcula essa raiz e a adiciona à lista.

### 5.16. Retorno das raízes

```python
return raizes
```

Ao final, a função retorna a lista `raizes`, que pode conter valores reais e complexos.

A ordem das raízes não é garantida e pode variar de acordo com o processo numérico.

---

## 6. Exemplo de utilização

Considere o polinômio:

\[
P(x)=x^3-6x^2+11x-6
\]

Sua lista de coeficientes, na ordem exigida pelo código, é:

```python
coeficientes = [-6, 11, -6, 1]
```

O exemplo de chamada é:

```python
coeficientes = [-6, 11, -6, 1]

raizes = bairstow(coeficientes)

print("Raízes encontradas:", raizes)
```

Para executar esse exemplo, é necessário que a função `bairstow` esteja definida no arquivo Python ou importada de outro módulo.

As raízes exatas desse polinômio são:

\[
x_1=1,\qquad x_2=2,\qquad x_3=3
\]

Esses valores permitem comparar as raízes calculadas pelo método com as raízes conhecidas.

Como o método é numérico, podem ocorrer pequenas diferenças devido ao arredondamento e à tolerância escolhida.

**Importante:** esse exemplo ilustra como chamar a função, mas não garante que a implementação convergirá para as raízes corretas em todos os casos. O resultado deve ser verificado numericamente.

---

## 7. Verificação das raízes encontradas

Uma forma de avaliar a qualidade das raízes calculadas é substituí-las no polinômio original.

Para uma raiz aproximada \(x_i\), calcula-se o resíduo:

\[
R_i=P(x_i)
\]

Se \(x_i\) for uma boa aproximação de uma raiz, espera-se que:

\[
|P(x_i)|\approx 0
\]

Em Python, essa verificação pode ser feita com:

```python
def avaliar_polinomio(coefs, x):
    resultado = 0

    for i, coef in enumerate(coefs):
        resultado += coef * x**i

    return resultado


coeficientes = [-6, 11, -6, 1]
raizes = bairstow(coeficientes)

for raiz in raizes:
    residuo = avaliar_polinomio(coeficientes, raiz)
    print(f"Raiz: {raiz}, resíduo: {residuo}")
```

A função auxiliar `avaliar_polinomio` utiliza a mesma convenção de coeficientes do método de Bairstow, isto é, do termo constante ao termo de maior grau.

Essa verificação também funciona para raízes complexas.

Um resíduo pequeno é um indicador útil, mas deve ser interpretado considerando a escala dos coeficientes e os erros de arredondamento.

---

## 8. Limitações e cuidados

### 8.1. Ordem dos coeficientes

A função espera os coeficientes em ordem crescente das potências de \(x\). Informá-los em outra ordem altera o polinômio representado.

### 8.2. Dependência das estimativas iniciais

Os parâmetros `r0` e `s0` podem afetar a convergência. Alguns polinômios podem exigir estimativas iniciais diferentes das fornecidas por padrão.

### 8.3. Número máximo de iterações

O parâmetro `max_it` limita o número de iterações para cada fator quadrático.

Atingir esse limite não significa necessariamente que o método convergiu. A implementação original não retorna um indicador explícito de convergência.

### 8.4. Determinante próximo de zero

Quando o determinante calculado está muito próximo de zero, o código aumenta `r` e `s` em `0.1` e tenta novamente.

Essa estratégia pode evitar uma divisão por um valor muito pequeno, mas não garante que o cálculo seguinte será estável ou que haverá convergência.

### 8.5. Tolerância absoluta

A tolerância `tol` é utilizada como limite absoluto para as correções ou os resíduos.

Como não existe uma normalização automática pela escala dos coeficientes, o comportamento pode variar conforme a magnitude dos valores presentes no polinômio.

### 8.6. Validação da entrada

A função não realiza validações explícitas para todos os casos especiais.

Por exemplo:

- Uma lista vazia não representa um polinômio válido para esse procedimento.
- Um polinômio constante não possui raízes determinadas por esse algoritmo.
- Um coeficiente líder igual a zero pode indicar que o grau real é menor do que o esperado.
- Divisões por coeficientes nulos podem provocar erros.
- Coeficientes não finitos, como `NaN` ou infinito, podem prejudicar os cálculos.

Por isso, é recomendável validar os coeficientes antes de chamar a função.

### 8.7. Consistência do polinômio restante

O vetor `b` pode ter sido calculado com estimativas anteriores de `r` e `s`, enquanto as raízes são calculadas com os parâmetros após a última atualização.

Essa diferença pode fazer com que o polinômio reduzido não corresponda exatamente ao fator utilizado para calcular as raízes. Recalcular a divisão sintética com os parâmetros finais é uma melhoria recomendável.

### 8.8. Verificação dos resultados

As raízes retornadas devem ser verificadas substituindo-as no polinômio original. Essa etapa ajuda a identificar resultados imprecisos ou problemas de convergência.

---

## 9. Complexidade computacional

Em cada iteração, o algoritmo percorre os coeficientes do polinômio para calcular os vetores `b` e `c`. O custo de uma iteração é, portanto, proporcional ao grau atual do polinômio.

Como o grau diminui em duas unidades após cada fator quadrático, o custo total depende de:

- Grau inicial do polinômio.
- Número de iterações realizadas para cada fator.
- Valor de `max_it`.
- Comportamento de convergência para as estimativas iniciais escolhidas.

O método pode ser adequado para exercícios acadêmicos e polinômios de tamanho moderado. Para aplicações que exigem maior confiabilidade, é recomendável adicionar validações, critérios de convergência mais robustos e testes sistemáticos dos resíduos.

---

## 10. Possíveis melhorias na implementação

Algumas melhorias que podem aumentar a confiabilidade do código são:

1. **Validar a entrada:** verificar se a lista contém coeficientes válidos e se o coeficiente de maior grau é diferente de zero.
2. **Informar a convergência:** retornar ou registrar se o método convergiu ou atingiu o limite de iterações.
3. **Recalcular os coeficientes restantes:** realizar novamente a divisão sintética usando os valores finais de `r` e `s`.
4. **Melhorar o tratamento do determinante:** substituir incrementos fixos por uma estratégia de reinicialização ou recuperação mais adequada.
5. **Utilizar critérios de parada relativos:** considerar a escala dos coeficientes, além de tolerâncias absolutas.
6. **Verificar os resíduos:** avaliar cada raiz no polinômio original.
7. **Criar testes automatizados:** comparar os resultados com polinômios cujas raízes são conhecidas.

Essas melhorias não fazem parte da implementação original apresentada, mas são sugestões para futuras versões.

---

## 11. Resumo do algoritmo

O método implementado pode ser resumido nas seguintes etapas:

1. Recebe os coeficientes do polinômio em ordem crescente das potências.
2. Inicializa os parâmetros \(r\) e \(s\).
3. Calcula o vetor `b` por recorrência.
4. Calcula o vetor auxiliar `c`.
5. Determina as correções `dr` e `ds`.
6. Atualiza os parâmetros até satisfazer um critério de parada ou atingir o limite de iterações.
7. Calcula as duas raízes do fator quadrático, reais ou complexas.
8. Reduz o grau do polinômio restante.
9. Repete o processo enquanto o grau for pelo menos 3.
10. Resolve diretamente o polinômio final de grau 2 ou 1.
11. Retorna a lista de raízes encontradas.

O Método de Bairstow é uma técnica numérica útil para calcular raízes de polinômios por meio da extração sucessiva de fatores quadráticos. A precisão dos resultados depende da convergência das iterações, das estimativas iniciais e da estabilidade dos cálculos.

Para obter resultados confiáveis, é importante verificar as raízes encontradas substituindo-as no polinômio original.
