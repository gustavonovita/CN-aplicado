# Simulação da Distribuição de Temperatura em uma Placa pelo Método de Gauss-Seidel

## 1. Visão geral do projeto

Este projeto implementa um método numérico para simular a distribuição de temperatura em uma placa bidimensional utilizando o **método iterativo de Gauss-Seidel**.

O código utiliza a linguagem Python e a biblioteca NumPy para armazenar e atualizar os valores de temperatura de uma malha numérica.

O objetivo é determinar a temperatura nos nós internos da placa a partir das condições de contorno estabelecidas e da equação de equilíbrio térmico utilizada pelo programa.

A simulação é executada para três tamanhos de malha:

- \(n=5\);
- \(n=10\);
- \(n=20\).

Para cada tamanho, o programa apresenta o número de iterações necessárias, a variação máxima de temperatura na última iteração e o resíduo máximo calculado. Para \(n=5\), também apresenta a temperatura do nó central e uma matriz de temperaturas reorganizada para visualização.

### Principais características

- Utilização do método iterativo de Gauss-Seidel;
- Representação da placa por uma matriz bidimensional;
- Condição de contorno de \(100\,^\circ\mathrm{C}\) em uma das bordas;
- Temperatura inicial de \(20\,^\circ\mathrm{C}\) para os nós;
- Critério de parada baseado na variação máxima de temperatura;
- Limite de 10.000 iterações;
- Cálculo do resíduo máximo após a simulação;
- Comparação de diferentes tamanhos de malha.

## 2. Modelo físico e interpretação

O código representa uma placa dividida em uma malha de pontos. Cada ponto interno possui uma temperatura que é atualizada a partir das temperaturas dos quatro vizinhos imediatos.

A equação utilizada é:

\[
T_{i,j}=
\frac{
T_{i-1,j}+T_{i+1,j}+T_{i,j-1}+T_{i,j+1}
}{4}
\]

Em que:

- \(T_{i,j}\) é a temperatura do nó atual;
- \(T_{i-1,j}\) é a temperatura do vizinho em uma direção;
- \(T_{i+1,j}\) é a temperatura do vizinho na direção oposta;
- \(T_{i,j-1}\) e \(T_{i,j+1}\) são as temperaturas dos outros dois vizinhos.

A equação determina que a temperatura de cada nó interno seja igual à média aritmética das temperaturas dos seus quatro vizinhos.

Essa formulação é frequentemente utilizada na solução numérica da equação de Laplace para problemas bidimensionais de condução de calor em regime estacionário, sob hipóteses adequadas.

### 2.1. Equação governante

Para uma placa homogênea, com propriedades constantes, sem geração interna de calor e em regime estacionário, a condução bidimensional pode ser descrita pela equação de Laplace:

\[
\frac{\partial^2 T}{\partial x^2}
+
\frac{\partial^2 T}{\partial y^2}=0
\]

A discretização por diferenças finitas, considerando espaçamentos iguais nas duas direções, resulta em:

\[
\frac{T_{i-1,j}-2T_{i,j}+T_{i+1,j}}{\Delta x^2}
+
\frac{T_{i,j-1}-2T_{i,j}+T_{i,j+1}}{\Delta y^2}
=0
\]

Se \(\Delta x=\Delta y\), a expressão pode ser reorganizada para:

\[
T_{i,j}=
\frac{
T_{i-1,j}+T_{i+1,j}+T_{i,j-1}+T_{i,j+1}
}{4}
\]

Essa é a relação utilizada no código.

**Limitação do modelo:** o programa não recebe dimensões físicas da placa nem os espaçamentos da malha. Portanto, considera implicitamente a formulação discreta com espaçamentos iguais. Também não modela explicitamente geração interna de calor ou variações das propriedades térmicas.

## 3. Condições iniciais e condições de contorno

As condições de contorno determinam as temperaturas nas bordas da placa.

### 3.1. Temperatura inicial

O programa inicializa toda a matriz com temperatura de \(20\,^\circ\mathrm{C}\):

```python
T = np.full((n + 2, n + 2), 20.0)
```

Isso significa que, inicialmente, todos os elementos da matriz, inclusive os que representam os nós internos e as bordas, recebem o valor 20.

Em seguida, o código modifica uma das bordas para estabelecer a condição de contorno de \(100\,^\circ\mathrm{C}\).

### 3.2. Borda a 100 °C

```python
T[1:n+1, 0] = 100.0
```

Essa instrução atribui \(100\,^\circ\mathrm{C}\) aos elementos cujo primeiro índice varia de 1 até \(n\) e cujo segundo índice é zero.

Na orientação utilizada para imprimir a matriz, essa borda corresponde à **base da placa**.

O uso de `1:n+1` é necessário porque o limite final de um intervalo de fatiamento em Python não é incluído. Assim, os índices selecionados são:

\[
1,2,\ldots,n
\]

O programa não atualiza explicitamente as demais bordas. Elas permanecem em \(20\,^\circ\mathrm{C}\), que é o valor inicial atribuído a todos os elementos dessas bordas.

### 3.3. Resumo das condições do código

| Região | Temperatura ou comportamento |
|---|---|
| Borda inferior na visualização | \(100\,^\circ\mathrm{C}\) |
| Demais bordas | \(20\,^\circ\mathrm{C}\) |
| Nós internos, no início | \(20\,^\circ\mathrm{C}\) |
| Nós internos, durante a simulação | Atualizados pelo método de Gauss-Seidel |

As condições de contorno são mantidas durante as iterações porque os laços de atualização percorrem somente os índices de 1 até \(n\). Os índices das bordas não são atualizados nesses laços.

## 4. Estrutura da malha

A matriz utilizada pelo programa possui dimensão:

\[
(n+2)\times(n+2)
\]

O acréscimo de dois elementos em cada dimensão permite representar os nós internos e as bordas.

Para uma malha com parâmetro \(n\):

- Existem \(n\times n\) nós internos atualizados pelo método;
- Existem duas posições adicionais em cada dimensão para representar os limites da matriz.

### 4.1. Exemplo para n = 5

Quando \(n=5\):

\[
\text{dimensão da matriz}=(5+2)\times(5+2)
\]

\[
\boxed{7\times7}
\]

A matriz possui 49 elementos. Os laços atualizam somente os 25 nós internos:

\[
5\times5=25
\]

Os demais elementos representam as bordas.

### 4.2. Exemplo para n = 10

Quando \(n=10\):

\[
\text{dimensão da matriz}=12\times12
\]

O programa atualiza:

\[
10\times10=100
\]

nós internos.

### 4.3. Exemplo para n = 20

Quando \(n=20\):

\[
\text{dimensão da matriz}=22\times22
\]

O número de nós internos atualizados é:

\[
20\times20=400
\]

### 4.4. Visualização conceitual

Para \(n=5\), a organização dos índices internos pode ser representada de maneira simplificada por:

```text
              Borda superior
       +-------------------------+
       |                         |
       |   T(1,5) ... T(5,5)     |
       |   T(1,4) ... T(5,4)     |
       |   T(1,3) ... T(5,3)     |
       |   T(1,2) ... T(5,2)     |
       |   T(1,1) ... T(5,1)     |
       |                         |
       +-------------------------+
              Borda inferior
                 100 °C
```

Esse desenho é apenas conceitual. A posição visual exata dos índices depende da transposição e da inversão aplicadas quando a matriz é impressa.

## 5. Método de Gauss-Seidel

O método de Gauss-Seidel é uma técnica iterativa para resolver sistemas de equações lineares.

Em vez de determinar todas as incógnitas simultaneamente, ele atualiza uma incógnita de cada vez, utilizando imediatamente os valores mais recentes disponíveis.

No problema da placa, cada atualização utiliza a média das temperaturas dos quatro vizinhos.

A fórmula implementada é:

\[
T_{i,j}^{(k+1)}=
\frac{
T_{i-1,j}+T_{i+1,j}+T_{i,j-1}+T_{i,j+1}
}{4}
\]

A notação \(k+1\) representa a nova atualização. Como o método é Gauss-Seidel, os valores dos vizinhos já atualizados na iteração atual são utilizados imediatamente.

### 5.1. Diferença em relação ao método de Jacobi

No método de Jacobi, todas as novas temperaturas de uma iteração são calculadas utilizando exclusivamente os valores da iteração anterior.

No método de Gauss-Seidel, cada valor atualizado pode ser utilizado nas próximas operações da mesma iteração.

Essa característica pode reduzir o número de iterações necessárias em comparação com Jacobi para determinados problemas, embora o desempenho dependa do sistema e da malha.

## 6. Código completo

```python
import numpy as np

def resolver_placa_gauss_seidel(n):

    # Em Python, os índices vão de 0 a n+1.
    T = np.full((n + 2, n + 2), 20.0)

    # T(i,1) = 100°C para i = 2 até n+1. Em Python: i = 1 até n.
    T[1:n+1, 0] = 100.0

    tol = 1e-6
    max_iter = 10000
    iteracao = 0
    erro = float('inf')

    while erro > tol and iteracao < max_iter:
        erro = 0.0

        for j in range(1, n + 1):
            for i in range(1, n + 1):
                antigo = T[i, j]

                T[i, j] = (
                    T[i-1, j] + T[i+1, j]
                    + T[i, j-1] + T[i, j+1]
                ) / 4.0

                diferenca = abs(T[i, j] - antigo)
                erro = max(erro, diferenca)

        iteracao += 1

    max_residuo = 0.0
    for j in range(1, n + 1):
        for i in range(1, n + 1):
            r_ij = 4 * T[i, j] - (
                T[i-1, j] + T[i+1, j]
                + T[i, j-1] + T[i, j+1]
            )
            max_residuo = max(max_residuo, abs(r_ij))

    return T, iteracao, erro, max_residuo


# Simulação principal para diferentes tamanhos de malha
tamanhos_malha = [5, 10, 20]

for n in tamanhos_malha:
    print(f"--- Simulação para n = {n} ---")
    T, iteracoes, erro_final, max_residuo = resolver_placa_gauss_seidel(n)

    print(f"Iterações necessárias: {iteracoes}")
    print(f"Erro final (variação máxima): {erro_final:.2e} °C")
    print(f"Resíduo máximo calculado: {max_residuo:.2e}")

    if n == 5:
        temp_centro = T[3, 3]
        print(f"Temperatura do nó central T(4,4) para n=5: {temp_centro:.4f} °C")
        print("\nMatriz de Temperaturas (Transposta para visualização com y=0 na base):")
        # Imprime a matriz com a borda inferior a 100°C fisicamente em baixo
        print(np.flipud(T.T))
    print("\n")
```

## 7. Explicação detalhada do código

### 7.1. Importação do NumPy

```python
import numpy as np
```

Importa a biblioteca NumPy e atribui a ela o apelido `np`.

A biblioteca é utilizada para criar a matriz de temperaturas, inicializar seus valores e realizar operações com arrays.

### 7.2. Definição da função

```python
def resolver_placa_gauss_seidel(n):
```

Define uma função chamada `resolver_placa_gauss_seidel`.

O parâmetro `n` determina a quantidade de nós internos em cada direção da malha.

A função realiza a simulação e retorna quatro resultados:

1. Matriz final de temperaturas;
2. Número de iterações;
3. Erro final baseado na variação máxima;
4. Resíduo máximo das equações discretizadas.

### 7.3. Criação da matriz de temperaturas

```python
T = np.full((n + 2, n + 2), 20.0)
```

Cria uma matriz quadrada com dimensão \((n+2)\times(n+2)\).

Todos os elementos são inicialmente iguais a \(20\,^\circ\mathrm{C}\).

A matriz é chamada `T` porque armazena os valores de temperatura.

### 7.4. Aplicação da condição de contorno

```python
T[1:n+1, 0] = 100.0
```

Atribui \(100\,^\circ\mathrm{C}\) aos elementos selecionados.

O índice zero corresponde à posição inicial da segunda dimensão. Na representação final impressa pelo programa, essa borda aparece na parte inferior.

A atribuição é feita antes do processo iterativo, e esses elementos não são alterados durante a atualização dos nós internos.

### 7.5. Definição da tolerância

```python
tol = 1e-6
```

Define a tolerância como:

\[
\boxed{\text{tol}=10^{-6}=0{,}000001}
\]

A tolerância determina o limite para a variação máxima de temperatura entre atualizações.

Quanto menor a tolerância, mais rigoroso é o critério de parada. Em geral, isso pode exigir mais iterações.

### 7.6. Número máximo de iterações

```python
max_iter = 10000
```

Estabelece um limite de 10.000 iterações.

Esse limite evita que o laço continue indefinidamente caso o critério de convergência não seja atingido.

### 7.7. Contador de iterações

```python
iteracao = 0
```

Inicializa o contador de iterações em zero.

O valor é incrementado depois de cada passagem completa pelos nós internos.

### 7.8. Inicialização do erro

```python
erro = float('inf')
```

Define o erro inicial como infinito.

Isso garante que a condição do `while` seja verdadeira no início, pois infinito é maior que a tolerância.

### 7.9. Laço principal

```python
while erro > tol and iteracao < max_iter:
```

O laço continua enquanto as duas condições forem satisfeitas:

- O erro for maior que a tolerância;
- O número de iterações for menor que o limite máximo.

O processo termina quando a variação máxima fica abaixo ou igual à tolerância, ou quando o limite de iterações é alcançado.

### 7.10. Reinicialização do erro em cada iteração

```python
erro = 0.0
```

No início de cada iteração completa, o erro é reiniciado em zero.

Durante a atualização dos nós, o programa identifica a maior variação absoluta observada nessa iteração.

### 7.11. Percurso das linhas e colunas internas

```python
for j in range(1, n + 1):
    for i in range(1, n + 1):
```

Percorre os nós internos da matriz.

Os índices variam de 1 até `n`, inclusive, pois o limite superior do `range` não é incluído.

Os laços não percorrem os índices zero e `n+1`, preservando as bordas.

A ordem de atualização importa: o laço percorre primeiro os valores de `j` e, para cada um deles, percorre os valores de `i`.

### 7.12. Armazenamento da temperatura anterior

```python
antigo = T[i, j]
```

Guarda a temperatura atual do nó antes de atualizá-lo.

Esse valor será utilizado para calcular a variação de temperatura produzida pela atualização.

### 7.13. Atualização da temperatura

```python
T[i, j] = (
    T[i-1, j] + T[i+1, j]
    + T[i, j-1] + T[i, j+1]
) / 4.0
```

Calcula a média das temperaturas dos quatro vizinhos imediatos.

Em termos matemáticos:

\[
T_{i,j}^{\text{novo}}=
\frac{
T_{i-1,j}+T_{i+1,j}+T_{i,j-1}+T_{i,j+1}
}{4}
\]

O resultado substitui imediatamente a temperatura anterior do nó.

Essa substituição imediata é a característica central do método de Gauss-Seidel: as próximas atualizações podem utilizar os valores novos já calculados na mesma iteração.

### 7.14. Cálculo da variação de temperatura

```python
diferenca = abs(T[i, j] - antigo)
```

Calcula a diferença absoluta entre a temperatura atualizada e a anterior:

\[
\Delta T_{i,j}=
\left|T_{i,j}^{\text{novo}}-T_{i,j}^{\text{antigo}}\right|
\]

O valor absoluto garante que a variação seja não negativa, independentemente de a temperatura ter aumentado ou diminuído.

### 7.15. Identificação da maior variação

```python
erro = max(erro, diferenca)
```

Compara a variação recém-calculada com o maior erro já registrado na iteração.

Ao final dos laços internos, `erro` representa a maior variação absoluta de temperatura observada entre todos os nós atualizados naquela iteração.

Matematicamente:

\[
E^{(k)}=
\max_{i,j}
\left|
T_{i,j}^{(k)}-T_{i,j}^{(k-1)}
\right|
\]

Esse valor é utilizado como critério de convergência.

### 7.16. Atualização do contador

```python
iteracao += 1
```

Incrementa o número de iterações depois que todos os nós internos foram percorridos uma vez.

Portanto, uma iteração corresponde a uma passagem completa pela malha interna, e não à atualização de apenas um nó.

### 7.17. Inicialização do resíduo máximo

```python
max_residuo = 0.0
```

Inicializa a variável que armazenará o maior resíduo absoluto encontrado na matriz.

O resíduo é calculado depois do término do processo iterativo.

### 7.18. Percurso dos nós para calcular o resíduo

```python
for j in range(1, n + 1):
    for i in range(1, n + 1):
```

Percorre novamente todos os nós internos, desta vez para avaliar o quanto cada um satisfaz a equação discreta de equilíbrio.

### 7.19. Cálculo do resíduo local

```python
r_ij = 4 * T[i, j] - (
    T[i-1, j] + T[i+1, j]
    + T[i, j-1] + T[i, j+1]
)
```

Calcula o resíduo da equação discreta no nó \((i,j)\):

\[
r_{i,j}=
4T_{i,j}-
\left(
T_{i-1,j}+T_{i+1,j}
+T_{i,j-1}+T_{i,j+1}
\right)
\]

Se a equação discreta for satisfeita exatamente, o resíduo será zero.

O resíduo mede a diferença entre os dois lados da equação de equilíbrio.

### 7.20. Atualização do resíduo máximo

```python
max_residuo = max(max_residuo, abs(r_ij))
```

Compara o módulo do resíduo do nó atual com o maior módulo de resíduo já encontrado.

Ao final, `max_residuo` representa:

\[
R_{\max}=\max_{i,j}|r_{i,j}|
\]

Esse valor fornece uma medida adicional da qualidade da solução numérica.

### 7.21. Retorno dos resultados

```python
return T, iteracao, erro, max_residuo
```

Retorna quatro informações:

| Variável retornada | Significado |
|---|---|
| `T` | Matriz final de temperaturas |
| `iteracao` | Número de iterações executadas |
| `erro` | Maior variação absoluta observada na última iteração |
| `max_residuo` | Maior resíduo absoluto das equações discretizadas |

O retorno de todos esses valores permite que o programa principal apresente não apenas as temperaturas, mas também informações sobre a convergência.

## 8. Programa principal

Depois de definir a função, o código executa a simulação para diferentes tamanhos de malha.

### 8.1. Lista de tamanhos

```python
tamanhos_malha = [5, 10, 20]
```

Define os valores de `n` utilizados nas simulações.

O objetivo é observar o comportamento do algoritmo quando a quantidade de nós internos aumenta.

### 8.2. Laço de execução das simulações

```python
for n in tamanhos_malha:
```

Percorre cada tamanho de malha definido na lista.

Para cada valor de `n`, o programa realiza uma nova simulação independente.

### 8.3. Identificação da simulação

```python
print(f"--- Simulação para n = {n} ---")
```

Imprime um cabeçalho para identificar o tamanho de malha correspondente aos resultados.

### 8.4. Chamada da função

```python
T, iteracoes, erro_final, max_residuo = resolver_placa_gauss_seidel(n)
```

Executa a função e distribui os quatro valores retornados nas variáveis correspondentes.

### 8.5. Impressão do número de iterações

```python
print(f"Iterações necessárias: {iteracoes}")
```

Apresenta quantas iterações completas foram executadas.

Esse valor pode ser utilizado para comparar o comportamento do algoritmo entre as diferentes malhas.

### 8.6. Impressão do erro final

```python
print(f"Erro final (variação máxima): {erro_final:.2e} °C")
```

Apresenta a variação máxima de temperatura observada na última iteração.

O formato `:.2e` exibe o número em notação científica com duas casas decimais.

Por exemplo, um número como \(0{,}000001\) pode ser apresentado como `1.00e-06`.

### 8.7. Impressão do resíduo máximo

```python
print(f"Resíduo máximo calculado: {max_residuo:.2e}")
```

Apresenta o maior resíduo absoluto calculado entre os nós internos.

Um resíduo menor indica que a equação discreta está sendo satisfeita com maior precisão, embora a interpretação também dependa da escala do problema.

### 8.8. Exibição exclusiva para n = 5

```python
if n == 5:
```

O bloco seguinte é executado somente na simulação cuja malha possui cinco nós internos em cada direção.

Essa condição evita imprimir a matriz completa para todos os tamanhos de malha.

### 8.9. Temperatura do nó central

```python
temp_centro = T[3, 3]
```

Armazena a temperatura do elemento de índices Python `(3, 3)`.

Como o código utiliza índices iniciados em zero e inclui uma borda em cada lado da matriz, esse elemento corresponde ao nó central da malha interna quando \(n=5\).

### 8.10. Impressão da temperatura central

```python
print(f"Temperatura do nó central T(4,4) para n=5: {temp_centro:.4f} °C")
```

Apresenta a temperatura central com quatro casas decimais.

A indicação `T(4,4)` usa a convenção de índices iniciados em um descrita no comentário do programa; no Python, o elemento acessado é `T[3,3]`.

### 8.11. Impressão da matriz reorganizada

```python
print(np.flipud(T.T))
```

Aplica duas operações à matriz:

1. `T.T`: transpõe a matriz, trocando linhas por colunas;
2. `np.flipud(...)`: inverte a ordem das linhas da matriz transposta.

A combinação é utilizada para apresentar a matriz com a borda a \(100\,^\circ\mathrm{C}\) visualmente na parte inferior.

Essa operação altera apenas a apresentação dos dados impressos. A matriz original `T` retornada pela função permanece na sua orientação original.

## 9. Entendendo os critérios de convergência

O programa utiliza dois indicadores diferentes: a variação máxima entre atualizações e o resíduo máximo da equação discreta.

Eles não representam a mesma coisa.

### 9.1. Erro de atualização

O erro utilizado no laço principal é:

\[
E^{(k)}=
\max_{i,j}
\left|T_{i,j}^{(k)}-T_{i,j}^{(k-1)}\right|
\]

Ele mede a maior mudança de temperatura em uma passagem completa pela malha.

O critério de parada é:

\[
E^{(k)}\leq 10^{-6}
\]

O laço continua enquanto `erro > tol`.

Esse indicador mede a mudança entre iterações, não a diferença exata entre a solução numérica e uma solução verdadeira desconhecida.

### 9.2. Resíduo máximo

O resíduo é calculado por:

\[
r_{i,j}=4T_{i,j}-
\left(
T_{i-1,j}+T_{i+1,j}
+T_{i,j-1}+T_{i,j+1}
\right)
\]

O maior resíduo absoluto é:

\[
R_{\max}=\max_{i,j}|r_{i,j}|
\]

Esse indicador mede o quanto a solução final satisfaz as equações discretizadas.

O resíduo não é utilizado no `while` para determinar quando o algoritmo deve parar. Ele é calculado apenas depois que o laço iterativo termina.

### 9.3. Por que os valores podem ser diferentes?

O erro de atualização e o resíduo medem aspectos distintos do processo numérico.

- O erro de atualização mede quanto os valores mudam entre iterações;
- O resíduo mede quanto a equação discreta é satisfeita na solução final.

Assim, não se deve interpretar automaticamente o erro final como sendo igual ao resíduo máximo.

### 9.4. O que acontece se o algoritmo atingir 10.000 iterações?

A condição do laço também verifica:

```python
iteracao < max_iter
```

Se o contador alcançar 10.000, o laço termina mesmo que o erro ainda esteja acima da tolerância.

A função, no entanto, retorna os resultados sem indicar explicitamente se a convergência foi alcançada.

Por isso, ao analisar os resultados, é importante conferir o valor de `erro_final` e verificar se ele é menor ou igual a `tol`.

## 10. Comparação entre os tamanhos de malha

O programa executa simulações com \(n=5\), \(n=10\) e \(n=20\).

| Parâmetro \(n\) | Dimensão da matriz | Nós internos atualizados |
|---:|---:|---:|
| 5 | \(7\times7\) | 25 |
| 10 | \(12\times12\) | 100 |
| 20 | \(22\times22\) | 400 |

Quando o valor de `n` aumenta, cresce o número de nós internos e, consequentemente, a quantidade de atualizações realizadas em cada iteração completa.

Uma malha maior representa a distribuição por meio de mais pontos internos. Entretanto, isso não significa automaticamente que o resultado seja mais preciso em relação ao problema físico contínuo: essa precisão também depende do refinamento da malha, do modelo adotado e das condições de contorno.

O código apresentado não calcula explicitamente o erro de discretização em relação a uma solução analítica.

## 11. Complexidade computacional

A complexidade do algoritmo depende do tamanho da malha e do número de iterações.

Em cada iteração, o programa percorre \(n^2\) nós internos. Para cada nó, realiza um número constante de operações aritméticas.

Portanto, o custo de uma iteração completa é:

\[
O(n^2)
\]

Se o algoritmo executar \(K\) iterações, o custo total da etapa iterativa é aproximadamente:

\[
O(Kn^2)
\]

O cálculo do resíduo percorre novamente os \(n^2\) nós, acrescentando um custo de ordem \(O(n^2)\).

O armazenamento da matriz de temperaturas exige espaço de ordem:

\[
O(n^2)
\]

Assim, o custo total depende tanto do tamanho da malha quanto do número de iterações necessárias para satisfazer o critério de parada.

## 12. Limitações e cuidados

### 12.1. Validação do parâmetro n

A função pressupõe que `n` seja um inteiro positivo. O código original não verifica explicitamente esse requisito.

### 12.2. Condições de contorno fixas

A temperatura de \(100\,^\circ\mathrm{C}\) e a temperatura inicial de \(20\,^\circ\mathrm{C}\) estão definidas diretamente no código.

Para simular outros cenários, seria necessário modificar esses valores ou tornar as condições de contorno parâmetros da função.

### 12.3. Espaçamento da malha

A fórmula utilizada pressupõe espaçamentos iguais nas duas direções.

Para espaçamentos diferentes, a equação discretizada precisaria ser adaptada.

### 12.4. Critério de parada

O programa usa a variação máxima de temperatura para interromper o processo. Ele não exige simultaneamente que o resíduo máximo seja inferior a uma tolerância específica.

### 12.5. Ausência de indicador explícito de convergência

O programa retorna o erro e o número de iterações, mas não retorna uma variável booleana que indique diretamente se o critério foi satisfeito antes do limite máximo.

### 12.6. Interpretação física

A simulação representa um modelo discreto de condução de calor. Para interpretar os resultados como um problema físico real, é necessário verificar se as hipóteses do modelo correspondem ao sistema estudado.

## 13. Resumo das variáveis principais

| Variável | Significado |
|---|---|
| `np` | Apelido utilizado para a biblioteca NumPy |
| `n` | Número de nós internos em cada direção |
| `T` | Matriz que armazena as temperaturas |
| `tol` | Tolerância para o critério de convergência |
| `max_iter` | Número máximo de iterações |
| `iteracao` | Contador de iterações completas |
| `erro` | Maior variação de temperatura na iteração atual |
| `antigo` | Temperatura do nó antes da atualização |
| `diferenca` | Variação absoluta da temperatura de um nó |
| `h` | Não é utilizada no código desta simulação |
| `max_residuo` | Maior resíduo absoluto entre os nós internos |
| `r_ij` | Resíduo da equação discreta em um nó |
| `temp_centro` | Temperatura armazenada em `T[3,3]` para \(n=5\) |
| `tamanhos_malha` | Lista de tamanhos de malha simulados |
| `iteracoes` | Número de iterações retornado pela função |
| `erro_final` | Erro de atualização retornado pela função |

## 14. Como executar o programa

### Requisitos

- Python 3;
- Biblioteca NumPy.

Para instalar o NumPy, se necessário, utilize:

```bash
pip install numpy
```

Salve o código em um arquivo Python, por exemplo:

```text
PPC surpresa.py
```

Execute o arquivo no terminal:

```bash
python "PPC surpresa.py"
```

O programa imprimirá os resultados das três simulações e a matriz de temperaturas para \(n=5\).

## 15. Conclusão

O programa demonstra a aplicação do método de Gauss-Seidel na solução iterativa de um problema bidimensional de distribuição de temperatura.

A malha é representada por uma matriz NumPy, e as temperaturas dos nós internos são atualizadas a partir da média dos quatro vizinhos. O processo continua até que a variação máxima fique abaixo da tolerância ou que o limite de iterações seja atingido.

Além das temperaturas calculadas, o código fornece o número de iterações, o erro de atualização e o resíduo máximo. Essas informações permitem avaliar o comportamento numérico da simulação e comparar diferentes tamanhos de malha.

Para análises futuras, melhorias úteis incluem tornar as condições de contorno configuráveis, validar o parâmetro `n`, informar explicitamente se a convergência foi alcançada e comparar os resultados com uma solução analítica ou outra referência numérica quando disponível.
