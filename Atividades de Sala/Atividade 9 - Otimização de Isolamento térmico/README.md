# Otimização da Espessura de Isolamento Térmico em Tubulações

## 1. Visão geral do projeto

Este projeto implementa um modelo numérico para determinar a **espessura ótima de isolamento térmico de uma tubulação**, considerando simultaneamente a perda de calor, o custo da energia térmica reposta e o custo de instalação do material isolante.

O programa combina conceitos de transferência de calor, análise econômica e métodos numéricos de otimização.

São utilizados dois métodos para encontrar a espessura que minimiza o custo total:

1. **Método da Razão Áurea:** procura o mínimo de uma função em um intervalo, reduzindo progressivamente a região de busca.
2. **Método da Interpolação Quadrática:** constrói uma parábola a partir de três pontos da função e utiliza o vértice dessa parábola como nova aproximação do mínimo.

Além de comparar os métodos, o programa:

- Calcula a taxa de perda de calor por unidade de comprimento;
- Estima o custo acumulado da energia térmica perdida;
- Calcula o custo de instalação do isolamento;
- Avalia o custo total para diferentes espessuras;
- Gera um gráfico do custo total em função da espessura;
- Compara os resultados dos dois métodos de otimização;
- Investiga a influência do preço da energia na espessura ótima;
- Calcula o raio crítico de isolamento térmico para uma geometria cilíndrica.

### Objetivo principal

Encontrar a espessura \(e\) que minimiza o custo total:

\[
\boxed{\min_{e\in[0,\;0{,}12]} C(e)}
\]

em que \(C(e)\) representa a soma do custo da energia térmica perdida durante o período de operação e do custo de instalação do isolamento.

O intervalo de busca corresponde a espessuras entre \(0\) e \(0{,}12\ \text{m}\), ou seja, entre 0 e 120 mm.

---

## 2. Organização do código

O programa pode ser dividido em seis partes principais.

| Parte | Funções ou comandos | Objetivo |
|---|---|---|
| 1. Dados do problema | Variáveis físicas e econômicas | Definir as condições da tubulação |
| 2. Modelo térmico | `calcular_q_linha(e)` | Calcular a perda de calor |
| 3. Modelo econômico | `calcular_custos(e, pE)` e `funcao_objetivo(e, pE)` | Calcular os custos e a função a minimizar |
| 4. Otimização | `otimizar_razao_aurea(...)` | Encontrar o mínimo pelo método da Razão Áurea |
| 5. Otimização | `otimizar_interp_quadratica(...)` | Encontrar o mínimo por interpolação quadrática |
| 6. Execução e análise | Bloco `if __name__ == "__main__":` | Gerar tabelas, comparar métodos e produzir gráficos |

A separação em funções facilita a consulta, pois cada parte desempenha uma tarefa específica.

---

## 3. Bibliotecas utilizadas

O código utiliza duas bibliotecas Python.

### 3.1. NumPy

```python
import numpy as np
```

O NumPy é utilizado em operações matemáticas e numéricas, incluindo:

- Cálculo de raízes quadradas;
- Cálculo do logaritmo natural;
- Geração de valores igualmente espaçados para os gráficos;
- Avaliação de valores absolutos;
- Criação de vetores numéricos.

Exemplos presentes no código:

```python
np.sqrt(5)
np.log(r2 / r1)
np.linspace(0, 0.12, 100)
np.abs(x_novo - x2)
```

### 3.2. Matplotlib

```python
import matplotlib.pyplot as plt
```

O Matplotlib é utilizado para gerar e salvar os gráficos da análise econômica.

O programa cria dois gráficos:

- `grafico_custo_e.png`;
- `grafico_e_pE.png`.

Os gráficos são salvos em arquivos de imagem no diretório de execução do programa.

### 3.3. Instalação das bibliotecas

Se as bibliotecas ainda não estiverem instaladas, utilize:

```bash
pip install numpy matplotlib
```

---

## 4. Dados físicos e econômicos

Os parâmetros do problema são definidos no início do código.

```python
r1 = 0.040
L = 20.0
Ti = 180.0
T_inf = 40.0
k = 0.045
h = 12.0
t = 30000.0
pE_padrao = 0.70
cI = 9000.0
```

Cada variável representa uma característica física ou econômica do sistema.

### 4.1. Tabela de parâmetros

| Variável | Valor | Unidade | Significado |
|---|---:|---|---|
| `r1` | 0,040 | m | Raio externo da tubulação sem isolamento |
| `L` | 20,0 | m | Comprimento da tubulação |
| `Ti` | 180,0 | °C | Temperatura na interface da tubulação com o isolamento |
| `T_inf` | 40,0 | °C | Temperatura ambiente |
| `k` | 0,045 | W/(m·K) | Condutividade térmica do isolante |
| `h` | 12,0 | W/(m²·K) | Coeficiente de convecção externa |
| `t` | 30.000 | h | Tempo acumulado de operação |
| `pE_padrao` | 0,70 | R$/kWh | Preço de referência da energia térmica reposta |
| `cI` | 9.000 | R$/m³ | Custo instalado do material isolante por volume |

### 4.2. Interpretação dos parâmetros

**Raio da tubulação (`r1`)**

O raio inicial da tubulação é:

\[
r_1=0{,}040\ \text{m}=40\ \text{mm}
\]

Esse valor representa o raio externo antes da aplicação do isolamento.

**Comprimento (`L`)**

O comprimento considerado é de 20 metros. Ele influencia tanto a perda total de calor quanto o volume de material isolante necessário.

**Temperaturas (`Ti` e `T_inf`)**

A diferença de temperatura é:

\[
\Delta T=T_i-T_\infty
\]

Substituindo os valores:

\[
\Delta T=180-40
\]

\[
\boxed{\Delta T=140\ ^\circ\text{C}}
\]

Para diferenças de temperatura, o valor numérico em graus Celsius é igual ao valor numérico em kelvin.

**Condutividade térmica (`k`)**

A condutividade representa a capacidade do material isolante de conduzir calor. No modelo utilizado, ela é considerada constante.

**Coeficiente de convecção (`h`)**

Representa a intensidade da transferência de calor entre a superfície externa do isolamento e o ambiente.

**Tempo de operação (`t`)**

O tempo acumulado é de 30.000 horas. Esse valor é utilizado para estimar o custo acumulado da energia térmica perdida.

**Preço da energia (`pE_padrao`)**

O preço de referência é de R$ 0,70 por kWh de energia térmica reposta.

**Custo de instalação (`cI`)**

O custo instalado é de R$ 9.000 por metro cúbico de isolante, considerando o modelo econômico adotado no código.

---

## 5. Geometria da tubulação com isolamento

A espessura do isolamento é representada por \(e\).

Quando o isolamento é aplicado, o raio externo passa a ser:

\[
\boxed{r_2=r_1+e}
\]

Em que:

- \(r_1\) é o raio externo original da tubulação;
- \(e\) é a espessura do isolamento;
- \(r_2\) é o raio externo após o isolamento.

Por exemplo, se:

\[
r_1=0{,}040\ \text{m}
\]

e:

\[
e=0{,}020\ \text{m}
\]

então:

\[
r_2=0{,}040+0{,}020
\]

\[
\boxed{r_2=0{,}060\ \text{m}}
\]

O raio externo passa a ser 60 mm.

Essa relação é utilizada tanto no cálculo da resistência térmica quanto no cálculo do volume do isolante.

---

## 6. Modelo físico: resistência térmica

O modelo de perda de calor considera duas resistências térmicas em série:

1. Resistência à condução através do isolamento;
2. Resistência à convecção entre a superfície externa e o ambiente.

A resistência total utilizada é:

\[
R_{\text{total}}=R_{\text{cond}}+R_{\text{conv}}
\]

A taxa de perda de calor por unidade de comprimento é calculada por:

\[
\boxed{
q'(e)=\frac{T_i-T_\infty}
{R_{\text{cond}}+R_{\text{conv}}}
}
\]

A unidade de \(q'\) é W/m.

### 6.1. Resistência térmica por condução

Para uma camada cilíndrica de isolamento, a resistência térmica de condução por unidade de comprimento é:

\[
\boxed{
R'_{\text{cond}}=
\frac{\ln(r_2/r_1)}{2\pi k}
}
\]

A expressão utiliza o logaritmo natural da razão entre os raios.

No código:

```python
R_cond = np.log(r2 / r1) / (2.0 * np.pi * k)
```

A resistência depende da espessura do isolamento porque o raio externo \(r_2\) varia com \(e\).

Quanto maior o raio externo, maior tende a ser a resistência à condução.

**Observação sobre unidades:** essa expressão representa a resistência térmica por unidade de comprimento, com unidade equivalente a m·K/W. Ela deve ser utilizada em conjunto com a taxa de calor por unidade de comprimento, em W/m.

### 6.2. Caso sem isolamento

Quando \(e=0\):

\[
r_2=r_1
\]

Logo:

\[
\ln\left(\frac{r_2}{r_1}\right)=\ln(1)=0
\]

Portanto:

\[
R'_{\text{cond}}=0
\]

O código trata explicitamente esse caso:

```python
if e == 0.0:
    R_cond = 0.0
else:
    R_cond = np.log(r2 / r1) / (2.0 * np.pi * k)
```

A condição evita calcular desnecessariamente a expressão de condução quando a espessura é exatamente zero.

### 6.3. Resistência térmica por convecção

A resistência térmica de convecção externa por unidade de comprimento é:

\[
\boxed{
R'_{\text{conv}}=\frac{1}{2\pi h r_2}
}
\]

No código:

```python
R_conv = 1.0 / (2.0 * np.pi * h * r2)
```

Essa resistência diminui quando o raio externo aumenta, pois a área externa de transferência de calor aumenta.

A área lateral de um cilindro é proporcional a:

\[
A=2\pi r_2L
\]

Assim, um raio maior representa uma área externa maior para a troca de calor com o ambiente.

### 6.4. Taxa de perda de calor

A expressão utilizada pelo programa é:

\[
q'(e)=
\frac{T_i-T_\infty}
{
\displaystyle
\frac{\ln(r_2/r_1)}{2\pi k}
+
\frac{1}{2\pi h r_2}
}
\]

Para \(e=0\), a resistência de condução é zero e o cálculo considera somente a resistência de convecção.

A taxa calculada é uma perda de calor por unidade de comprimento. Para obter a potência térmica total da tubulação, multiplica-se por seu comprimento:

\[
\boxed{\dot Q=q'(e)L}
\]

A unidade de \(\dot Q\) é W.

---

## 7. Função `calcular_q_linha(e)`

Essa função calcula a taxa de perda de calor por unidade de comprimento para uma espessura de isolamento específica.

### Código

```python
def calcular_q_linha(e):
    """
    Calcula a taxa de perda de calor q'(e) em W/m.
    """
    r2 = r1 + e

    if e == 0.0:
        R_cond = 0.0
    else:
        R_cond = np.log(r2 / r1) / (2.0 * np.pi * k)

    R_conv = 1.0 / (2.0 * np.pi * h * r2)
    q_linha = (Ti - T_inf) / (R_cond + R_conv)

    return q_linha
```

### Explicação linha por linha

**1. Definição da função**

```python
def calcular_q_linha(e):
```

Recebe a espessura do isolamento, em metros.

**2. Cálculo do raio externo**

```python
r2 = r1 + e
```

Determina o raio externo da tubulação isolada.

**3. Tratamento do caso sem isolamento**

```python
if e == 0.0:
    R_cond = 0.0
```

Se a espessura for zero, a resistência de condução da camada isolante é definida como zero.

**4. Cálculo da resistência de condução**

```python
else:
    R_cond = np.log(r2 / r1) / (2.0 * np.pi * k)
```

Para espessuras positivas, calcula a resistência térmica cilíndrica do isolamento.

**5. Cálculo da resistência de convecção**

```python
R_conv = 1.0 / (2.0 * np.pi * h * r2)
```

Calcula a resistência à transferência de calor entre a superfície externa e o ambiente.

**6. Cálculo da taxa de perda de calor**

```python
q_linha = (Ti - T_inf) / (R_cond + R_conv)
```

Divide a diferença de temperatura pela soma das resistências térmicas por unidade de comprimento.

**7. Retorno**

```python
return q_linha
```

Retorna o valor de \(q'\), em W/m.

### Utilização isolada da função

```python
e = 0.02
q_linha = calcular_q_linha(e)

print(f"Perda de calor por unidade de comprimento: {q_linha:.4f} W/m")
```

Esse exemplo permite calcular a perda de calor para uma espessura de 20 mm, utilizando os parâmetros definidos no início do programa.

---

## 8. Modelo econômico: custo da energia

O custo da energia térmica perdida é calculado considerando a taxa de perda de calor, o comprimento da tubulação, o tempo de operação e o preço da energia.

A potência térmica perdida é:

\[
\dot Q=q'(e)L
\]

Como \(q'\) está em W/m e \(L\) em m, \(\dot Q\) fica em W.

A energia acumulada durante \(t\) horas, em kWh, é:

\[
E_{\text{perdida}}=
\frac{q'(e)Lt}{1000}
\]

O fator 1000 é utilizado para converter Wh em kWh.

O custo de energia é:

\[
\boxed{
C_E(e)=\frac{p_E L q'(e)t}{1000}
}
\]

Em que:

- \(p_E\) é o preço da energia, em R$/kWh;
- \(L\) é o comprimento da tubulação, em m;
- \(q'(e)\) é a taxa de perda de calor, em W/m;
- \(t\) é o tempo de operação, em h.

O resultado é expresso em reais.

### Interpretação

Se a perda de calor diminui, a energia térmica necessária para repor essa perda também diminui.

Consequentemente, o custo de energia tende a diminuir quando a espessura do isolamento aumenta, desde que o modelo térmico apresente redução da perda de calor nesse intervalo.

---

## 9. Modelo econômico: custo de instalação do isolamento

O custo de instalação depende do volume de material isolante aplicado à tubulação.

O volume de uma camada cilíndrica é calculado pela diferença entre o volume do cilindro externo e o volume ocupado pela tubulação original:

\[
V_{\text{isolante}}=
\pi r_2^2L-\pi r_1^2L
\]

Colocando os termos em evidência:

\[
\boxed{
V_{\text{isolante}}=
\pi L(r_2^2-r_1^2)
}
\]

Como:

\[
r_2=r_1+e
\]

também podemos escrever:

\[
V_{\text{isolante}}=
\pi L\left[(r_1+e)^2-r_1^2\right]
\]

O custo de instalação é:

\[
\boxed{
C_I(e)=c_I\pi L\left[(r_1+e)^2-r_1^2\right]
}
\]

Em que \(c_I\) representa o custo instalado por unidade de volume, em R$/m³.

O resultado é expresso em reais.

### Interpretação

Quanto maior a espessura, maior o volume de material necessário. Portanto, o custo de instalação cresce com a espessura.

Esse crescimento cria um compromisso econômico: aumentar o isolamento pode reduzir a perda de energia, mas exige maior investimento inicial.

---

## 10. Função `calcular_custos(e, pE)`

Essa função reúne os cálculos de custo de energia e de instalação, retornando os dois componentes e sua soma.

### Código

```python
def calcular_custos(e, pE):
    """
    Retorna os componentes de custo (CE, CI) e o Custo Total C(e) em reais.
    """
    q_linha = calcular_q_linha(e)

    # Custo de energia
    CE = (pE * L * q_linha * t) / 1000.0

    # Custo de instalação
    CI = cI * np.pi * L * ((r1 + e)**2 - r1**2)

    return CE, CI, CE + CI
```

### Explicação linha por linha

**1. Definição da função**

```python
def calcular_custos(e, pE):
```

Recebe dois parâmetros:

- `e`: espessura do isolamento, em metros;
- `pE`: preço da energia térmica, em R$/kWh.

**2. Cálculo da perda de calor**

```python
q_linha = calcular_q_linha(e)
```

Chama a função térmica para obter a perda de calor por unidade de comprimento.

**3. Cálculo do custo de energia**

```python
CE = (pE * L * q_linha * t) / 1000.0
```

Aplica a expressão:

\[
C_E=\frac{p_E L q't}{1000}
\]

**4. Cálculo do custo de instalação**

```python
CI = cI * np.pi * L * ((r1 + e)**2 - r1**2)
```

Calcula o volume da camada de isolamento e multiplica pelo custo instalado por metro cúbico.

**5. Retorno dos três valores**

```python
return CE, CI, CE + CI
```

Retorna uma tupla com:

1. `CE`: custo de energia;
2. `CI`: custo de instalação;
3. `CE + CI`: custo total.

### Exemplo de utilização

```python
e = 0.02
pE = 0.70

CE, CI, C_total = calcular_custos(e, pE)

print(f"Custo de energia: R$ {CE:.2f}")
print(f"Custo de instalação: R$ {CI:.2f}")
print(f"Custo total: R$ {C_total:.2f}")
```

Esse exemplo calcula os três valores para uma espessura de 20 mm e preço de energia de R$ 0,70/kWh.

---

## 11. Função objetivo `funcao_objetivo(e, pE)`

Os métodos de otimização precisam de uma função numérica que possam avaliar em diferentes valores de espessura.

Essa função é chamada de **função objetivo**.

No projeto, a função objetivo é o custo total:

\[
\boxed{C(e)=C_E(e)+C_I(e)}
\]

### Código

```python
def funcao_objetivo(e, pE):
    """
    Retorna estritamente o Custo Total, servindo como nossa f(x) a minimizar.
    """
    return calcular_custos(e, pE)[2]
```

O índice `[2]` seleciona o terceiro elemento da tupla retornada por `calcular_custos`, isto é, o custo total.

A função objetivo não retorna separadamente os custos de energia e instalação. Ela retorna apenas o valor que os métodos de otimização devem minimizar.

### Exemplo

```python
custo = funcao_objetivo(0.02, 0.70)
print(custo)
```

Esse comando avalia o custo total para \(e=0{,}02\ \text{m}\) e \(p_E=0{,}70\ \text{R$/kWh}\).

---

## 12. Método da Razão Áurea

O método da Razão Áurea é um método de otimização unidimensional utilizado para localizar o mínimo de uma função em um intervalo.

Ele não necessita do cálculo da derivada da função objetivo. Em vez disso, compara os valores da função em dois pontos internos e elimina progressivamente a região que não contém o mínimo procurado, sob a hipótese de unimodalidade no intervalo.

### 12.1. Constante da Razão Áurea

O código calcula:

```python
tau = (np.sqrt(5) - 1.0) / 2.0
```

Matematicamente:

\[
\boxed{
\tau=\frac{\sqrt{5}-1}{2}
\approx0{,}618034
}
\]

Essa constante permite definir dois pontos internos do intervalo de busca.

### 12.2. Cálculo dos pontos internos

Para um intervalo \([a,b]\), o código utiliza:

\[
x_1=b-\tau(b-a)
\]

\[
x_2=a+\tau(b-a)
\]

Os pontos estão dentro do intervalo de busca.

Calculam-se os valores da função:

\[
f_1=C(x_1)
\]

\[
f_2=C(x_2)
\]

Esses valores são comparados para decidir qual parte do intervalo pode ser descartada.

### 12.3. Regra de atualização

Para minimizar uma função unimodal:

- Se \(f_1<f_2\), o mínimo está na região preservada ao descartar a parte direita do intervalo;
- Caso contrário, descarta-se a parte esquerda.

A comparação permite reduzir o intervalo de busca sem precisar calcular derivadas.

### 12.4. Critério de parada

O código verifica:

```python
if (b - a) < tol:
```

Quando a largura do intervalo fica menor que a tolerância, o método retorna o ponto médio:

\[
\boxed{x_{\text{opt}}=\frac{a+b}{2}}
\]

A tolerância padrão é:

\[
\text{tol}=10^{-6}
\]

Como a variável otimizada é a espessura, a tolerância está em metros.

### 12.5. Implementação completa

```python
def otimizar_razao_aurea(a, b, pE, tol=1e-6, max_iter=200):
    tau = (np.sqrt(5) - 1.0) / 2.0

    x1 = b - tau * (b - a)
    x2 = a + tau * (b - a)

    f1 = funcao_objetivo(x1, pE)
    f2 = funcao_objetivo(x2, pE)
    avaliacoes = 2

    for k_iter in range(1, max_iter + 1):
        if (b - a) < tol:
            x_opt = (a + b) / 2.0
            return (
                x_opt,
                funcao_objetivo(x_opt, pE),
                k_iter,
                avaliacoes,
                "Convergiu"
            )

        if f1 < f2:
            b = x2
            x2 = x1
            f2 = f1

            x1 = b - tau * (b - a)
            f1 = funcao_objetivo(x1, pE)

        else:
            a = x1
            x1 = x2
            f1 = f2

            x2 = a + tau * (b - a)
            f2 = funcao_objetivo(x2, pE)

        avaliacoes += 1

    x_opt = (a + b) / 2.0

    return (
        x_opt,
        funcao_objetivo(x_opt, pE),
        max_iter,
        avaliacoes,
        "Máx iterações"
    )
```

### 12.6. Explicação dos parâmetros

| Parâmetro | Significado |
|---|---|
| `a` | Limite inferior do intervalo |
| `b` | Limite superior do intervalo |
| `pE` | Preço da energia |
| `tol` | Tolerância para a largura do intervalo |
| `max_iter` | Número máximo de iterações |

### 12.7. Valores retornados

A função retorna cinco elementos:

1. `x_opt`: espessura ótima estimada;
2. Custo total correspondente;
3. Número de iterações;
4. Número de avaliações contabilizadas da função objetivo;
5. Situação final do método.

A situação pode ser `"Convergiu"` ou `"Máx iterações"`.

**Observação de implementação:** o contador `avaliacoes` não inclui a avaliação adicional de `funcao_objetivo(x_opt, pE)` realizada no retorno. Portanto, o número apresentado não contabiliza todas as chamadas da função objetivo.

### 12.8. Exemplo de utilização

```python
e_otimo, custo_min, iteracoes, avaliacoes, status = (
    otimizar_razao_aurea(
        0.0,
        0.12,
        0.70
    )
)

print(f"Espessura ótima: {e_otimo:.6f} m")
print(f"Custo mínimo: R$ {custo_min:.2f}")
print(f"Iterações: {iteracoes}")
print(f"Avaliações contabilizadas: {avaliacoes}")
print(f"Situação: {status}")
```

O método procura o menor custo no intervalo entre 0 e 0,12 m.

---

## 13. Método da Interpolação Quadrática

O método da Interpolação Quadrática utiliza três pontos para construir uma parábola que aproxima localmente a função objetivo.

A ideia é utilizar o vértice dessa parábola como uma nova estimativa do ponto de mínimo.

Diferentemente da Razão Áurea, esse método utiliza informações sobre os valores da função em três pontos para estimar diretamente a localização do mínimo.

### 13.1. Escolha dos pontos iniciais

O código inicia com:

```python
otimizar_interp_quadratica(0.04, 0.07, 0.11, pE_padrao)
```

Portanto:

\[
x_1=0{,}04\ \text{m}
\]

\[
x_2=0{,}07\ \text{m}
\]

\[
x_3=0{,}11\ \text{m}
\]

Os pontos devem ser distintos e, para a lógica de atualização implementada, devem estar ordenados:

\[
x_1<x_2<x_3
\]

O código fornecido não realiza uma validação completa dessas condições.

### 13.2. Valores da função

Calculam-se:

\[
f_1=C(x_1)
\]

\[
f_2=C(x_2)
\]

\[
f_3=C(x_3)
\]

Esses valores são utilizados para determinar a parábola interpoladora.

### 13.3. Primeira diferença dividida

O código calcula:

```python
d1 = (f2 - f1) / (x2 - x1)
```

Matematicamente:

\[
d_1=\frac{f_2-f_1}{x_2-x_1}
\]

Essa expressão representa a inclinação da reta que passa pelos pontos \((x_1,f_1)\) e \((x_2,f_2)\).

### 13.4. Coeficiente quadrático

O coeficiente quadrático é calculado por:

```python
A = ((f3 - f2) / (x3 - x2) - d1) / (x3 - x1)
```

Matematicamente:

\[
\boxed{
A=
\frac{
\displaystyle\frac{f_3-f_2}{x_3-x_2}
-
\displaystyle\frac{f_2-f_1}{x_2-x_1}
}{
x_3-x_1
}
}
\]

Esse coeficiente representa a curvatura da parábola interpoladora.

### 13.5. Interpretação do coeficiente A

Uma parábola pode ser escrita como:

\[
p(x)=Ax^2+Bx+C
\]

Para possuir um mínimo, ela precisa abrir para cima, o que ocorre quando:

\[
A>0
\]

Se \(A<0\), a parábola abre para baixo e seu vértice corresponde a um máximo, não a um mínimo.

No código, existe uma proteção:

```python
if A <= 1e-12:
    return x2, f2, k_iter, avaliacoes, "Proteção ativada (A <= 0)"
```

Essa condição encerra o método quando o coeficiente é não positivo ou muito pequeno. Nesse caso, o código retorna o ponto atual \(x_2\) e seu custo, sem continuar a interpolação.

### 13.6. Coeficiente linear

O código calcula:

```python
B = d1 - A * (x1 + x2)
```

A expressão decorre da representação quadrática e das condições de interpolação nos pontos iniciais.

O coeficiente constante \(C\) não precisa ser calculado explicitamente para encontrar o vértice da parábola.

### 13.7. Vértice da parábola

Para:

\[
p(x)=Ax^2+Bx+C
\]

o vértice ocorre em:

\[
\boxed{x_v=-\frac{B}{2A}}
\]

O código utiliza:

```python
x_novo = -B / (2.0 * A)
```

Esse valor é a nova aproximação candidata para o mínimo.

### 13.8. Proteção contra pontos fora do intervalo

O código verifica:

```python
if x_novo <= x1 or x_novo >= x3 or np.abs(x_novo - x2) < 1e-10:
    return x2, f2, k_iter, avaliacoes, "Proteção (Fora dos limites/Coincidente)"
```

A condição impede a continuação quando:

- A nova aproximação está fora do intervalo delimitado por \(x_1\) e \(x_3\), incluindo os limites;
- A nova aproximação está praticamente coincidente com \(x_2\).

Essa proteção evita algumas situações em que o método não consegue produzir progresso útil.

Se a proteção for ativada, a função retorna o ponto atual \(x_2\), e não a nova aproximação rejeitada.

### 13.9. Avaliação do novo ponto

```python
f_novo = funcao_objetivo(x_novo, pE)
avaliacoes += 1
```

Calcula o custo total na nova espessura e incrementa o contador de avaliações da função objetivo.

### 13.10. Atualização dos três pontos

O método compara \(x_{\text{novo}}\) com \(x_2\) e compara o custo correspondente com \(f_2\).

Se o novo ponto estiver à esquerda de \(x_2\):

- Quando \(f_{\text{novo}}<f_2\), o novo ponto passa a ser o melhor ponto central;
- Caso contrário, o ponto novo substitui o limite inferior.

Se o novo ponto estiver à direita de \(x_2\):

- Quando \(f_{\text{novo}}<f_2\), o novo ponto passa a ser o melhor ponto central;
- Caso contrário, o ponto novo substitui o limite superior.

Essas regras mantêm três pontos que delimitam a região de busca e preservam o melhor ponto conhecido segundo as comparações realizadas.

### 13.11. Implementação completa

```python
def otimizar_interp_quadratica(
    x1, x2, x3, pE, tol=1e-6, max_iter=200
):
    f1 = funcao_objetivo(x1, pE)
    f2 = funcao_objetivo(x2, pE)
    f3 = funcao_objetivo(x3, pE)
    avaliacoes = 3

    for k_iter in range(1, max_iter + 1):
        if (x3 - x1) < tol:
            return x2, f2, k_iter, avaliacoes, "Convergiu"

        d1 = (f2 - f1) / (x2 - x1)
        A = ((f3 - f2) / (x3 - x2) - d1) / (x3 - x1)

        if A <= 1e-12:
            return x2, f2, k_iter, avaliacoes, "Proteção ativada (A <= 0)"

        B = d1 - A * (x1 + x2)
        x_novo = -B / (2.0 * A)

        if (
            x_novo <= x1
            or x_novo >= x3
            or np.abs(x_novo - x2) < 1e-10
        ):
            return (
                x2,
                f2,
                k_iter,
                avaliacoes,
                "Proteção (Fora dos limites/Coincidente)"
            )

        f_novo = funcao_objetivo(x_novo, pE)
        avaliacoes += 1

        if x_novo < x2:
            if f_novo < f2:
                x3, f3 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x1, f1 = x_novo, f_novo
        else:
            if f_novo < f2:
                x1, f1 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x3, f3 = x_novo, f_novo

    return x2, f2, max_iter, avaliacoes, "Máx iterações"
```

### 13.12. Parâmetros da função

| Parâmetro | Significado |
|---|---|
| `x1` | Primeiro ponto inicial |
| `x2` | Segundo ponto inicial |
| `x3` | Terceiro ponto inicial |
| `pE` | Preço da energia |
| `tol` | Tolerância para a largura do intervalo |
| `max_iter` | Número máximo de iterações |

### 13.13. Valores retornados

A função retorna:

1. Espessura estimada;
2. Custo correspondente;
3. Número de iterações;
4. Número de avaliações contabilizadas;
5. Situação final.

A situação pode indicar convergência, ativação de uma proteção numérica ou alcance do número máximo de iterações.

**Importante:** o retorno da função não significa necessariamente que o método tenha convergido. Se uma proteção for ativada, o programa retorna o melhor ponto central disponível segundo a implementação, mas a espessura retornada não está necessariamente dentro da tolerância especificada.

---

## 14. Comparação dos métodos de otimização

Os dois métodos procuram minimizar a mesma função objetivo, mas utilizam estratégias diferentes.

| Característica | Razão Áurea | Interpolação Quadrática |
|---|---|---|
| Tipo | Redução de intervalo | Interpolação parabólica |
| Pontos iniciais | Dois limites | Três pontos |
| Utiliza derivadas? | Não | Não |
| Estratégia | Descarta partes do intervalo | Estima o vértice de uma parábola |
| Principal vantagem | Busca sistemática e robusta para funções unimodais | Pode aproximar o mínimo rapidamente |
| Principal limitação | Pode exigir muitas avaliações | Sensível à escolha dos pontos e à qualidade da parábola |
| Proteções numéricas no código | Limite de iterações e tolerância | Limite de iterações, curvatura e limites do intervalo |

### Qual método é melhor?

Não existe uma resposta universal.

A Razão Áurea costuma ser uma escolha confiável quando a função é unimodal no intervalo de busca.

A Interpolação Quadrática pode ser eficiente quando os pontos iniciais representam adequadamente a região do mínimo e a parábola fornece uma aproximação útil.

Por isso, comparar o resultado, o custo mínimo, o número de iterações e as avaliações da função objetivo ajuda a compreender o comportamento dos métodos.

---

## 15. Exploração numérica da função objetivo

Antes de executar a otimização, o programa calcula os custos para algumas espessuras previamente selecionadas.

### Código

```python
espessuras_teste = [
    0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12
]

print(f"{'e (m)':<8} | {'CE (R$)':<12} | {'CI (R$)':<12} | {'C Total (R$)'}")

for e in espessuras_teste:
    CE, CI, C_tot = calcular_custos(e, pE_padrao)
    print(f"{e:<8.2f} | {CE:<12.2f} | {CI:<12.2f} | {C_tot:.2f}")
```

### Como funciona?

**1. Definição das espessuras**

A lista contém sete valores, desde 0 até 0,12 m.

**2. Impressão do cabeçalho**

O cabeçalho identifica as colunas da tabela.

**3. Laço `for`**

O programa percorre cada espessura.

**4. Cálculo dos custos**

A função `calcular_custos` retorna o custo de energia, o custo de instalação e o custo total.

**5. Formatação**

Os valores são impressos com duas casas decimais, facilitando a comparação.

### Por que essa etapa é importante?

Essa exploração permite observar o comportamento econômico do modelo antes de analisar o resultado dos métodos de otimização.

Em particular, ela ajuda a verificar se:

- O custo de energia varia com a espessura;
- O custo de instalação aumenta com a quantidade de material;
- O custo total apresenta uma região de mínimo no intervalo analisado.

A tabela é uma amostragem da função. Ela não garante, sozinha, que o menor valor esteja exatamente em um dos pontos testados.

---

## 16. Gráfico do custo total em função da espessura

O primeiro gráfico representa:

\[
C(e)=C_E(e)+C_I(e)
\]

em função da espessura \(e\).

### Código

```python
e_plot = np.linspace(0, 0.12, 100)
c_plot = [
    funcao_objetivo(e, pE_padrao)
    for e in e_plot
]

plt.figure(figsize=(8, 5))
plt.plot(e_plot, c_plot, 'k-', label='Custo Total C(e)')
plt.title('Custo Total vs Espessura de Isolamento (pE = 0.70 R$/kWh)')
plt.xlabel('Espessura e (m)')
plt.ylabel('Custo Total (R$)')
plt.grid(True)
plt.legend()
plt.savefig('grafico_custo_e.png')
```

### Explicação

**1. Geração dos pontos de espessura**

```python
e_plot = np.linspace(0, 0.12, 100)
```

Cria 100 valores igualmente espaçados no intervalo de 0 a 0,12 m.

**2. Avaliação da função objetivo**

```python
c_plot = [
    funcao_objetivo(e, pE_padrao)
    for e in e_plot
]
```

Calcula o custo total para cada espessura.

**3. Criação da figura**

```python
plt.figure(figsize=(8, 5))
```

Define o tamanho da figura.

**4. Desenho da curva**

```python
plt.plot(e_plot, c_plot, 'k-', label='Custo Total C(e)')
```

Desenha a curva do custo total em função da espessura.

**5. Título e rótulos**

Os comandos `plt.title`, `plt.xlabel` e `plt.ylabel` identificam o gráfico e suas unidades.

**6. Grade e legenda**

```python
plt.grid(True)
plt.legend()
```

Adicionam uma grade e uma legenda para facilitar a leitura.

**7. Salvamento**

```python
plt.savefig('grafico_custo_e.png')
```

Salva o gráfico no arquivo `grafico_custo_e.png`.

### Como interpretar o gráfico?

A curva representa o equilíbrio entre dois componentes:

- O custo de energia, associado à perda de calor;
- O custo de instalação, associado ao volume de isolamento.

A espessura ótima corresponde ao ponto em que o custo total é mínimo dentro do intervalo considerado.

O gráfico auxilia na interpretação visual, enquanto os métodos numéricos fornecem uma estimativa quantitativa da localização do mínimo.

---

## 17. Comparação numérica dos métodos

O programa executa os dois métodos com o preço de energia de referência.

### Código

```python
a_opt, a_custo, a_iter, a_evals, a_status = (
    otimizar_razao_aurea(
        0.0, 0.12, pE_padrao
    )
)

q_opt, q_custo, q_iter, q_evals, q_status = (
    otimizar_interp_quadratica(
        0.04, 0.07, 0.11, pE_padrao
    )
)
```

### Método da Razão Áurea

O intervalo de busca é:

\[
[0,\ 0{,}12]\ \text{m}
\]

O método procura o mínimo nesse intervalo.

### Método da Interpolação Quadrática

Os pontos iniciais são:

\[
0{,}04,\quad 0{,}07,\quad 0{,}11\ \text{m}
\]

A função utiliza esses pontos para construir a parábola interpoladora.

### Impressão dos resultados

```python
print(f"{'Resultado':<25} | {'Razão Áurea':<18} | {'Interpolação Quadrática'}")
print(f"{'-'*70}")
print(f"{'Espessura ótima (m)':<25} | {a_opt:<18.6f} | {q_opt:.6f}")
print(f"{'Custo mínimo (R$)':<25} | {a_custo:<18.2f} | {q_custo:.2f}")
print(f"{'Número de iterações':<25} | {a_iter:<18} | {q_iter}")
print(f"{'Avaliações de C(e)':<25} | {a_evals:<18} | {q_evals}")
print(f"{'Critério/Situação':<25} | {a_status:<18} | {q_status}")
```

A tabela permite comparar:

- A espessura ótima estimada;
- O custo mínimo correspondente;
- O número de iterações;
- As avaliações contabilizadas da função objetivo;
- A situação final de cada método.

**Atenção à interpretação:** o número de avaliações contabilizadas pela Razão Áurea subestima o total de chamadas porque não inclui a avaliação do custo no ponto médio retornado. Já a interpolação quadrática também pode encerrar antecipadamente por uma proteção, situação que deve ser distinguida da convergência efetiva.

---

## 18. Investigação paramétrica: influência do preço da energia

O programa investiga como a espessura ótima muda quando o preço da energia térmica varia.

### Código

```python
precos_E = [0.35, 0.70, 1.05, 1.40]
e_otimos_grafico = []

print(
    f"{'pE (R$/kWh)':<12} | {'e_ótimo (m)':<12} | "
    f"{'CE (R$)':<10} | {'CI (R$)':<10} | {'Custo Mínimo (R$)'}"
)

for pE in precos_E:
    e_otimo, custo_min, _, _, _ = otimizar_razao_aurea(
        0.0, 0.12, pE
    )

    CE, CI, _ = calcular_custos(e_otimo, pE)
    e_otimos_grafico.append(e_otimo)

    print(
        f"{pE:<12.2f} | {e_otimo:<12.4f} | "
        f"{CE:<10.2f} | {CI:<10.2f} | {custo_min:.2f}"
    )
```

### 18.1. Preços utilizados

A lista contém os seguintes preços:

\[
p_E\in\{0{,}35;\ 0{,}70;\ 1{,}05;\ 1{,}40\}
\]

Todos estão expressos em R$/kWh.

### 18.2. Procedimento

Para cada preço de energia, o programa:

1. Executa o método da Razão Áurea;
2. Obtém a espessura ótima estimada;
3. Calcula os custos de energia e instalação nessa espessura;
4. Armazena a espessura para o gráfico;
5. Imprime os resultados.

### 18.3. Interpretação econômica

Quando o preço da energia aumenta, a perda de calor passa a ter maior impacto econômico.

Por isso, em um modelo como este, pode ser vantajoso utilizar uma espessura maior para reduzir a perda de energia.

Entretanto, o custo de instalação também cresce com a espessura. O resultado ótimo depende do equilíbrio entre esses componentes e do intervalo de busca.

A relação exata deve ser analisada a partir dos resultados numéricos, em vez de presumir uma variação específica para todos os cenários.

---

## 19. Gráfico da influência do preço da energia

O segundo gráfico mostra a relação entre o preço da energia e a espessura ótima encontrada.

### Código

```python
plt.figure(figsize=(8, 5))
plt.plot(precos_E, e_otimos_grafico, 'ro-')
plt.title('Influência do Preço da Energia na Espessura Ótima')
plt.xlabel('Preço da Energia pE (R$/kWh)')
plt.ylabel('Espessura Ótima do Isolamento (m)')
plt.grid(True)
plt.savefig('grafico_e_pE.png')
```

### Interpretação

- O eixo horizontal representa o preço da energia;
- O eixo vertical representa a espessura ótima estimada;
- Cada ponto corresponde a uma simulação;
- A linha conecta os pontos para facilitar a visualização da tendência.

O gráfico é útil para analisar a sensibilidade econômica do isolamento.

Ele não constitui uma demonstração de que a relação seja linear: o desenho apenas conecta os resultados calculados para os preços selecionados.

O arquivo produzido pelo programa é:

`grafico_e_pE.png`

---

## 20. Verificação da consistência física: raio crítico de isolamento

O código termina calculando o raio crítico de isolamento cilíndrico.

### Código

```python
rc = k / h

print(f"Raio crítico calculado: {rc:.5f} m ({rc*1000:.2f} mm)")
print(f"Raio da tubulação (r1): {r1:.5f} m ({r1*1000:.2f} mm)")
```

### 20.1. Fórmula do raio crítico

Para um cilindro com condução radial e convecção externa, o raio crítico é:

\[
\boxed{r_c=\frac{k}{h}}
\]

Em que:

- \(k\) é a condutividade térmica do isolante;
- \(h\) é o coeficiente de convecção externa.

### 20.2. Substituição dos valores

Utilizando os parâmetros do programa:

\[
r_c=\frac{0{,}045}{12}
\]

\[
\boxed{r_c=0{,}00375\ \text{m}}
\]

Convertendo para milímetros:

\[
r_c=0{,}00375\times1000
\]

\[
\boxed{r_c=3{,}75\ \text{mm}}
\]

O raio original da tubulação é:

\[
r_1=0{,}040\ \text{m}=40\ \text{mm}
\]

Logo:

\[
r_1>r_c
\]

### 20.3. Significado físico

O raio crítico indica a região em que o aumento do raio externo pode inicialmente aumentar a perda de calor em sistemas cilíndricos, devido ao crescimento da área de convecção externa.

No modelo ideal de resistência térmica utilizado, como o raio original é maior que o raio crítico, aumentar o raio externo por meio do isolamento reduz a perda de calor.

Isso não significa que a espessura economicamente ótima seja infinita: o custo de instalação aumenta e limita o benefício econômico de adicionar mais material.

O raio crítico é uma verificação física complementar. Ele não determina sozinho a espessura ótima econômica.

---

## 21. Fluxo completo de execução

O funcionamento geral do programa pode ser resumido da seguinte maneira:

1. Importa as bibliotecas necessárias;
2. Define os parâmetros físicos e econômicos;
3. Calcula a perda de calor para uma espessura;
4. Calcula os custos de energia e instalação;
5. Define o custo total como função objetivo;
6. Avalia custos para espessuras selecionadas;
7. Gera o gráfico do custo total;
8. Executa a otimização pela Razão Áurea;
9. Executa a otimização pela Interpolação Quadrática;
10. Imprime a comparação dos resultados;
11. Repete a otimização para diferentes preços de energia;
12. Gera o gráfico da espessura ótima em função do preço;
13. Calcula e imprime o raio crítico.

---

## 22. Complexidade computacional

O custo computacional dos métodos depende principalmente do número de avaliações da função objetivo.

A função objetivo envolve um número fixo de operações matemáticas, como logaritmos, multiplicações e divisões. Assim, considerando o custo dessas operações como constante, cada avaliação possui complexidade \(O(1)\).

### Razão Áurea

Se forem realizadas \(N\) iterações, o custo é aproximadamente:

\[
O(N)
\]

O método reduz o intervalo em uma proporção constante a cada iteração, quando aplicado a uma função unimodal no intervalo de busca.

### Interpolação Quadrática

Se forem realizadas \(N\) iterações, o custo também é aproximadamente:

\[
O(N)
\]

Cada iteração realiza um número constante de operações e, normalmente, uma nova avaliação da função objetivo.

Entretanto, o método pode encerrar antecipadamente por suas proteções numéricas.

### Geração dos gráficos

O primeiro gráfico avalia a função em 100 pontos. Portanto, seu custo é aproximadamente proporcional a 100 avaliações da função objetivo.

A investigação paramétrica executa quatro otimizações pela Razão Áurea, uma para cada preço de energia.

---

## 23. Limitações e cuidados importantes

### 23.1. Intervalo de busca limitado

Os métodos não procuram espessuras fora dos intervalos fornecidos.

Na Razão Áurea, o intervalo é de 0 a 0,12 m. Portanto, o resultado é o mínimo estimado dentro desse intervalo, e não necessariamente o mínimo global sem restrições.

### 23.2. Hipótese de unimodalidade

A Razão Áurea pressupõe que a função possua um único mínimo relevante no intervalo.

Se a função apresentar vários mínimos locais, a estratégia pode não encontrar o mínimo global.

### 23.3. Pontos iniciais da interpolação quadrática

Os pontos iniciais devem ser distintos e ordenados para que as expressões utilizadas sejam válidas e as regras de atualização funcionem conforme esperado.

O código não valida explicitamente todos esses requisitos.

### 23.4. Proteções numéricas

A interpolação quadrática possui verificações para curvatura pequena, aproximações fora do intervalo e pontos praticamente coincidentes.

Essas verificações evitam algumas situações problemáticas, mas podem fazer o método terminar antes de atingir a tolerância.

### 23.5. Contabilização das avaliações

O contador de avaliações da Razão Áurea não inclui a chamada final da função objetivo para calcular o custo no ponto médio retornado.

Portanto, os contadores não devem ser interpretados sem considerar esse detalhe da implementação.

### 23.6. Modelo econômico simplificado

O custo total é a soma do custo instalado do isolamento e do custo acumulado da energia perdida.

O modelo não inclui explicitamente inflação, taxa de desconto, manutenção, substituição do isolante, variações temporais do preço de energia ou outros custos operacionais.

Além disso, o tempo de operação é considerado fixo e o preço de energia é considerado constante durante esse período.

### 23.7. Hipóteses do modelo térmico

O cálculo utiliza propriedades constantes e resistências térmicas por unidade de comprimento para condução radial e convecção externa.

O modelo não inclui explicitamente resistência térmica da parede da tubulação, resistência de contato, radiação térmica ou variações da temperatura ao longo do comprimento.

Esses efeitos podem ser importantes em uma aplicação real, dependendo das condições de operação.

### 23.8. Ausência de validação dos parâmetros

O código não valida explicitamente valores como espessura negativa, raio não positivo, condutividade não positiva, coeficiente de convecção não positivo ou parâmetros econômicos incompatíveis.

Para aplicações gerais, é recomendável incluir essas verificações.

---

## 24. Variáveis e funções para consulta rápida

### Variáveis principais

| Nome | Descrição |
|---|---|
| `r1` | Raio original da tubulação |
| `L` | Comprimento da tubulação |
| `Ti` | Temperatura na interface |
| `T_inf` | Temperatura ambiente |
| `k` | Condutividade térmica do isolante |
| `h` | Coeficiente de convecção |
| `t` | Tempo acumulado de operação |
| `pE_padrao` | Preço de referência da energia |
| `cI` | Custo instalado do isolante por volume |
| `e` | Espessura do isolamento |
| `r2` | Raio externo após a aplicação do isolante |
| `CE` | Custo da energia térmica perdida |
| `CI` | Custo de instalação |
| `C_tot` | Custo total |
| `rc` | Raio crítico de isolamento |

### Funções principais

| Função | Finalidade |
|---|---|
| `calcular_q_linha(e)` | Calcula a perda de calor por unidade de comprimento |
| `calcular_custos(e, pE)` | Calcula os custos de energia, instalação e total |
| `funcao_objetivo(e, pE)` | Retorna o custo total para uma espessura |
| `otimizar_razao_aurea(a, b, pE, ...)` | Otimiza o custo total por redução do intervalo |
| `otimizar_interp_quadratica(x1, x2, x3, pE, ...)` | Otimiza o custo total por interpolação quadrática |

---

## 25. Como executar o projeto

### Requisitos

- Python 3;
- NumPy;
- Matplotlib.

Instale as dependências com:

```bash
pip install numpy matplotlib
```

Salve o código em um arquivo Python, por exemplo:

```text
Otimização de isolamento térmico.py
```

Execute no terminal:

```bash
python "Otimização de isolamento térmico.py"
```

O programa imprimirá:

- A tabela de custos para diferentes espessuras;
- A comparação entre os métodos de otimização;
- Os resultados da investigação paramétrica;
- O raio crítico calculado.

Também serão salvos dois gráficos no diretório de execução:

```text
grafico_custo_e.png
grafico_e_pE.png
```

O código não chama explicitamente `plt.show()`. Por isso, a forma principal de visualização prevista pela implementação é a abertura dos arquivos PNG salvos.

---

## 26. Conclusão

O projeto combina transferência de calor, avaliação econômica e otimização numérica para estudar a espessura de isolamento térmico de uma tubulação.

A função térmica calcula a perda de calor por unidade de comprimento usando resistências de condução e convecção. A função econômica combina o custo da energia térmica perdida durante o período de operação com o custo de instalação do isolamento.

A função objetivo resultante é minimizada por dois métodos: Razão Áurea e Interpolação Quadrática. A comparação permite observar como diferentes estratégias numéricas encontram uma espessura de baixo custo.

A análise paramétrica amplia o estudo ao investigar a influência do preço da energia sobre a espessura ótima, enquanto o cálculo do raio crítico oferece uma verificação física complementar.

O resultado deve ser interpretado dentro das hipóteses, dos parâmetros e dos limites de busca do modelo. Para aplicações de engenharia, é importante validar as condições físicas, os custos considerados e a convergência efetiva dos métodos numéricos antes de utilizar a espessura calculada como decisão de projeto.
