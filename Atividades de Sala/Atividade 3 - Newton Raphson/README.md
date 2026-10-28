# Método de Newton–Raphson para Sistemas Não Lineares: Bifurcação de Tubulação

Atividade prática de **Cálculo Numérico Aplicado**. O método de Newton–Raphson para sistemas não lineares é implementado em Python (apenas **NumPy**) para determinar como uma vazão total se divide entre dois ramos paralelos de uma tubulação, considerando perdas distribuídas (Darcy–Weisbach) e perdas localizadas (acessórios e válvula).

## Sumário

- [Problema físico](#problema-físico)
- [Formulação matemática](#formulação-matemática)
- [Método de Newton–Raphson](#método-de-newtonraphson)
- [Implementação](#implementação)
- [Como executar](#como-executar)
- [Resultados](#resultados)
- [Verificação física](#verificação-física)
- [Discussão](#discussão)
- [Estrutura do repositório](#estrutura-do-repositório)

## Problema físico

Uma vazão volumétrica total conhecida $Q$ chega a um ponto **A**, onde o escoamento se divide em dois ramos paralelos que voltam a se encontrar em um ponto **B**. O ramo 2 possui uma válvula parcialmente aberta. As vazões $Q_1$ e $Q_2$ em cada ramo são inicialmente desconhecidas.

| Parâmetro | Ramo 1 | Ramo 2 |
|---|---|---|
| Comprimento $L$ (m) | 30 | 20 |
| Diâmetro $D$ (m) | 0,10 | 0,08 |
| Fator de atrito $f$ | 0,020 | 0,022 |
| Perda localizada $K$ | 1,5 | 1,0 |
| Válvula $K_v$ | n/a | 8,0 |

Dados gerais: $Q = 0{,}050\ \text{m}^3/\text{s}$ e $g = 9{,}81\ \text{m/s}^2$. Os fatores de atrito são considerados conhecidos e constantes.

As perdas de carga são dadas por:

$$h_f = f\,\frac{L}{D}\,\frac{U^2}{2g}, \qquad h_m = K\,\frac{U^2}{2g}, \qquad U = \frac{4Q}{\pi D^2}$$

## Formulação matemática

Substituindo $U_i = 4Q_i/(\pi D_i^2)$ nas perdas, cada ramo fica com $h_{L,i} = R_i Q_i^2$, onde:

$$R_1 = \left(f_1\frac{L_1}{D_1} + K_1\right)\frac{8}{g\pi^2 D_1^4} \approx 6197{,}01\ \text{s}^2/\text{m}^5$$

$$R_2 = \left(f_2\frac{L_2}{D_2} + K_2 + K_v\right)\frac{8}{g\pi^2 D_2^4} \approx 29\,250{,}23\ \text{s}^2/\text{m}^5$$

Duas condições físicas precisam ser satisfeitas ao mesmo tempo:

1. **Conservação da vazão** na bifurcação: $Q = Q_1 + Q_2$.
2. **Igualdade das perdas de carga** entre A e B pelos dois caminhos: $h_{L,1} = h_{L,2}$.

Isso leva ao sistema não linear:

$$
\begin{cases}
F_1(Q_1,Q_2) = Q_1 + Q_2 - Q = 0 \\[4pt]
F_2(Q_1,Q_2) = R_1 Q_1^2 - R_2 Q_2^2 = 0
\end{cases}
$$

## Método de Newton–Raphson

Em cada iteração $k$, o sistema é linearizado em torno da aproximação atual $(Q_1^{(k)}, Q_2^{(k)})$ e se busca a correção $(\Delta Q_1, \Delta Q_2)$ tal que $Q_i^{(k+1)} = Q_i^{(k)} + \Delta Q_i$. A correção vem de:

$$
J^{(k)}
\begin{bmatrix}\Delta Q_1 \\ \Delta Q_2\end{bmatrix}
= -
\begin{bmatrix}F_1^{(k)} \\ F_2^{(k)}\end{bmatrix}
$$

A matriz Jacobiana, obtida derivando as duas equações à mão, é:

$$
J =
\begin{bmatrix}
\dfrac{\partial F_1}{\partial Q_1} & \dfrac{\partial F_1}{\partial Q_2} \\[8pt]
\dfrac{\partial F_2}{\partial Q_1} & \dfrac{\partial F_2}{\partial Q_2}
\end{bmatrix}
=
\begin{bmatrix}
1 & 1 \\
2R_1Q_1 & -2R_2Q_2
\end{bmatrix}
$$

Como o sistema é $2\times2$, as correções são obtidas pela **regra de Cramer**. Para $a\,\Delta Q_1 + b\,\Delta Q_2 = r_1$ e $c\,\Delta Q_1 + d\,\Delta Q_2 = r_2$:

$$D = ad - bc, \qquad \Delta Q_1 = \frac{r_1 d - b r_2}{D}, \qquad \Delta Q_2 = \frac{a r_2 - r_1 c}{D}$$

Aqui, $a = 1$, $b = 1$, $c = 2R_1Q_1$, $d = -2R_2Q_2$, $r_1 = -F_1$ e $r_2 = -F_2$. Note que $D = -2R_2Q_2 - 2R_1Q_1 < 0$ para vazões positivas, então o sistema nunca é singular nesse regime.

**Parâmetros do método:**

- **Estimativa inicial:** $Q_1^{(0)} = 0{,}030$ e $Q_2^{(0)} = 0{,}020\ \text{m}^3/\text{s}$. Ela já satisfaz a conservação da vazão, mas não a igualdade das perdas.
- **Critério de parada:** $\max(|\Delta Q_1|, |\Delta Q_2|) < 10^{-8}\ \text{m}^3/\text{s}$.
- **Limite de iterações:** 100, para evitar laço infinito caso o método não convirja.

## Implementação

O código está em [`newton_raphson_nao_linear.py`](newton_raphson_nao_linear.py) e usa apenas NumPy. Ele está organizado em quatro blocos:

1. **Parâmetros do problema:** vazão total, gravidade, tolerância, número máximo de iterações e os dados de cada ramo.
2. **Cálculo de $R_1$ e $R_2$:** as constantes de resistência, com `np.pi` para o valor de $\pi$.
3. **Funções do sistema:**
   - `funcoes_F1_F2(Q1, Q2)` devolve os resíduos $F_1$ e $F_2$.
   - `matriz_jacobiana(Q1, Q2)` devolve os quatro elementos da Jacobiana.
4. **Laço de Newton–Raphson:** a cada iteração, avalia os resíduos e a Jacobiana, calcula o determinante $D$, resolve por Cramer, imprime o histórico, atualiza as vazões e testa a convergência.

Trecho central do laço:

```python
D = J11 * J22 - J12 * J21
r1 = -F1
r2 = -F2

# Regra de Cramer para as correções dQ1 e dQ2:
dQ1 = ((r1 * J22) - (J12 * r2)) / D
dQ2 = ((J11 * r2) - (r1 * J21)) / D

Q1_novo = Q1 + dQ1
Q2_novo = Q2 + dQ2

erro = max(abs(dQ1), abs(dQ2))
if erro < tol:
    break
```

## Como executar

Requisitos: Python 3.8+ e **NumPy**.

```bash
pip install numpy
python newton_raphson_nao_linear.py
```

## Resultados

Saída do programa (iteração, vazões atuais, resíduos e correções):

```
Iteração 1, Q1 = 0.03, Q2 = 0.02, F1 = 0.0, F2 = -6.12278, dQ1 = 0.00397111, dQ2 = -0.00397111
Iteração 2, Q1 = 0.0339711, Q2 = 0.0160289, F1 = 0.0, F2 = -0.363543, dQ1 = 0.00026756, dQ2 = -0.00026756
Iteração 3, Q1 = 0.0342387, Q2 = 0.0157613, F1 = 0.0, F2 = -0.00165034, dQ1 = 1.22574e-06, dQ2 = -1.22574e-06
Iteração 4, Q1 = 0.0342399, Q2 = 0.0157601, F1 = 6.93889e-18, F2 = -3.46361e-08, dQ1 = 2.5726e-11, dQ2 = -2.5726e-11
Convergiu em 4 iterações
Vazões finais:Q1 = 0.034240 m**3/s, Q2 = 0.015760 m**3/s
```

Em resumo:

| k | $Q_1$ (m³/s) | $Q_2$ (m³/s) | $F_1$ | $F_2$ | $\Delta Q_1$ | $\Delta Q_2$ |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | 0,030000 | 0,020000 | 0 | −6,12278 | +3,97×10⁻³ | −3,97×10⁻³ |
| 2 | 0,033971 | 0,016029 | 0 | −0,363543 | +2,68×10⁻⁴ | −2,68×10⁻⁴ |
| 3 | 0,034239 | 0,015761 | 0 | −1,65×10⁻³ | +1,23×10⁻⁶ | −1,23×10⁻⁶ |
| 4 | 0,034240 | 0,015760 | ~10⁻¹⁷ | −3,46×10⁻⁸ | +2,57×10⁻¹¹ | −2,57×10⁻¹¹ |

**Solução convergida em 4 iterações:**

| Grandeza | Ramo 1 | Ramo 2 |
|---|---|---|
| Vazão $Q_i$ (m³/s) | 0,034240 | 0,015760 |
| Parcela da vazão total | 68,48 % | 31,52 % |

Duas observações sobre o histórico:

- Como $F_1 = 0$ na estimativa inicial, vale $\Delta Q_1 = -\Delta Q_2$ em todas as iterações. A correção preserva a conservação da vazão e apenas redistribui o escoamento entre os ramos.
- A convergência é **quadrática**: a correção cai de $10^{-3}$ para $10^{-4}$, $10^{-6}$ e $10^{-11}$.

## Verificação física

As velocidades médias, a conservação da vazão e a igualdade das perdas podem ser conferidas acrescentando o trecho abaixo ao final do script, após o laço:

```python
# Verificação física
U1 = 4 * Q1_novo / (np.pi * D1**2)
U2 = 4 * Q2_novo / (np.pi * D2**2)
print(f"U1 = {U1:.3f} m/s, U2 = {U2:.3f} m/s")
print(f"Q1 + Q2 = {Q1_novo + Q2_novo:.10f} m**3/s (Q = {Q})")
print(f"hL1 = {R1 * Q1_novo**2:.6f} m, hL2 = {R2 * Q2_novo**2:.6f} m")
```

Saída:

```
U1 = 4.360 m/s, U2 = 3.135 m/s
Q1 + Q2 = 0.0500000000 m**3/s (Q = 0.05)
hL1 = 7.265197 m, hL2 = 7.265197 m
```

| Verificação | Resultado |
|---|---|
| Velocidade média $U_1$ | 4,360 m/s |
| Velocidade média $U_2$ | 3,135 m/s |
| $Q_1 + Q_2$ | 0,0500000000 m³/s (igual a $Q$) |
| $h_{L,1} = R_1 Q_1^2$ | 7,265197 m |
| $h_{L,2} = R_2 Q_2^2$ | 7,265197 m |

A solução satisfaz simultaneamente a conservação da vazão e a igualdade das perdas de carga. Há também uma checagem analítica: de $R_1Q_1^2 = R_2Q_2^2$ vem $Q_1/Q_2 = \sqrt{R_2/R_1} = 2{,}1726$, e com $Q_1 + Q_2 = Q$ obtém-se $Q_2 = Q/(1 + 2{,}1726) = 0{,}01576\ \text{m}^3/\text{s}$, o mesmo valor do Newton.

## Discussão

**Qual ramo transporta mais vazão?**
O ramo 1, com cerca de 68,5 % do total. Como $R_1 \ll R_2$ ($6197$ contra $29\,250\ \text{s}^2/\text{m}^5$), o escoamento prefere o caminho de menor resistência, e vale $Q_1/Q_2 = \sqrt{R_2/R_1}$. O ramo 2 tem diâmetro menor (e $R \propto D^{-4}$), mais acessórios e a válvula.

**Quantas iterações foram necessárias?**
Quatro, com o critério $10^{-8}\ \text{m}^3/\text{s}$.

**Por que o sistema vira linear a cada iteração?**
O método substitui $\mathbf{F}$ por sua expansão de Taylor de primeira ordem, $\mathbf{F}(\mathbf{Q}+\Delta\mathbf{Q}) \approx \mathbf{F}(\mathbf{Q}) + J\,\Delta\mathbf{Q}$, desprezando os termos de ordem superior em $\Delta\mathbf{Q}$. Os termos não lineares ($Q^2$) só entram no cálculo de $\mathbf{F}$ e $J$, que são apenas números em cada iteração. Assim, as incógnitas $\Delta Q_1$ e $\Delta Q_2$ aparecem apenas de forma linear.

**Qual o papel da Jacobiana?**
Ela generaliza a derivada $f'(x)$ do caso unidimensional. Cada elemento mede a sensibilidade de uma equação a uma incógnita, e juntos definem o plano tangente que aproxima o sistema. Resolver o sistema linear com $J$ dá a correção que anula os resíduos nessa aproximação.

**E se mudar a estimativa inicial?**
Basta alterar `Q1` e `Q2` antes do laço. Em testes com outras estimativas (usando o mesmo algoritmo), o método convergiu sempre para a mesma solução quando as vazões iniciais eram positivas, mudando apenas o número de iterações:

| Estimativa $(Q_1^{(0)}, Q_2^{(0)})$ | Solução $(Q_1, Q_2)$ | Iterações |
|:-:|:-:|:-:|
| (0,030; 0,020) | (0,034240; 0,015760) | 4 |
| (0,025; 0,025) | (0,034240; 0,015760) | 4 |
| (0,050; 0,000) | (0,034240; 0,015760) | 5 |
| (0,010; 0,040) | (0,034240; 0,015760) | 5 |
| (0,200; 0,200) | (0,034240; 0,015760) | 7 |
| (0,100; −0,050) | (0,092641; **−0,042641**) | 4 |

A última linha mostra uma **segunda raiz matemática**, com escoamento invertido no ramo 2, sem sentido físico para o problema da figura. Isso ocorre porque $F_2$ usa $Q^2$ e perde o sinal da perda de carga (que seria $RQ|Q|$). Por isso, a escolha de uma estimativa inicial positiva importa.

## Estrutura do repositório

```
.
├── README.md                        # este arquivo
└── newton_raphson_nao_linear.py     # implementação do método
```
