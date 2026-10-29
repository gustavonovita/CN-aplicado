# Matplotlib.pyplot — Guia de consulta em Python

Este arquivo reúne as principais funções do módulo **`matplotlib.pyplot`**, com explicações e exemplos práticos para criação e personalização de gráficos.

O Matplotlib é uma biblioteca de visualização para Python. O módulo `pyplot` fornece uma interface simples para criar gráficos de linha, dispersão, barras, histogramas, pizza e diversos outros tipos de visualização.

## Instalação

```bash
python -m pip install numpy matplotlib
```

Se estiver usando um ambiente virtual:

```bash
python -m venv .venv
```

No Windows:

```bash
.venv\Scripts\activate
```

No Linux ou macOS:

```bash
source .venv/bin/activate
```


## Importação

```python
import matplotlib.pyplot as plt
```

O apelido `plt` é a convenção mais utilizada para o módulo `matplotlib.pyplot`.

# Matplotlib.pyplot

A convenção de importação é:

```python
import matplotlib.pyplot as plt
```

O fluxo mais comum é:

1. Preparar os dados.
2. Criar a figura e os eixos.
3. Desenhar os dados.
4. Adicionar título, legendas e rótulos.
5. Exibir ou salvar o gráfico.

Exemplo básico:

```python
import matplotlib.pyplot as plt

x = [1, 2, 3, 4]
y = [1, 4, 9, 16]

plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Exemplo básico")
plt.show()
```

> Em scripts Python, `plt.show()` normalmente abre a janela do gráfico. Em notebooks, o gráfico costuma aparecer automaticamente, mas ainda é comum usar `plt.show()` para deixar o código explícito.

---

## Gráficos de linha

A função principal é `plt.plot()`.

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 100)
y = x ** 2

plt.plot(x, y)
plt.show()
```

### Cor, estilo e marcador

```python
plt.plot(
    x,
    y,
    color="blue",
    linestyle="--",
    linewidth=2,
    marker="o",
    markersize=4,
    label="y = x²"
)
```

Atalhos de estilo comuns:

```python
plt.plot(x, y, "r--")  # vermelho e tracejado
plt.plot(x, y, "bo-")  # azul, círculos e linha contínua
```

Cores frequentes:

- `"r"`: vermelho.
- `"g"`: verde.
- `"b"`: azul.
- `"k"`: preto.
- `"y"`: amarelo.
- `"m"`: magenta.
- `"c"`: ciano.

Estilos de linha:

- `"-"`: contínua.
- `"--"`: tracejada.
- `"-."`: traço e ponto.
- `":"`: pontilhada.

### Mais de uma linha

```python
x = np.linspace(0, 2 * np.pi, 200)

plt.plot(x, np.sin(x), label="seno")
plt.plot(x, np.cos(x), label="cosseno")

plt.title("Funções trigonométricas")
plt.xlabel("Ângulo em radianos")
plt.ylabel("Valor")
plt.legend()
plt.grid(True)
plt.show()
```

---

## Gráficos de dispersão

Use `plt.scatter()` para visualizar a relação entre duas variáveis.

```python
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)
x = rng.random(50)
y = 2 * x + rng.normal(0, 0.1, 50)

plt.scatter(x, y)
plt.xlabel("Variável X")
plt.ylabel("Variável Y")
plt.title("Gráfico de dispersão")
plt.show()
```

### Tamanho e cor dos pontos

```python
tamanhos = rng.integers(20, 200, size=50)
cores = y

plt.scatter(
    x,
    y,
    s=tamanhos,
    c=cores,
    cmap="viridis",
    alpha=0.7,
    edgecolors="black"
)
plt.colorbar(label="Valor de Y")
plt.show()
```

Parâmetros úteis:

- `s`: tamanho dos pontos.
- `c`: cor dos pontos.
- `cmap`: mapa de cores.
- `alpha`: transparência, entre `0` e `1`.
- `edgecolors`: cor da borda.

---

## Gráficos de barras

Use `plt.bar()` para comparar categorias.

```python
categorias = ["Python", "Java", "C++", "JavaScript"]
quantidades = [80, 55, 40, 70]

plt.bar(categorias, quantidades, color="steelblue")
plt.xlabel("Linguagem")
plt.ylabel("Quantidade")
plt.title("Popularidade por linguagem")
plt.show()
```

### Barras horizontais

```python
plt.barh(categorias, quantidades, color="seagreen")
plt.xlabel("Quantidade")
plt.ylabel("Linguagem")
plt.show()
```

### Exibir valores sobre as barras

```python
fig, ax = plt.subplots()
barras = ax.bar(categorias, quantidades, color="orange")

ax.bar_label(barras, padding=3)
ax.set_ylim(0, max(quantidades) + 15)
ax.set_title("Valores nas barras")

plt.show()
```

### Barras agrupadas

```python
categorias = ["A", "B", "C"]
valores_2024 = [10, 15, 12]
valores_2025 = [14, 13, 18]

posicoes = np.arange(len(categorias))
largura = 0.35

fig, ax = plt.subplots()
ax.bar(posicoes - largura / 2, valores_2024, largura, label="2024")
ax.bar(posicoes + largura / 2, valores_2025, largura, label="2025")

ax.set_xticks(posicoes)
ax.set_xticklabels(categorias)
ax.set_ylabel("Valor")
ax.set_title("Comparação entre anos")
ax.legend()

plt.show()
```

---

## Histogramas

Use `plt.hist()` para observar a distribuição de uma variável numérica.

```python
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)
dados = rng.normal(loc=70, scale=10, size=1000)

plt.hist(dados, bins=20, color="purple", alpha=0.75, edgecolor="black")
plt.xlabel("Valor")
plt.ylabel("Frequência")
plt.title("Distribuição dos dados")
plt.show()
```

Parâmetros importantes:

- `bins`: quantidade de intervalos ou os limites dos intervalos.
- `density=True`: normaliza o histograma para representar densidade.
- `alpha`: transparência.
- `edgecolor`: cor das bordas.

Comparando duas distribuições:

```python
amostra_a = rng.normal(50, 8, 500)
amostra_b = rng.normal(70, 8, 500)

plt.hist(amostra_a, bins=20, alpha=0.6, label="Grupo A")
plt.hist(amostra_b, bins=20, alpha=0.6, label="Grupo B")
plt.legend()
plt.show()
```

---

## Gráficos de pizza

Use `plt.pie()` para mostrar proporções de um total.

```python
categorias = ["Python", "Java", "C++", "JavaScript"]
valores = [40, 25, 15, 20]

plt.pie(
    valores,
    labels=categorias,
    autopct="%.1f%%",
    startangle=90
)
plt.title("Distribuição percentual")
plt.show()
```

Parâmetros úteis:

- `labels`: nomes das categorias.
- `autopct`: formato dos percentuais.
- `startangle`: ângulo inicial.
- `explode`: separa uma ou mais fatias.
- `shadow`: adiciona sombra.

Exemplo com destaque:

```python
explode = [0.1, 0, 0, 0]

plt.pie(
    valores,
    labels=categorias,
    autopct="%.1f%%",
    explode=explode,
    shadow=True
)
plt.show()
```

---

## Subplots

`plt.subplots()` permite criar vários gráficos em uma mesma figura.

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 100)

fig, eixos = plt.subplots(2, 2, figsize=(10, 7))

# Gráfico superior esquerdo
eixos[0, 0].plot(x, x)
eixos[0, 0].set_title("y = x")

# Gráfico superior direito
eixos[0, 1].plot(x, x ** 2, color="red")
eixos[0, 1].set_title("y = x²")

# Gráfico inferior esquerdo
eixos[1, 0].plot(x, np.sqrt(x), color="green")
eixos[1, 0].set_title("y = √x")

# Gráfico inferior direito
eixos[1, 1].plot(x, np.log1p(x), color="purple")
eixos[1, 1].set_title("y = log(1 + x)")

fig.suptitle("Várias funções")
fig.tight_layout()
plt.show()
```

### `plt.subplots()` versus `plt.subplot()`

A abordagem recomendada é `plt.subplots()`, pois ela retorna explicitamente a figura e os eixos:

```python
fig, ax = plt.subplots()
```

Quando há apenas um eixo, `ax` representa esse eixo. Quando há vários, `ax` geralmente é um array de eixos.

---

## Personalização dos gráficos

### Título, eixos e legenda

```python
fig, ax = plt.subplots()

ax.plot([1, 2, 3], [2, 4, 3], label="Dados")
ax.set_title("Título do gráfico")
ax.set_xlabel("Eixo X")
ax.set_ylabel("Eixo Y")
ax.legend()

plt.show()
```

### Grade

```python
plt.grid(True)
```

Com maior controle:

```python
ax.grid(
    True,
    linestyle="--",
    alpha=0.5,
    color="gray"
)
```

### Limites dos eixos

```python
plt.xlim(0, 10)
plt.ylim(-1, 1)
```

Ou utilizando o objeto `Axes`:

```python
ax.set_xlim(0, 10)
ax.set_ylim(-1, 1)
```

### Marcas dos eixos

```python
ax.set_xticks([0, 2, 4, 6, 8, 10])
ax.set_yticks([-1, -0.5, 0, 0.5, 1])
```

### Escala logarítmica

```python
ax.set_xscale("log")
ax.set_yscale("log")
```

Também existem atalhos:

```python
plt.xscale("log")
plt.yscale("log")
```

### Anotações

```python
x = np.array([1, 2, 3, 4])
y = np.array([2, 5, 3, 7])

fig, ax = plt.subplots()
ax.plot(x, y, marker="o")
ax.annotate(
    "Maior valor",
    xy=(4, 7),
    xytext=(3, 5),
    arrowprops={"arrowstyle": "->"}
)
plt.show()
```

### Linhas de referência

```python
fig, ax = plt.subplots()
ax.plot(x, y)
ax.axhline(y=4, color="red", linestyle="--", label="Média")
ax.axvline(x=2, color="green", linestyle=":", label="Limite")
ax.legend()
plt.show()
```

### Formatar números no eixo

```python
from matplotlib.ticker import FuncFormatter

fig, ax = plt.subplots()
ax.plot([1, 2, 3], [1000, 2500, 4000])
ax.yaxis.set_major_formatter(FuncFormatter(lambda valor, pos: f"R$ {valor:,.0f}"))
plt.show()
```

---

## Tamanho, proporção e layout

### Tamanho da figura

```python
plt.figure(figsize=(10, 6))
plt.plot([1, 2, 3], [1, 4, 9])
plt.show()
```

O `figsize` é informado em polegadas: `(largura, altura)`.

### Ajustar o layout

```python
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_title("Título")
ax.set_xlabel("Eixo X")
ax.set_ylabel("Eixo Y")
fig.tight_layout()
plt.show()
```

Outra opção é criar a figura com:

```python
fig, ax = plt.subplots(constrained_layout=True)
```

---

## Salvar gráficos

Use `plt.savefig()` ou `fig.savefig()` antes de `plt.show()`.

```python
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot([1, 2, 3], [1, 4, 9])
ax.set_title("Gráfico salvo")
fig.tight_layout()
fig.savefig("grafico.png", dpi=300, bbox_inches="tight")
plt.show()
```

Formatos comuns:

```python
fig.savefig("grafico.png", dpi=300)
fig.savefig("grafico.pdf")
fig.savefig("grafico.svg")
```

Parâmetros úteis:

- `dpi=300`: resolução adequada para impressão.
- `bbox_inches="tight"`: reduz espaços vazios nas bordas.
- `transparent=True`: fundo transparente.

---

## Cores e mapas de cores

### Cores nomeadas

```python
plt.plot(x, y, color="darkorange")
```

### Código hexadecimal

```python
plt.plot(x, y, color="#2E86AB")
```

### Mapa de cores

Mapas de cores são úteis quando a cor representa uma variável numérica.

```python
plt.scatter(
    x,
    y,
    c=y,
    cmap="plasma"
)
plt.colorbar(label="Valor")
plt.show()
```

Mapas conhecidos:

- `viridis`
- `plasma`
- `inferno`
- `magma`
- `cividis`
- `coolwarm`
- `gray`

---

## Exemplo completo

O exemplo abaixo cria dados de vendas mensais, calcula estatísticas com NumPy e produz três gráficos.

```python
import numpy as np
import matplotlib.pyplot as plt

# Gerador reproduzível
rng = np.random.default_rng(42)

# Dados
meses = np.arange(1, 13)
vendas = np.array([120, 135, 128, 150, 160, 155, 170, 180, 175, 190, 205, 220])

# Estatísticas
media = np.mean(vendas)
mediana = np.median(vendas)
maior_venda = np.max(vendas)
mes_maior_venda = meses[np.argmax(vendas)]

print(f"Média: {media:.2f}")
print(f"Mediana: {mediana:.2f}")
print(f"Maior venda: {maior_venda}")
print(f"Mês com maior venda: {mes_maior_venda}")

# Criação da figura com três gráficos
fig, eixos = plt.subplots(1, 3, figsize=(16, 4))

# 1. Linha das vendas
eixos[0].plot(meses, vendas, marker="o", color="royalblue", label="Vendas")
eixos[0].axhline(media, color="red", linestyle="--", label=f"Média: {media:.1f}")
eixos[0].set_title("Vendas por mês")
eixos[0].set_xlabel("Mês")
eixos[0].set_ylabel("Quantidade")
eixos[0].legend()
eixos[0].grid(True, alpha=0.3)

# 2. Barras
eixos[1].bar(meses, vendas, color="seagreen")
eixos[1].set_title("Barras de vendas")
eixos[1].set_xlabel("Mês")
eixos[1].set_ylabel("Quantidade")

# 3. Histograma
eixos[2].hist(vendas, bins=5, color="darkorange", edgecolor="black")
eixos[2].set_title("Distribuição das vendas")
eixos[2].set_xlabel("Quantidade")
eixos[2].set_ylabel("Frequência")

fig.suptitle("Análise mensal de vendas", fontsize=16)
fig.tight_layout()
fig.savefig("analise_vendas.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## Exemplo: regressão linear simples

O NumPy pode calcular uma reta aproximada com `np.polyfit()`.

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.array([1, 2, 3, 4, 5, 6])
y = np.array([2, 4, 5, 8, 9, 12])

# Grau 1 significa uma reta: y = a*x + b
a, b = np.polyfit(x, y, deg=1)
y_estimado = a * x + b

print(f"Coeficiente angular: {a:.2f}")
print(f"Intercepto: {b:.2f}")

plt.scatter(x, y, label="Dados observados")
plt.plot(x, y_estimado, color="red", label="Reta ajustada")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Regressão linear simples")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

Para uma regressão linear mais completa, com métricas como `R²`, normalmente são utilizadas bibliotecas como `scikit-learn` ou `scipy`.

---

## Exemplo: imagem como matriz NumPy

Uma imagem pode ser representada como um array NumPy:

- Imagem em tons de cinza: array 2D `(altura, largura)`.
- Imagem colorida RGB: array 3D `(altura, largura, 3)`.

```python
import numpy as np
import matplotlib.pyplot as plt

imagem = np.random.default_rng(42).random((100, 100))

plt.imshow(imagem, cmap="gray")
plt.colorbar(label="Intensidade")
plt.title("Imagem representada por um array")
plt.axis("off")
plt.show()
```

Criando uma imagem RGB simples:

```python
imagem_rgb = np.zeros((100, 100, 3))
imagem_rgb[:, :, 0] = 1  # Canal vermelho

plt.imshow(imagem_rgb)
plt.axis("off")
plt.show()
```

---

## Erros comuns relacionados ao Matplotlib

### Esquecer os parênteses em `plt.show()`

Errado:

```python
plt.show
```

Correto:

```python
plt.show()
```

### Usar arrays com tamanhos diferentes

Os valores de `x` e `y` usados em `plot()` geralmente precisam ter o mesmo tamanho:

```python
x = np.array([1, 2, 3])
y = np.array([10, 20])

# Erro: tamanhos incompatíveis
# plt.plot(x, y)
```

> Se o exemplo usar `np`, importe também `numpy` com `import numpy as np`.

### Sobrepor gráficos acidentalmente

Em códigos longos, pode ser útil criar uma nova figura:

```python
plt.figure()
plt.plot(x, y)
```

Ou fechar figuras antigas:

```python
plt.close("all")
```

### Matplotlib.pyplot

| Função | Uso |
|---|---|
| `plt.plot()` | Gráfico de linha |
| `plt.scatter()` | Gráfico de dispersão |
| `plt.bar()` | Barras verticais |
| `plt.barh()` | Barras horizontais |
| `plt.hist()` | Histograma |
| `plt.pie()` | Gráfico de pizza |
| `plt.subplots()` | Cria figura e eixos |
| `plt.figure()` | Cria uma figura |
| `plt.title()` | Define título |
| `plt.xlabel()` / `plt.ylabel()` | Define rótulos dos eixos |
| `plt.legend()` | Exibe legenda |
| `plt.grid()` | Exibe grade |
| `plt.xlim()` / `plt.ylim()` | Define limites dos eixos |
| `plt.xticks()` / `plt.yticks()` | Configura marcações |
| `plt.xscale()` / `plt.yscale()` | Define escala dos eixos |
| `plt.colorbar()` | Adiciona barra de cores |
| `plt.axhline()` | Linha horizontal de referência |
| `plt.axvline()` | Linha vertical de referência |
| `plt.tight_layout()` | Ajusta o espaçamento |
| `plt.savefig()` | Salva o gráfico |
| `plt.show()` | Exibe o gráfico |
| `plt.close()` | Fecha a figura |

---

## Boas práticas

1. Identifique os eixos com `set_xlabel()` e `set_ylabel()`.
2. Use `label` e `legend()` quando houver mais de uma série de dados.
3. Use `tight_layout()` antes de salvar ou exibir gráficos complexos.
4. Salve figuras importantes com resolução apropriada usando `dpi=300`.
5. Evite gráficos de pizza quando houver muitas categorias; barras geralmente facilitam a comparação.
6. Prefira `fig, ax = plt.subplots()` em códigos maiores, pois essa forma facilita controlar cada eixo.

## Referências oficiais

- [Documentação do Matplotlib](https://matplotlib.org/stable/)
- [Referência do `matplotlib.pyplot`](https://matplotlib.org/stable/api/pyplot_summary.html)

