import numpy as np
import matplotlib.pyplot as plt

# ==================================================
# Função da equação
# ==================================================
def f(H):
    return (0.471405 * (20 * H)**(5/3)) / ((20 + 2 * H)**(2/3)) - 5


# ==================================================
# Método da Bissecção
# ==================================================
a = 0.0
b = 1.5

tol = 1e-6
max_iter = 100

erros = []
iteracoes = []

xr_ant = None

print("Iteração |       H       |     Erro (%)")
print("-" * 40)

for i in range(1, max_iter + 1):

    xr = (a + b) / 2

    if xr_ant is not None:
        erro = abs((xr - xr_ant) / xr) * 100
    else:
        erro = 100.0

    iteracoes.append(i)
    erros.append(erro)

    print(f"{i:8d} | {xr:12.6f} | {erro:12.6f}")

    if f(a) * f(xr) < 0:
        b = xr
    else:
        a = xr

    if erro < tol:
        break

    xr_ant = xr

# ==================================================
# Resultado
# ==================================================
print("\nRaiz encontrada:")
print(f"H = {xr:.6f} m")

# ==================================================
# Salvar arquivo de saída
# ==================================================
with open("saida.txt", "w", encoding="utf-8") as arq:
    arq.write("Iteracao\tH\tErro(%)\n")

    for n, e in zip(iteracoes, erros):
        arq.write(f"{n}\t{xr:.8f}\t{e:.8f}\n")

    arq.write("\n")
    arq.write(f"Raiz encontrada: H = {xr:.8f} m\n")

# ==================================================
# Gráfico Erro × Iteração
# ==================================================
plt.figure(figsize=(8,5))
plt.semilogy(iteracoes, erros, 'bo-')
plt.grid(True)
plt.xlabel("Número de Iterações")
plt.ylabel("Erro Relativo (%)")
plt.title("Método da Bissecção")
plt.tight_layout()
plt.show()