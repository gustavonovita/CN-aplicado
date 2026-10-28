import numpy as np
import matplotlib.pyplot as plt

# Parâmetros iniciais[cite: 14]
b = 4.0     
z = 1.5   
g = 9.81  
max_iter = 100
tol = 1e-8

def polinomio_f(H, Q):
    a6 = g * (z**3) * (H**6)
    a5 = 3 * g * b * (z**2) * (H**5) 
    a4 = 3 * g * (b**2) * z * (H**4)
    a3 = g * (b**3) * (H**3) 
    a1 = -2 * (Q**2) * z * H   
    a0 = -(Q**2) * b           
    return a1 + a3 + a4 + a5 + a6 + a0  

def metodo_muller(Q, H0_val, H1_val, H2_val, imprimir_passos=False):
    H0 = float(H0_val)
    H1 = float(H1_val)
    H2 = float(H2_val)

    for k in range(1, max_iter + 1):
        f0 = polinomio_f(H0, Q)
        f1 = polinomio_f(H1, Q)
        f2 = polinomio_f(H2, Q)

        h0 = H1 - H0
        h1 = H2 - H1

        delta_0 = (f1 - f0) / h0
        delta_1 = (f2 - f1) / h1

        a = (delta_1 - delta_0) / (h1 + h0)
        b_m = a * h1 + delta_1
        c = f2

        D = np.emath.sqrt(b_m**2 - 4 * a * c)

        if abs(b_m + D) > abs(b_m - D):
            E = b_m + D
        else:
            E = b_m - D

        Delta_H = - (2 * c / E)
        H3 = H2 + Delta_H

        erro = abs(H3 - H2)
        f_H3 = polinomio_f(H3, Q).real

        # Imprime exatamente as 5 exigências contidas na imagem durante a execução[cite: 1]
        if imprimir_passos:
            print(f"Iteração k={k:02d} | Aproximações: H0={H0.real:.4f}, H1={H1.real:.4f}, H2={H2.real:.4f} | H3={H3.real:.6f} | f(H3)={f_H3:.4e} | erro ={erro:.4e}")

        if erro < tol:
            return H3.real

        H0, H1, H2 = H1, H2, H3
        
    return H2.real

# =========================================================
# Execução Principal
# =========================================================
if __name__ == "__main__":
    # Validação inicial exigida para Q = 20 m³/s com impressão detalhada ativada[cite: 14]
    print("--- Validação do Método (Q = 20 m³/s) ---")
    Hc_valido = metodo_muller(20.0, 0.5, 1.0, 1.5, imprimir_passos=True)
    print(f"-> Raiz encontrada: Hc = {Hc_valido:.6f} m\n")

    # Parte 5: Varredura de vazões, arquivo e gráfico[cite: 14]
    vazoes = np.array([5, 10, 15, 20, 25, 30, 35, 40], dtype=float)
    profundidades_criticas = []

    print("--- Calculando Profundidades Críticas (Varredura) ---")
    for Q_atual in vazoes:
        # Chama a função de forma silenciosa para não poluir a tela
        Hc = metodo_muller(Q_atual, 0.5, 1.0, 1.5, imprimir_passos=False)
        profundidades_criticas.append(Hc)
        print(f"Vazão Q = {Q_atual:02.0f} m³/s -> Hc = {Hc:.6f} m")

    # Salvar resultados em um arquivo .txt de duas colunas (Q e Hc)[cite: 14]
    dados_exportacao = np.column_stack((vazoes, profundidades_criticas))
    np.savetxt("resultados_Q_Hc.txt", dados_exportacao, fmt="%.2f %.6f", header="Q(m3/s) Hc(m)")
    print("\nResultados numéricos salvos no arquivo 'resultados_Q_Hc.txt'.")

    # Construir e salvar o gráfico Hc x Q[cite: 14]
    plt.figure(figsize=(8, 5))
    plt.plot(vazoes, profundidades_criticas, marker='o', linestyle='-', color='b')
    plt.title("Variação da Profundidade Crítica com a Vazão")
    plt.xlabel("Vazão Volumétrica Q (m³/s)")
    plt.ylabel("Profundidade Crítica Hc (m)")
    plt.grid(True)
    plt.savefig("grafico_Hc_x_Q.png")
    
    print("Gráfico finalizado e salvo como 'grafico_Hc_x_Q.png'.")