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
                
                T[i, j] = (T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1]) / 4.0
                
                diferenca = abs(T[i, j] - antigo)
                erro = max(erro, diferenca)
                
        iteracao += 1

    max_residuo = 0.0
    for j in range(1, n + 1):
        for i in range(1, n + 1):
            r_ij = 4 * T[i, j] - (T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1])
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