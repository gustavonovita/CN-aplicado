def eliminacao_gaussiana_ingenua(A, b, tol=1e-15):
    n = len(b)
    M = [[float(val) for val in linha] for linha in A]
    v = [float(val) for val in b]

    # Eliminação progressiva
    for k in range(n - 1):
        if abs(M[k][k]) < tol:
            raise ZeroDivisionError(f"Pivô nulo ou muito próximo de zero na linha {k}.")
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k + 1, n):
                M[i][j] -= m * M[k][j]
            v[i] -= m * v[k]
            M[i][k] = 0.0

    if abs(M[n - 1][n - 1]) < tol:
        raise ZeroDivisionError("Último pivô nulo ou muito próximo de zero.")

    # Substituição retroativa
    x = [0.0] * n
    x[n - 1] = v[n - 1] / M[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = v[i]
        for j in range(i + 1, n):
            soma -= M[i][j] * x[j]
        x[i] = soma / M[i][i]

    return x