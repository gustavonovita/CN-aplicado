def eliminacao_gaussiana_pivoteamento_parcial(A, b, tol=1e-15):
    n = len(b)
    M = [[float(val) for val in linha] for linha in A]
    v = [float(val) for val in b]

    # Eliminação progressiva com pivoteamento parcial
    for k in range(n - 1):
        # Localiza a linha p com o maior elemento em módulo na coluna k
        p = max(range(k, n), key=lambda i: abs(M[i][k]))

        if abs(M[p][k]) < tol:
            raise ZeroDivisionError(f"Matriz singular ou quase singular na coluna {k}.")

        # Permuta as linhas k e p em M e em v, se necessário
        if p != k:
            M[k], M[p] = M[p], M[k]
            v[k], v[p] = v[p], v[k]

        # Elimina os elementos abaixo do pivô
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k + 1, n):
                M[i][j] -= m * M[k][j]
            v[i] -= m * v[k]
            M[i][k] = 0.0

    if abs(M[n - 1][n - 1]) < tol:
        raise ZeroDivisionError("Matriz singular: último pivô nulo ou muito próximo de zero.")

    # Substituição retroativa
    x = [0.0] * n
    x[n - 1] = v[n - 1] / M[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = v[i]
        for j in range(i + 1, n):
            soma -= M[i][j] * x[j]
        x[i] = soma / M[i][i]

    return x