def inversao_lu(A, tol=1e-15):
    n = len(A)
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[float(val) for val in linha] for linha in A]

    # Decomposição A = L * U (sem pivoteamento)
    for k in range(n - 1):
        if abs(U[k][k]) < tol:
            raise ZeroDivisionError(f"Pivô nulo ou muito próximo de zero na linha {k}.")
        for i in range(k + 1, n):
            fator = U[i][k] / U[k][k]
            L[i][k] = fator
            U[i][k] = 0.0
            for j in range(k + 1, n):
                U[i][j] -= fator * U[k][j]

    if abs(U[n - 1][n - 1]) < tol:
        raise ZeroDivisionError("Matriz singular: último pivô nulo ou muito próximo de zero.")

    # Resolução das n colunas da matriz inversa
    A_inv = [[0.0] * n for _ in range(n)]
    for c in range(n):
        # Substituição progressiva: L * d = e_c
        d = [0.0] * n
        for i in range(n):
            soma = 1.0 if i == c else 0.0
            for j in range(i):
                soma -= L[i][j] * d[j]
            d[i] = soma

        # Substituição retroativa: U * x_c = d
        for i in range(n - 1, -1, -1):
            soma = d[i]
            for j in range(i + 1, n):
                soma -= U[i][j] * A_inv[j][c]
            A_inv[i][c] = soma / U[i][i]

    return A_inv


def inversao_lu_pivoteamento(A, tol=1e-15):
    n = len(A)
    LU = [[float(val) for val in linha] for linha in A]
    P = list(range(n))

    # Decomposição PA = LU compacta (com pivoteamento parcial)
    for k in range(n - 1):
        p = max(range(k, n), key=lambda i: abs(LU[i][k]))
        if abs(LU[p][k]) < tol:
            raise ZeroDivisionError(f"Matriz singular ou quase singular na coluna {k}.")
        if p != k:
            LU[k], LU[p] = LU[p], LU[k]
            P[k], P[p] = P[p], P[k]

        for i in range(k + 1, n):
            LU[i][k] /= LU[k][k]
            for j in range(k + 1, n):
                LU[i][j] -= LU[i][k] * LU[k][j]

    if abs(LU[n - 1][n - 1]) < tol:
        raise ZeroDivisionError("Matriz singular: último pivô nulo ou muito próximo de zero.")

    # Resolução das n colunas da matriz inversa
    A_inv = [[0.0] * n for _ in range(n)]
    for c in range(n):
        # Substituição progressiva: L * d = P * e_c
        d = [0.0] * n
        for i in range(n):
            soma = 1.0 if P[i] == c else 0.0
            for j in range(i):
                soma -= LU[i][j] * d[j]
            d[i] = soma

        # Substituição retroativa: U * x_c = d
        for i in range(n - 1, -1, -1):
            soma = d[i]
            for j in range(i + 1, n):
                soma -= LU[i][j] * A_inv[j][c]
            A_inv[i][c] = soma / LU[i][i]

    return A_inv