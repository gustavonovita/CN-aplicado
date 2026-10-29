def gauss_seidel(A, b, x0=None, tol=1e-6, max_it=100):
    n = len(b)
    x = [0.0] * n if x0 is None else [float(val) for val in x0]

    for i in range(n):
        if abs(A[i][i]) < 1e-15:
            raise ZeroDivisionError(f"Elemento diagonal nulo ou muito próximo de zero na linha {i}.")

    for it in range(1, max_it + 1):
        maior_erro = 0.0

        for i in range(n):
            x_ant = x[i]
            soma = float(b[i])

            for j in range(i):
                soma -= A[i][j] * x[j]
            for j in range(i + 1, n):
                soma -= A[i][j] * x[j]

            x[i] = soma / A[i][i]
            erro = abs(x[i] - x_ant)

            if erro > maior_erro:
                maior_erro = erro

        if maior_erro < tol:
            return x, it

    return x, max_it