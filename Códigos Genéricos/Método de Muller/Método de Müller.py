def muller(f, x0, x1, x2, tol=1e-6, max_it=100):
    fx0, fx1, fx2 = f(x0), f(x1), f(x2)
    if abs(fx0) < tol:
        return x0, 0
    if abs(fx1) < tol:
        return x1, 0
    if abs(fx2) < tol:
        return x2, 0

    for it in range(1, max_it + 1):
        h0 = x1 - x0
        h1 = x2 - x1
        if abs(h0) < 1e-15 or abs(h1) < 1e-15 or abs(h1 + h0) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: pontos muito próximos.")

        d0 = (fx1 - fx0) / h0
        d1 = (fx2 - fx1) / h1

        a = (d1 - d0) / (h1 + h0)
        b = a * h1 + d1
        c = fx2

        raiz_disc = (b**2 - 4 * a * c) ** 0.5
        den_pos = b + raiz_disc
        den_neg = b - raiz_disc
        denominador = den_pos if abs(den_pos) >= abs(den_neg) else den_neg

        if abs(denominador) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: denominador nulo na fórmula quadrática.")

        dx = -2 * c / denominador
        x3 = x2 + dx
        fx3 = f(x3)
        erro = abs(dx)

        if erro < tol or abs(fx3) < tol:
            return x3, it

        x0, fx0 = x1, fx1
        x1, fx1 = x2, fx2
        x2, fx2 = x3, fx3

    return x2, max_it