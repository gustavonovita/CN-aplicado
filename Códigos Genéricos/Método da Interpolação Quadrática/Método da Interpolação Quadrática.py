def interpolacao_quadratica(f, x0, x1, x2, tol=1e-6, max_it=100, minimizar=True):
    fx0, fx1, fx2 = f(x0), f(x1), f(x2)
    x3_ant = x1
    x3 = x1

    for it in range(1, max_it + 1):
        num = fx0 * (x1**2 - x2**2) + fx1 * (x2**2 - x0**2) + fx2 * (x0**2 - x1**2)
        den = 2.0 * (fx0 * (x1 - x2) + fx1 * (x2 - x0) + fx2 * (x0 - x1))

        if abs(den) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: os três pontos são colineares.")

        x3 = num / den
        fx3 = f(x3)
        erro = abs(x3 - x3_ant)

        if it > 1 and erro < tol:
            return x3, it

        melhor = (fx3 < fx1) if minimizar else (fx3 > fx1)
        if x3 > x1:
            if melhor:
                x0, fx0 = x1, fx1
                x1, fx1 = x3, fx3
            else:
                x2, fx2 = x3, fx3
        else:
            if melhor:
                x2, fx2 = x1, fx1
                x1, fx1 = x3, fx3
            else:
                x0, fx0 = x3, fx3

        x3_ant = x3

    return x3, max_it