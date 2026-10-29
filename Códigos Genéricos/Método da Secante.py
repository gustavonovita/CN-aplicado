def secante(f, x0, x1, tol=1e-6, max_it=100):
    fx0, fx1 = f(x0), f(x1)
    if abs(fx0) < tol:
        return x0, 0
    if abs(fx1) < tol:
        return x1, 0

    for it in range(1, max_it + 1):
        if abs(fx1 - fx0) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: f(x1) e f(x0) são muito próximos.")

        x_novo = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        fx_novo = f(x_novo)
        erro = abs(x_novo - x1)

        if erro < tol or abs(fx_novo) < tol:
            return x_novo, it

        x0, fx0 = x1, fx1
        x1, fx1 = x_novo, fx_novo

    return x1, max_it