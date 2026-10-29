def newton_raphson(f, df, x0, tol=1e-6, max_it=100):
    x_atual = x0
    fx = f(x_atual)
    if abs(fx) < tol:
        return x_atual, 0

    for it in range(1, max_it + 1):
        dfx = df(x_atual)
        if abs(dfx) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: derivada f'(x) muito próxima de zero.")

        x_novo = x_atual - fx / dfx
        fx_novo = f(x_novo)
        erro = abs(x_novo - x_atual)

        if erro < tol or abs(fx_novo) < tol:
            return x_novo, it

        x_atual, fx = x_novo, fx_novo

    return x_atual, max_it