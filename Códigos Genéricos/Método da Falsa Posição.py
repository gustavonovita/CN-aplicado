def falsa_posicao(f, a, b, tol=1e-6, max_it=100):
    fa, fb = f(a), f(b)
    if abs(fa) < tol:
        return a, 0
    if abs(fb) < tol:
        return b, 0
    if fa * fb > 0:
        raise ValueError(f"O intervalo [{a}, {b}] não apresenta mudança de sinal.")

    xr_ant = a
    for it in range(1, max_it + 1):
        if abs(fb - fa) < 1e-15:
            raise ZeroDivisionError("Divisão por zero: f(b) e f(a) são muito próximos.")

        xr = (a * fb - b * fa) / (fb - fa)
        fxr = f(xr)
        erro = abs(xr - xr_ant) if it > 1 else abs(b - a)

        if abs(fxr) < tol or (it > 1 and erro < tol):
            return xr, it

        if fa * fxr < 0:
            b, fb = xr, fxr
        else:
            a, fa = xr, fxr

        xr_ant = xr

    return xr, max_it