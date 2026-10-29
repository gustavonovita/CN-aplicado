def bissecao(f, a, b, tol=1e-6, max_it=100):
    fa, fb = f(a), f(b)
    if abs(fa) < tol:
        return a, 0
    if abs(fb) < tol:
        return b, 0
    if fa * fb > 0:
        raise ValueError(f"O intervalo [{a}, {b}] não apresenta mudança de sinal.")

    xr_ant = a
    for it in range(1, max_it + 1):
        xr = (a + b) / 2.0
        fxr = f(xr)
        erro = abs(xr - xr_ant) if it > 1 else (b - a) / 2.0

        if abs(fxr) < tol or (b - a) / 2.0 < tol or (it > 1 and erro < tol):
            return xr, it

        if fa * fxr < 0:
            b, fb = xr, fxr
        else:
            a, fa = xr, fxr

        xr_ant = xr

    return xr, max_it