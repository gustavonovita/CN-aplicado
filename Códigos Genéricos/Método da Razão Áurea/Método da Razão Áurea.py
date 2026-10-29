def razao_aurea(f, a, b, tol=1e-6, max_it=100, minimizar=True):
    phi = (5.0 ** 0.5 - 1.0) / 2.0
    x1 = b - phi * (b - a)
    x2 = a + phi * (b - a)
    fx1, fx2 = f(x1), f(x2)

    for it in range(1, max_it + 1):
        if (fx1 < fx2) if minimizar else (fx1 > fx2):
            b = x2
            x2, fx2 = x1, fx1
            x1 = b - phi * (b - a)
            fx1 = f(x1)
        else:
            a = x1
            x1, fx1 = x2, fx2
            x2 = a + phi * (b - a)
            fx2 = f(x2)

        if abs(b - a) < tol:
            return (a + b) / 2.0, it

    return (a + b) / 2.0, max_it