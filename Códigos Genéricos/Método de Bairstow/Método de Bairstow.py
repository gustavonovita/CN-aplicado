def bairstow(a, r0=0.1, s0=0.1, tol=1e-6, max_it=100):
    # a = [a0, a1, ..., an] onde a[i] multiplica x^i
    coefs = [float(val) for val in a]
    raizes = []
    r, s = float(r0), float(s0)

    while len(coefs) - 1 >= 3:
        n = len(coefs) - 1
        b = [0.0] * (n + 1)
        c = [0.0] * (n + 1)

        for _ in range(max_it):
            b[n] = coefs[n]
            b[n - 1] = coefs[n - 1] + r * b[n]
            for i in range(n - 2, -1, -1):
                b[i] = coefs[i] + r * b[i + 1] + s * b[i + 2]

            c[n] = b[n]
            c[n - 1] = b[n - 1] + r * c[n]
            for i in range(n - 2, 0, -1):
                c[i] = b[i] + r * c[i + 1] + s * c[i + 2]

            det = c[2] * c[2] - c[1] * c[3]
            if abs(det) < 1e-15:
                r += 0.1
                s += 0.1
                continue

            dr = (-b[1] * c[2] + b[0] * c[3]) / det
            ds = (-b[0] * c[2] + b[1] * c[1]) / det
            r += dr
            s += ds

            if max(abs(dr), abs(ds)) < tol or max(abs(b[0]), abs(b[1])) < tol:
                break

        disc = r * r + 4.0 * s
        if disc >= 0.0:
            raizes.append((r + disc**0.5) / 2.0)
            raizes.append((r - disc**0.5) / 2.0)
        else:
            imag = (-disc) ** 0.5 / 2.0
            raizes.append(complex(r / 2.0, imag))
            raizes.append(complex(r / 2.0, -imag))

        coefs = b[2:]

    grau_restante = len(coefs) - 1
    if grau_restante == 2:
        r_fin = -coefs[1] / coefs[2]
        s_fin = -coefs[0] / coefs[2]
        disc = r_fin * r_fin + 4.0 * s_fin
        if disc >= 0.0:
            raizes.append((r_fin + disc**0.5) / 2.0)
            raizes.append((r_fin - disc**0.5) / 2.0)
        else:
            imag = (-disc) ** 0.5 / 2.0
            raizes.append(complex(r_fin / 2.0, imag))
            raizes.append(complex(r_fin / 2.0, -imag))
    elif grau_restante == 1:
        raizes.append(-coefs[0] / coefs[1])

    return raizes