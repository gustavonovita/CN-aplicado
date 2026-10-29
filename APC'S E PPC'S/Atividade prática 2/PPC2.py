import numpy as np
import matplotlib.pyplot as plt


def bairstow_step(a, r, s, tol=1e-8, max_iter=100, blowup=1e12):
    a = np.asarray(a, dtype=float)
    n = len(a) - 1

    b = np.zeros(n + 1)
    c = np.zeros(n + 1)

    for iter_count in range(1, max_iter + 1):
        b[n] = a[n]
        b[n - 1] = a[n - 1] + r * b[n]
        for i in range(n - 2, -1, -1):
            b[i] = a[i] + r * b[i + 1] + s * b[i + 2]

        c[n] = b[n]
        c[n - 1] = b[n - 1] + r * c[n]
        for i in range(n - 2, 0, -1):
            c[i] = b[i] + r * c[i + 1] + s * c[i + 2]

        det = c[1] * c[3] - c[2] ** 2

        if abs(det) < 1e-300:
            r += 1e-3
            s += 1e-3
            continue

        dr = (-b[0] * c[3] + c[2] * b[1]) / det
        ds = (-c[1] * b[1] + c[2] * b[0]) / det

        r_new = r + dr
        s_new = s + ds

        if (not np.isfinite(r_new)) or (not np.isfinite(s_new)) \
                or abs(r_new) > blowup or abs(s_new) > blowup:
            return r, s, b, iter_count, False

        tol_r = abs(dr / r_new) if r_new != 0 else abs(dr)
        tol_s = abs(ds / s_new) if s_new != 0 else abs(ds)

        r, s = r_new, s_new

        if tol_r < tol and tol_s < tol:
            b[n] = a[n]
            b[n - 1] = a[n - 1] + r * b[n]
            for i in range(n - 2, -1, -1):
                b[i] = a[i] + r * b[i + 1] + s * b[i + 2]
            return r, s, b, iter_count, True

    return r, s, b, max_iter, False


def get_roots_from_quadratic(r, s):
    delta = r ** 2 + 4 * s
    if delta > 0:
        sq = np.sqrt(delta)
        return [(r + sq) / 2, (r - sq) / 2]
    elif delta == 0:
        return [r / 2, r / 2]
    else:
        real_part = r / 2
        imag_part = np.sqrt(-delta) / 2
        return [complex(real_part, imag_part), complex(real_part, -imag_part)]


def bairstow_full(a, r0=1.0, s0=1.0, tol=1e-8, max_iter=200, verbose=False):
    roots = []
    a_curr = np.array(a, dtype=float)
    r_guess, s_guess = r0, s0

    while len(a_curr) > 3:
        r, s, b, iters, ok = bairstow_step(a_curr, r_guess, s_guess, tol, max_iter)
        if not ok:
            raise RuntimeError(
                f"Bairstow nao convergiu a partir de r0={r0}, s0={s0} "
                f"(grau atual = {len(a_curr) - 1}, iteracoes = {iters})."
            )
        novas_raizes = get_roots_from_quadratic(r, s)
        roots.extend(novas_raizes)
        if verbose:
            print(f"  [grau {len(a_curr)-1} -> {len(a_curr)-3}] "
                  f"r={r:.6f}  s={s:.6f}  iters={iters}  "
                  f"raizes: {novas_raizes[0]:.6f}, {novas_raizes[1]:.6f}")
        a_curr = b[2:]
        r_guess, s_guess = r, s

    if len(a_curr) == 3:
        r_final = -a_curr[1] / a_curr[2]
        s_final = -a_curr[0] / a_curr[2]
        roots.extend(get_roots_from_quadratic(r_final, s_final))
    elif len(a_curr) == 2:
        roots.append(-a_curr[0] / a_curr[1])

    return roots


def analisar_convergencia(a, r0_vals, s0_vals, tol=1e-8, max_iter=100):
    print(f"{'r0':>8s} {'s0':>8s}   {'status':>10s}   {'iteracoes':>10s}")

    matriz_iter = np.full((len(s0_vals), len(r0_vals)), np.nan)
    n_ok = 0
    n_total = len(r0_vals) * len(s0_vals)

    for i, s0 in enumerate(s0_vals):
        for j, r0 in enumerate(r0_vals):
            _, _, _, iters, ok = bairstow_step(a, r0, s0, tol=tol, max_iter=max_iter)
            if ok:
                matriz_iter[i, j] = iters
                n_ok += 1
                print(f"{r0:8.3f} {s0:8.3f}   {'convergiu':>10s}   {iters:10d}")
            else:
                print(f"{r0:8.3f} {s0:8.3f}   {'DIVERGIU':>10s}   {'--':>10s}")

    print(f"\nTaxa de convergencia: {100*n_ok/n_total:.1f}%  "
          f"({n_ok}/{n_total} chutes testados)")
    if n_ok > 0:
        print(f"Media de iteracoes (apenas casos convergidos): "
              f"{np.nanmean(matriz_iter):.2f}")

    return matriz_iter


def plotar_mapa_convergencia(matriz_iter, r0_vals, s0_vals,
                              titulo="Convergencia do Bairstow vs. chute inicial",
                              arquivo="mapa_convergencia.png"):
    plt.figure(figsize=(7, 6))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="black")
    plt.imshow(matriz_iter, origin="lower", aspect="auto", cmap=cmap,
               extent=[r0_vals.min(), r0_vals.max(),
                       s0_vals.min(), s0_vals.max()])
    plt.colorbar(label="Iteracoes ate a convergencia (preto = divergiu)")
    plt.xlabel(r"$r_0$")
    plt.ylabel(r"$s_0$")
    plt.title(titulo)
    plt.tight_layout()
    plt.savefig(arquivo, dpi=150)
    plt.show()


def gerar_fractal_bairstow(a, resolucao=300, r_lim=(-3, 3), s_lim=(-3, 3),
                            max_iter=50, arquivo="fractal_bairstow.png"):
    r_vals = np.linspace(r_lim[0], r_lim[1], resolucao)
    s_vals = np.linspace(s_lim[0], s_lim[1], resolucao)

    fractal = np.zeros((resolucao, resolucao))

    for i, s0 in enumerate(s_vals):
        for j, r0 in enumerate(r_vals):
            _, _, _, iters, ok = bairstow_step(a, r0, s0, max_iter=max_iter)
            fractal[i, j] = iters if ok else 0

    plt.figure(figsize=(8, 6))
    plt.imshow(fractal, extent=[r_lim[0], r_lim[1], s_lim[0], s_lim[1]],
               origin='lower', cmap='nipy_spectral')
    plt.colorbar(label='Iteracoes para Convergencia (0 = divergiu)')
    plt.title('Fractal de Bairstow - Polinomio da APC2')
    plt.xlabel(r'$r_0$ inicial')
    plt.ylabel(r'$s_0$ inicial')
    plt.tight_layout()
    plt.savefig(arquivo, dpi=150)
    plt.show()


def interpretar_autovalores(autovalores):
    print("\n--- Interpretacao fisica dos autovalores ---")
    todos_estaveis = True
    algum_oscilatorio = False

    for lam in autovalores:
        lam = complex(lam)
        sigma = lam.real
        omega_d = lam.imag
        omega_n = abs(lam)
        zeta = -sigma / omega_n if omega_n > 1e-12 else None

        if sigma >= -1e-9:
            todos_estaveis = False
        if abs(omega_d) > 1e-9:
            algum_oscilatorio = True

        print(f"  lambda = {lam:.6f}")
        print(f"      sigma (decaimento)     = {sigma:.6f} 1/s")
        if abs(omega_d) > 1e-9:
            periodo = 2 * np.pi / abs(omega_d)
            print(f"      omega_d (oscilacao)    = {omega_d:.6f} rad/s")
            print(f"      omega_n = |lambda|     = {omega_n:.6f} rad/s")
            print(f"      zeta (amortecimento)   = {zeta:.4f}")
            print(f"      periodo de oscilacao   = {periodo:.4f} s")
        else:
            print("      modo nao oscilatorio (raiz real)")
        print()

    print(f"Estabilidade global : sistema "
          f"{'ESTAVEL' if todos_estaveis else 'INSTAVEL'} "
          f"(todas as partes reais negativas -> vibracoes decaem)."
          if todos_estaveis else
          "Estabilidade global : sistema INSTAVEL "
          "(ha pelo menos um polo com parte real >= 0).")
    if algum_oscilatorio:
        print("Comportamento       : resposta OSCILATORIA (subamortecida), "
              "evidenciada pela parte imaginaria nao nula dos autovalores.")


if __name__ == "__main__":

    print("=" * 70)
    print("VALIDACAO: Polinomio de 7a Ordem")
    print("=" * 70)
    raizes_conhecidas = [1, -1, 2, -2, 3, complex(-1, 1), complex(-1, -1)]
    poly_7_numpy = np.poly(raizes_conhecidas)
    a_validacao = poly_7_numpy[::-1].real

    raizes_encontradas = bairstow_full(a_validacao, r0=0.5, s0=0.5, verbose=True)

    print("\nRaizes conhecidas originais:")
    for raiz in sorted(raizes_conhecidas, key=lambda x: x.real if isinstance(x, complex) else x):
        print(f"  {np.round(raiz, 4)}")

    print("\nRaizes encontradas pelo metodo de Bairstow:")
    for raiz in sorted(raizes_encontradas, key=lambda x: x.real if isinstance(x, complex) else x):
        print(f"  {np.round(raiz, 4)}")

    conhecidas_c = sorted([complex(z) for z in raizes_conhecidas], key=lambda z: (round(z.real, 6), round(z.imag, 6)))
    encontradas_c = sorted([complex(z) for z in raizes_encontradas], key=lambda z: (round(z.real, 6), round(z.imag, 6)))
    erro_max = max(abs(k - f) for k, f in zip(conhecidas_c, encontradas_c))
    print(f"\nErro absoluto maximo entre raiz conhecida e encontrada: {erro_max:.3e}")
    print("Validacao:", "APROVADA" if erro_max < 1e-6 else "REPROVADA")

    print("\n" + "=" * 70)
    print("ANALISE DE CONVERGENCIA PARA DIFERENTES (r0, s0)")
    print("=" * 70)
    r0_vals = np.linspace(-3, 3, 9)
    s0_vals = np.linspace(-3, 3, 9)
    matriz_iter = analisar_convergencia(a_validacao, r0_vals, s0_vals)
    plotar_mapa_convergencia(matriz_iter, r0_vals, s0_vals,
                              arquivo="mapa_convergencia_validacao.png")

    print("\n" + "=" * 70)
    print("SOLUCAO DO POLINOMIO CARACTERISTICO DA APC2")
    print("=" * 70)
    a_apc2 = [100.0, 40.0, 34.0, 6.0, 1.0]

    autovalores = bairstow_full(a_apc2, r0=1.0, s0=-1.0, verbose=True)

    print("\nAutovalores do sistema (raizes do polinomio caracteristico):")
    for lam in autovalores:
        print(f"  lambda = {np.round(lam, 5)}")

    interpretar_autovalores(autovalores)

    print("\n" + "=" * 70)
    print("FRACTAL DE BAIRSTOW - polinomio da APC2")
    print("=" * 70)
    gerar_fractal_bairstow(a_apc2, resolucao=300, r_lim=(-6, 6), s_lim=(-8, 4),
                            arquivo="fractal_bairstow_apc2.png")