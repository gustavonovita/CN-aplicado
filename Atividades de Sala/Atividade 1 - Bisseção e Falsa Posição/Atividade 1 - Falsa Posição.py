# Método da Falsa posição
H_u = 1
H_l = 0.5
tol = 1e-6
max_iter = 100


def f(H):
    numerador = (20*H)**(5/3)
    denominador = (20 + 2*H)**(2/3)
    return 0.471405 * (numerador / denominador) - 5


for k in range(1, max_iter + 1):
    H_r = H_u - (f(H_u) * ((H_u - H_l) / (f(H_u) - f(H_l))))

    if abs(f(H_r)) < tol:
        print(f"Raiz encontrada: H = {H_r:.6f} após {k} iterações.")
        break

    print(f"Iteração {k}: Limite inferior = {H_l:.6f}, Limite superior = {H_u:.6f}, Hr = {H_r:.6f}, f(Hr) = {f(H_r):.6f}")

    if f(H_l) * f(H_r) < 0:
        H_u = H_r
    else:
        H_l = H_r