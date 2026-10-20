"""
Analysis of Algorithms - Assignment 07
Exercise 1: Forward elimination (Gaussian elimination)
Exercise 2: Backward substitution and its Theta(n^2) running time
"""
import random


def forward_elimination(A, b):
    """
    Transform the augmented system [A | b] into an upper triangular system [U | c].
    Uses partial pivoting (swap with the row that has the largest pivot) for stability.
    Returns new matrices; the input is not modified.
    """
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]  # augmented matrix
    for k in range(n - 1):
        # Partial pivoting
        p = max(range(k, n), key=lambda r: abs(M[r][k]))
        if abs(M[p][k]) < 1e-12:
            raise ValueError("Matrix is singular")
        M[k], M[p] = M[p], M[k]
        # Eliminate the entries below the pivot
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]
    U = [row[:n] for row in M]
    c = [row[n] for row in M]
    return U, c


def backward_substitution(U, c, counter=None):
    """
    Solve U x = c where U is upper triangular.
    x[i] = (c[i] - sum_{j=i+1}^{n-1} U[i][j] * x[j]) / U[i][i],  for i = n-1, ..., 0
    """
    n = len(U)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = c[i]
        for j in range(i + 1, n):
            s -= U[i][j] * x[j]
            if counter is not None:
                counter[0] += 1        # one multiplication-subtraction
        x[i] = s / U[i][i]
        if counter is not None:
            counter[0] += 1            # one division
    return x


def gauss_solve(A, b):
    U, c = forward_elimination(A, b)
    return backward_substitution(U, c)


def mat_vec(A, x):
    return [sum(a * xi for a, xi in zip(row, x)) for row in A]


def print_matrix(M, name):
    print(f"{name} =")
    for row in M:
        print("   [" + "  ".join(f"{v:8.3f}" for v in row) + " ]")


def random_system(n):
    A = [[random.uniform(-10, 10) for _ in range(n)] for _ in range(n)]
    x_true = [random.uniform(-10, 10) for _ in range(n)]
    return A, mat_vec(A, x_true), x_true


if __name__ == "__main__":
    print("EXERCISE 1: Forward elimination")
    A = [[2.0, 1.0, -1.0],
         [-3.0, -1.0, 2.0],
         [-2.0, 1.0, 2.0]]
    b = [8.0, -11.0, -3.0]
    print_matrix(A, "A")
    print("b =", b)
    U, c = forward_elimination(A, b)
    print("\nAfter forward elimination (upper triangular system):")
    print_matrix(U, "U")
    print("c =", [round(v, 3) for v in c])

    print()
    print("EXERCISE 2: Backward substitution")
    x = backward_substitution(U, c)
    print("x =", [round(v, 6) for v in x])
    print("Expected solution: [2, 3, -1]")
    print("Check A x =", [round(v, 6) for v in mat_vec(A, x)], " b =", b)

    random.seed(7)
    ok = 0
    for _ in range(200):
        n = random.randint(2, 30)
        A2, b2, xt = random_system(n)
        xs = gauss_solve(A2, b2)
        if max(abs(u - v) for u, v in zip(xs, xt)) < 1e-6:
            ok += 1
    print(f"\n{ok}/200 random systems solved correctly")

    print("\nRunning time of backward substitution (operation count):")
    print(f"{'n':>6} {'operations':>12} {'n(n+1)/2':>10} {'ops/n^2':>8}")
    for n in [10, 50, 100, 200, 400, 800]:
        A2, b2, _ = random_system(n)
        U2, c2 = forward_elimination(A2, b2)
        counter = [0]
        backward_substitution(U2, c2, counter)
        print(f"{n:>6} {counter[0]:>12} {n*(n+1)//2:>10} {counter[0]/n**2:>8.3f}")
    print("-> operations = n(n+1)/2 exactly, and ops/n^2 tends to 1/2, so T(n) is in Theta(n^2)")
