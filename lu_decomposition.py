"""
Analysis of Algorithms - Assignment 07
Exercise 4: Solve a system of equations using LU decomposition based on Gaussian elimination
"""
import random


def lu_decomposition(A):
    """
    Compute P A = L U using Gaussian elimination with partial pivoting.
    - U is the upper triangular matrix produced by forward elimination.
    - L stores the elimination factors (multipliers) below a diagonal of ones.
    - P is the row permutation, stored as a list of row indices.
    """
    n = len(A)
    U = [row[:] for row in A]
    L = [[0.0] * n for _ in range(n)]
    P = list(range(n))
    for k in range(n - 1):
        p = max(range(k, n), key=lambda r: abs(U[r][k]))
        if abs(U[p][k]) < 1e-12:
            raise ValueError("Matrix is singular")
        if p != k:
            U[k], U[p] = U[p], U[k]
            L[k], L[p] = L[p], L[k]
            P[k], P[p] = P[p], P[k]
        for i in range(k + 1, n):
            factor = U[i][k] / U[k][k]
            L[i][k] = factor                 # the multiplier goes into L
            for j in range(k, n):
                U[i][j] -= factor * U[k][j]  # same step as forward elimination
    for i in range(n):
        L[i][i] = 1.0
    return P, L, U


def forward_substitution(L, b):
    """Solve L y = b (L lower triangular with ones on the diagonal)."""
    n = len(L)
    y = [0.0] * n
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))
    return y


def backward_substitution(U, y):
    """Solve U x = y (U upper triangular)."""
    n = len(U)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return x


def lu_solve(P, L, U, b):
    """Solve A x = b using P A = L U:  L y = P b,  then  U x = y."""
    Pb = [b[i] for i in P]
    y = forward_substitution(L, Pb)
    return backward_substitution(U, y)


def mat_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def mat_vec(A, x):
    return [sum(a * xi for a, xi in zip(row, x)) for row in A]


def print_matrix(M, name):
    print(f"{name} =")
    for row in M:
        print("   [" + "  ".join(f"{v:8.3f}" for v in row) + " ]")


if __name__ == "__main__":
    print("EXERCISE 4: Solving a system with LU decomposition")
    A = [[2.0, -1.0, -2.0],
         [-4.0, 6.0, 3.0],
         [-4.0, -2.0, 8.0]]
    b = [-1.0, 13.0, -6.0]
    print_matrix(A, "A")
    print("b =", b)

    P, L, U = lu_decomposition(A)
    print("\nPermutation P (row order):", P)
    print_matrix(L, "L")
    print_matrix(U, "U")
    PA = [A[i] for i in P]
    LU = mat_mul(L, U)
    err = max(abs(PA[i][j] - LU[i][j]) for i in range(3) for j in range(3))
    print("Check P A = L U -> max error:", f"{err:.2e}")

    x = lu_solve(P, L, U, b)
    print("\nSolution x =", [round(v, 6) for v in x])
    print("Check A x =", [round(v, 6) for v in mat_vec(A, x)], " b =", b)

    # Reuse the same LU factorization for several right-hand sides
    print("\nReusing the same LU for other right-hand sides:")
    for b2 in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [5.0, 5.0, 5.0]):
        x2 = lu_solve(P, L, U, b2)
        print("   b =", b2, "-> x =", [round(v, 4) for v in x2])

    random.seed(4)
    ok = 0
    for _ in range(200):
        n = random.randint(2, 30)
        A2 = [[random.uniform(-10, 10) for _ in range(n)] for _ in range(n)]
        xt = [random.uniform(-10, 10) for _ in range(n)]
        P2, L2, U2 = lu_decomposition(A2)
        xs = lu_solve(P2, L2, U2, mat_vec(A2, xt))
        if max(abs(u - v) for u, v in zip(xs, xt)) < 1e-6:
            ok += 1
    print(f"\n{ok}/200 random systems solved correctly with LU")
