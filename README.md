# Analysis of Algorithms – Assignment 07

Yachay Tech – School of Mathematical and Computational Sciences

## Contents

| File | Exercise | Description |
|---|---|---|
| `gaussian_elimination.py` | 1, 2 | Forward elimination (with partial pivoting) and backward substitution, with an operation count showing Θ(n²) |
| `heap.py` | 3 | Bottom-up construction of a max-heap from a random array |
| `lu_decomposition.py` | 4 | Solving A x = b with LU decomposition (P A = L U) obtained from Gaussian elimination |

The Θ(n²) analysis of backward substitution (Exercise 2) is also included in the PDF submission.

## Requirements

- Python 3.8+
- No external libraries are needed.

## How to run

```bash
python gaussian_elimination.py
python heap.py
python lu_decomposition.py
```

## Summary of results

**Exercise 1 – Forward elimination.** Transforms the augmented system [A | b] into an upper triangular system [U | c]. Partial pivoting is used to avoid division by zero and improve numerical stability.

**Exercise 2 – Backward substitution.** Solves U x = c from the last row to the first: x[i] = (c[i] − Σ U[i][j]·x[j]) / U[i][i]. Row i needs n − 1 − i multiply-subtract operations plus one division, so the total is n(n+1)/2 operations, which is Θ(n²). The program counts the operations and confirms that operations / n² → 1/2.

**Exercise 3 – Heap construction.** Starting from the last internal node (index n/2 − 1) down to the root, each node is sifted down until the max-heap property holds. This runs in Θ(n). The program verifies the result on 1000 random arrays.

**Exercise 4 – LU decomposition.** During Gaussian elimination, each multiplier used to eliminate an entry is stored in L, and the resulting upper triangular matrix is U, so that P A = L U. The system is then solved in two triangular steps: L y = P b (forward substitution) and U x = y (backward substitution). The same factorization can be reused for several right-hand sides.
