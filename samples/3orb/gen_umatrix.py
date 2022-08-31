import numpy as np
from itertools import product

U = 1.0
Up = 0.5
J = 0.25
norb = 3

Umat = np.zeros((norb, 2, norb, 2), dtype=float)

for i in range(norb):
    Umat[i, 0, i, 1] = Umat[i, 1, i, 0] = U

for i, j in product(range(norb), repeat=2):
    if i != j:
        Umat[i, 0, j, 1] = Umat[i, 1, j, 0] = Up
        Umat[i, 0, j, 0] = Umat[i, 1, j, 1] = Up - J

Umat = Umat.reshape((2*norb, 2*norb))

np.savetxt("u.in", Umat, fmt="%.4f")
print("'u.in' generated")
