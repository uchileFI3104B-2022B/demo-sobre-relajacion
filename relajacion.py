"""
Resuelve la ecuacion:
nabla^2 phi = -rho(x, y) = -2 (2 - x^2 - y^2)
C.B: phi(+/-1, y) = 0
     phi(x, +/-1) = 0

Utiliza el metodo de la sobre-relajacion sucesiva.
"""

import numpy as np


def rho(i, j, h):
    x_i = -1 + h * i
    y_j = -1 + h * j
    return -2 * (2 - x_i**2 - y_j**2)


# Setup
Lx = Ly = 2
N_pasos = 5
h = Lx / (N_pasos - 1)
w = 1

phi = np.zeros((N_pasos, N_pasos))

for i in range(1, N_pasos-1):
    for j in range(1, N_pasos-1):
        phi[i, j] = (1-w) * phi[i, j] + w / 4 * (phi[i+1, j] + phi[i-1, j] +
                                                 phi[i, j+1] + phi[i, j-1] -
                                                 h**2 * rho(i, j, h))