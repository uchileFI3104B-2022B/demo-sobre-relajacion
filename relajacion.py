"""
Resuelve la ecuacion:
nabla^2 phi = -rho(x, y) = -2 (2 - x^2 - y^2)
C.B: phi(+/-1, y) = 0
     phi(x, +/-1) = 0

Utiliza el metodo de la sobre-relajacion sucesiva.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


def muestra_phi(phi, precision=3):
    with np.printoptions(precision=precision, suppress=True):
        print(phi[::-1, :])


def rho(i, j, h):
    x_i = -1 + h * i
    y_j = -1 + h * j
    return -2 * (2 - x_i**2 - y_j**2)


def convergio(phi, phi_anterior, rtol=0.1):
    not_zero = phi_anterior != 0
    dif_relativa = ((phi_anterior[not_zero] - phi[not_zero]) /
                    phi_anterior[not_zero])
    return np.fabs(dif_relativa).max() < rtol


def una_iteracion(phi, N_pasos, h, w=1):
    for i in range(1, N_pasos-1):
        for j in range(1, N_pasos-1):
            phi[i, j] = (1-w) * phi[i, j] + w/4 * (phi[i+1, j] + phi[i-1, j] +
                                                   phi[i, j+1] + phi[i, j-1] -
                                                   h**2 * rho(i, j, h))
    return phi


# Setup
Lx = Ly = 2
N_pasos = 31
h = Lx / (N_pasos - 1)
w = 1.8  # Parametro de sobre-relajacion

phi = np.zeros((N_pasos, N_pasos))

# Iteracion
phi = una_iteracion(phi, N_pasos, h, w=w)
phi_anterior = phi.copy()
phi = una_iteracion(phi, N_pasos, h, w=w)

counter = 2
while not convergio(phi, phi_anterior, rtol=1e-10) and counter < 2000:
    phi_anterior = phi.copy()
    phi = una_iteracion(phi, N_pasos, h, w=w)
    counter += 1

print("counter = {}".format(counter))
i_central = int((N_pasos-1) / 2)
print("phi(0, 0) = {}".format(phi[i_central, i_central]))

plt.figure(1)
plt.clf()

plt.imshow(phi, origin='bottom', extent=[-1, 1, -1, 1])
plt.contour(phi, origin='lower', extent=[-1, 1, -1, 1], colors='k')


fig = plt.figure(2)
fig.clf()
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(-1, 1, N_pasos)
y = np.linspace(-1, 1, N_pasos)
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, phi, rstride=1, cstride=1)

plt.show()
