import numpy as np


def SchroE(P, V, dx, hbar, m):
    dP = np.zeros_like(P)

    dP[1] = 1j * hbar / (2 * m) * (P[0] - 2 * P[1] +
                                   P[3]) / (dx**2) - 1j * V[1] * P[1] / hbar
    dP[-2] = 1j * hbar / (2 * m) * (P[-1] - 2 * P[-2] + P[-3]) / (
        dx**2) - 1j * V[-2] * P[-2] / hbar

    for i in range(2, len(P) - 2):
        dP[i] = 1j * hbar / (
            2 * m) * (-P[i + 2] + 16 * P[i + 1] - 30 * P[i] + 16 * P[i - 1] -
                      P[i - 2]) / (12 * dx**2) - 1j * V[i] * P[i] / hbar

    return dP
