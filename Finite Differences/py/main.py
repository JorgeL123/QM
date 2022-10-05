# libs
import numpy as np
import matplotlib.pyplot as plt

# constants
dx = 0.05
dt = 0.01
hbar = 1
m = 125
a = 10

# time and space
x = np.arange(-2 * a, 2 * a, dx)
t = np.arange(0, 300, dt)
c = 10
# wave psi
P = (2 / np.pi)**(1 / 4) * np.exp(-(x + a / 2)**2) * np.exp(1j * c * x)
T = np.arange(0, 500, 1)
n = 0
V0 = 0.3
r = hbar * c / (2 * m) * dt / dx

# define potential
V = np.zeros_like(x)

V[np.abs(x) < a / 5] = V0
# V[np.logical_and(np.abs(x) <= a / 5, np.abs(x) >= a / 15)] = V0
# V = x**2 / 40 - 0.5
# V = (x**4 / 25 - x**2 + 4) / 15 - 0.3

# ---

for k in range(len(t)):
    if t[k] == T[n]:
        plt.clf()
        plt.plot(x, V, 'k', linewidth=1.1)
        plt.plot(x, np.abs(P), 'b', linewidth=1.5)
        plt.plot(x, np.real(P), 'r', linewidth=1.5)
        plt.axis([-a, a, -1.1, 1.1])
        plt.xlabel("x")
        plt.ylabel('P')
        plt.legend(['V', f'|P(x,{t[k]})|', f'Re P(x,{t[k]})'])
        plt.show()
        n += 1
        break
