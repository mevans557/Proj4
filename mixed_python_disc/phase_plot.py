from matplotlib import pyplot as plt
import numpy as np

k_I = 0.002
N = 500
k_S = 0.1

I = np.linspace(0, 500, 1000)

d = k_I * (N - I) * I - k_S * I
plt.plot(I, np.zeros(1000), 'k')
plt.plot(I, d, 'g')
plt.scatter(0, 0, facecolors='none', edgecolors='green')
plt.plot((N - (k_S/k_I)), 0, 'go')
plt.show()
