import numpy as np
from matplotlib import pyplot as plt


def iterate(t, I, N, kI, alpha):
    Is = np.array([])
    ts = np.array([])
    while ((t < 100) and (alpha != 0)):
        ts = np.append(ts, t)
        Is = np.append(Is, I)
        r1 = np.random.rand()
        tau = (1/alpha)*np.log(1/r1)
        t = t+tau
        I += 1
        alpha = kI*I*(N-I)
    ts = np.append(ts, t)
    Is = np.append(Is, I)
    return ts, Is


t0 = 0
I0 = 1
N0 = 10000
kI0 = 0.00002
alpha0 = kI0*I0*(N0-I0)
iters = 10

tcollect = []
Icollect = []


for i in range(iters):
    tmpt, tmpI = iterate(t0, I0, N0, kI0, alpha0)
    tcollect.append(tmpt)
    Icollect.append(tmpI)


for i in range(iters):
    plt.plot(tcollect[i], Icollect[i])


# td = np.arange(0, 5, 0.1)
# Id = N0/(1+(N0-1)*np.exp(-N0*kI0*td))
# plt.plot(td, Id, linestyle='dashed', color='black', label="Deterministic")


plt.title("SSA of SI model Infecteds over Time")
plt.xlabel("Time t")
plt.ylabel("Infecteds I(t)")
plt.legend()
plt.show()
