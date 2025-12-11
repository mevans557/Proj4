import numpy as np
from matplotlib import pyplot as plt


def iterate(t, I, N, kI, kR, alpha):
    Is = np.array([])
    ts = np.array([])
    while ((t < 100) and (alpha != 0)):
        ts = np.append(ts, t)
        Is = np.append(Is, I)
        r1 = np.random.rand()
        r2 = np.random.rand()
        tau = (1/alpha)*np.log(1/r1)
        if (kR*I/alpha > r2):
            I -= 1
        else:
            I += 1
        t = t+tau
        alpha = kR*I + kI*I*(N-I)
    ts = np.append(ts, t)
    Is = np.append(Is, I)
    return ts, Is


t0 = 0
I0 = 1
N0 = 10000
kI0 = 0.00002
kR0 = 0.05
alpha0 = kR0*I0 + kI0*I0*(N0-I0)
iters = 10


tcollect = []
Icollect = []
Scollect = []

for i in range(iters):
    tmpt, tmpI = iterate(t0, I0, N0, kI0, kR0, alpha0)
    tcollect.append(tmpt)
    Icollect.append(tmpI)
    Scollect.append(N0-tmpI)


for i in range(iters):
    plt.plot(tcollect[i], Icollect[i])
    # plt.plot(tcollect[i], Scollect[i])


td = np.arange(0, 10, 0.1)
Id = (N0-kR0/kI0)/(1+(N0-kR0/kI0-1)*np.exp(kR0*td-N0*kI0*td))

plt.plot(td, Id, linestyle='dashed', color='black', label="Deterministic")


plt.title("SSA of SIS model Infecteds over Time")
plt.xlabel("Time t")
plt.ylabel("Infecteds I(t)")
plt.legend()
plt.show()
