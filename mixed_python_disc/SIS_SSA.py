import numpy as np
from matplotlib import pyplot as plt
from multiprocessing import Process, Pipe, current_process
from griddify import griddify


def iterate(t, I, N, kI, kR, alpha, END, conn):
    Is = np.array([])
    ts = np.array([])
    while ((t < END) and (alpha != 0)):
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
    conn.send((ts, Is))


t0 = 0
I0 = 1
N0 = 500
kI0 = 0.002
kR0 = 0.1
alpha0 = kR0*I0 + kI0*I0*(N0-I0)
iters = 20
END = 30
h = 0.01

tcollect = []
Icollect = []
Scollect = []
conns = []
processes = []

for i in range(iters):
    par_conn, child_conn = Pipe()
    process = Process(target=iterate, args=(
        t0, I0, N0, kI0, kR0, alpha0, END, child_conn,))
    processes.append(process)
    conns.append(par_conn)
    process.start()

for conn in conns:
    ts, Is = conn.recv()
    tcollect.append(ts)
    Icollect.append(Is)

for proc in processes:
    proc.join()

# for i in range(iters):
    # plt.plot(tcollect[i], Icollect[i])
    # plt.plot(tcollect[i], Scollect[i])


# Deterministic Mass-Action
td = np.arange(0, END, h)
Id = (N0-kR0/kI0)/(1+(N0-kR0/kI0-1)*np.exp(kR0*td-N0*kI0*td))


# RK4 on Moment closure bit
Ms = []
Vs = []
M = I0
V = 0


def f(M, V):
    # ODE for M
    return kI0*(N0 - M)*M - kR0*M - kI0*V


def g(M, V):
    # ODE for V
    return (kR0 + kI0*N0)*M + (2*kI0*N0 - 2*kR0 - kI0)*(V + M**2) - 6*kI0*M*V \
        - 2*kI0*(M**3)


# Actual RK4
# for t in td:
#     Ms.append(M)
#     Vs.append(V)
#     k1m = f(M, V)
#     k1v = g(M, V)
#     k2m = f((M + h*k1m/2), (V + h*k1v/2))
#     k2v = g((M + h*k1m/2), (V + h*k1v/2))
#     k3m = f((M + h*k2m/2), (V + h*k2v/2))
#     k3v = g((M + h*k2m/2), (V + h*k2v/2))
#     k4m = f((M + h*k3m), (V + h*k3v))
#     k4v = g((M + h*k3m), (V + h*k3v))
#     M = M + h*(k1m)  # + 2*k2m + 2*k3m + k4m)
#     V = V + h*(k1v)  # + 2*k2v + 2*k3v + k4v)
#     print(M, V)


# GRID
tgrid = np.arange(0, END, END/(N0*4))
Igrid = []
for i in range(iters):
    Igrid.append(griddify(np.array(tcollect[i]), tgrid, np.array(Icollect[i])))

Igrida = np.array(Igrid)
avs = np.sum(Igrida, 0)/iters
sds = np.sqrt(np.sum(Igrida**2, 0)/iters - avs**2)

plt.plot(tgrid, avs, "g")
plt.plot(tgrid, avs + sds, "--r")
plt.plot(tgrid, avs - sds, "--r")
plt.plot(td, Id, "--b")
# plt.plot(td, Ms, "--b")
# plt.fill_between(td, Ms - np.sqrt(Vs), Ms + np.sqrt(Vs), color=(0, 0, 0, 0.2))
plt.fill_between(tgrid, avs - sds, avs + sds, color=(1, 0, 0, 0.2))


plt.show()
