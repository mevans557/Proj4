import numpy as np
from matplotlib import pyplot as plt
from multiprocessing import Process, Pipe, current_process


def iterate(t, I, N, kI, kR, alpha, conn):
    Is = np.array([])
    ts = np.array([])
    while ((t < 20) and (alpha != 0)):
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
N0 = 100
kI0 = 0.002
kR0 = 0.1
alpha0 = kR0*I0 + kI0*I0*(N0-I0)
iters = 20


tcollect = []
Icollect = []
Scollect = []
conns = []
processes = []

for i in range(iters):
    par_conn, child_conn = Pipe()
    process = Process(target=iterate, args=(
        t0, I0, N0, kI0, kR0, alpha0, child_conn,))
    processes.append(process)
    conns.append(par_conn)
    process.start()

for conn in conns:
    ts, Is = conn.recv()
    tcollect.append(ts)
    Icollect.append(Is)

for proc in processes:
    proc.join()

for i in range(iters):
    plt.plot(tcollect[i], Icollect[i])
    # plt.plot(tcollect[i], Scollect[i])


td = np.arange(0, 30, 0.1)
Id = (N0-kR0/kI0)/(1+(N0-kR0/kI0-1)*np.exp(kR0*td-N0*kI0*td))

plt.plot(td, Id, linestyle='dashed', color='black', label="Deterministic")


plt.title("SSA of SIS model Infecteds over Time")
plt.xlabel("Time t")
plt.ylabel("Infecteds I(t)")
plt.legend()
plt.show()
