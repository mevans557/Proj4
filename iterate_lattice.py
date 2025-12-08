import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import networkx as nx

from graph_handling import make_perc_graph
from single_iterates import Gillespie_iterateSIS

SIZE = 20  # Side length of lattice
prob = 0.3  # Probability of edge existence
iters = 1   # Number of iterations to carry out
kInf = 0.5  # Chance of an infection along an edge in some time unit
KUninf = 0.3  # Chance of an infected becoming susceptible again in some time unit


# Track all data over iterations
tcollect = []
Icollect = []
Scollect = []


# Iterate
for i in range(iters):
    # Initialise SI vectors
    suscepts = np.ones(SIZE**2)
    infects = np.zeros(SIZE**2)
    # Begin with an infected
    infects[SIZE**2//2 + SIZE//2] = 1
    suscepts[SIZE**2//2 + SIZE//2] = 0
    # Make a lattice
    newLattice = make_perc_graph(SIZE, prob)
    tmpt, tmpI, tmpS = Gillespie_iterateSIS(newLattice, suscepts, infects,
                                            kInf, KUninf)
    tcollect.append(tmpt)
    Icollect.append(tmpI)
    print(i)
    # Scollect.append(tmpS)


# Plot data
for i in range(iters):
    plt.plot(tcollect[i], Icollect[i])
    # plt.plot(tcollect[i], Scollect[i])

# tlog = np.arange(0, 120, 0.1)
# Ilog = 9890/(1+np.e**(-(tlog-50)/12))
# plt.plot(tlog, Ilog, label="logistic")

plt.title("Graph SIS model Infecteds over Time")
plt.xlabel("Time t")
plt.ylabel("Infecteds I(t)")
plt.legend()
plt.show()
