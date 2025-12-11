import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import networkx as nx

from graph_handling import make_perc_graph
from single_discrete_lattice_iterates import Gillespie_iterateSIS

SIZE = 20  # Side length of lattice
prob = 0.3  # Probability of edge existence
iters = 1   # Number of iterations to carry out
kInf = 0.5  # Chance of an infection along an edge in some time unit
KUninf = 0.3  # Chance of an infected becoming susceptible again in some time unit


def get_data(iteror, lattice, sus, inf, kI, kU, tc, Ic, Sc=[]):
    '''
    Runs a single iteration of a simulation algorithm
    params:
    iteror: simulation algorithm single iteration function to use
    lattice: lattice to pass to iteror
    sus: initial susceptible array to pass to iteror
    inf: intiial infected array to pass to iteror
    kI: chance of infection to pass to iteror
    kU: chance of non-immune recovery to pass to iteror
    tc: list to collect times in
    Ic: list to collect infected pops in
    Sc: list to collect susceptible pops in
    no return, appends to tc, Ic, Sc
    '''
    tmpt, tmpI, tmpS = iteror(lattice, sus, inf, kI, kU)
    tc.append(tmpt)
    Ic.append(tmpI)
    Sc.append(tmpS)


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
    get_data(Gillespie_iterateSIS, newLattice, suscepts, infects, kInf, KUninf,
             tcollect, Icollect)

# Plot data
for i in range(iters):
    plt.plot(tcollect[i], Icollect[i])
    # plt.plot(tcollect[i], Scollect[i])

plt.title("Graph SIS model Infecteds over Time")
plt.xlabel("Time t")
plt.ylabel("Infecteds I(t)")
plt.legend()
plt.show()
