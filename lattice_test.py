import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import networkx as nx


from graph_handling import draw_lattice, make_perc_graph


SIZE = 100  # Side length of lattice
prob = 0.3  # Probability of edge existence
iters = 1   # Number of iterations to carry out
kInf = 0.5  # Chance of an infection along an edge in some time unit
KUninf = 0.3  # Chance of an infected becoming susceptible again in some time unit

# For passing to draw_lattice()
state_dict = {0: 'H', 1: 'S', 2: 'I'}


def naive_iterateSI(M, Si, Ii):
    '''
    Naively time-step iterates through graph-based SI
    Params:
    M: An adjacency matrix of a graph
    Si: Initial susceptible population vector
    Ii: Initial infected population vector
    Returns a tuple of time, infected, susceptible arrays
    '''
    # Track populations over time
    Is = np.array([np.sum(Ii)])
    Ss = np.array([np.sum(Si)])
    ts = np.array([0])
    while 1 in Si:
        # Ensure exit after no change
        oldIi = np.copy(Ii)
        # Set up state vectors to display correctly. NOTE: Only used for display
        statevec = Si + 2*Ii
        states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
        # Infect neighbours
        Ii = Ii + ((M @ Ii) * Si)
        # These style of expression reset all entries to 1 or 0
        Ii[Ii > 1] = 1
        Ii[Ii < 0] = 0
        # Make sure vectors sum to 1
        Si = Si - Ii
        Si[Si < 0] = 0
        Si[Si > 1] = 1
        # Track populations
        ts = np.append(ts, ts[-1]+1)
        Is = np.append(Is, np.sum(Ii))
        Ss = np.append(Ss, np.sum(Si))
        # Draw the graph
        # draw_lattice(M, states)
        # break loop if no updates occurred
        if (np.array_equal(Ii, oldIi)):
            break
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
    # draw_lattice(M, states)
    # Return the tracked data
    return ts, Is, Ss


def Gillespie_iterateSI(M, Si, Ii, kI):
    '''
    Uses Gillespie-style algorithm to iterate through graph-based SI
    Params:
    M: An adjacency matrix of a graph
    Si: Initial susceptible population vector
    Ii: Initial infected population vector
    kI: Probability of an infection on an edge in some time unit
    Returns a tuple of time, infected, susceptible arrays
    '''
    # Draw initial graph
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
    # draw_lattice(M, states)
    # Track populations over time
    Is = np.array([np.sum(Ii)])
    Ss = np.array([np.sum(Si)])
    ts = np.array([0])
    while 1 in Si:
        # Ensure exit after no change
        oldIi = np.copy(Ii)
        # Set up state vectors to display correctly. NOTE: Only used for display
        statevec = Si + 2*Ii
        states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
        # Get infectable nodes
        Iable = (M @ Ii) * Si
        # Calculate next reaction time
        if (np.sum(Iable) == 0):
            break
        alpha = kI * np.sum(Iable)
        r1 = np.random.rand()
        tau = 1/alpha * np.log(1/r1)
        # Select a reaction to occur
        change = np.random.choice(np.nonzero(Iable)[0])
        Ii[change] = 1
        Si[change] = 0
        # Track populations
        ts = np.append(ts, ts[-1] + tau)
        Is = np.append(Is, np.sum(Ii))
        Ss = np.append(Ss, np.sum(Si))
        # Draw the graph
        draw_lattice(M, states)
        # break loop if no updates occurred
        if (np.array_equal(Ii, oldIi)):
            break
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
    draw_lattice(M, states)
    # Return the tracked data
    return ts, Is, Ss


def Gillespie_iterateSIS(M, Si, Ii, kI, kS):
    '''
    Uses Gillespie-style algorithm to iterate through graph-based SI
    Params:
    M: An adjacency matrix of a graph
    Si: Initial susceptible population vector
    Ii: Initial infected population vector
    kI: Probability of an infection on an edge in some time unit
    kS: Probability of a node going back from I -> S in some time unit
    Returns a tuple of time, infected, susceptible arrays
    '''
    # Draw initial graph
    # statevec = Si + 2*Ii
    # states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
    # draw_lattice(M, states)
    # Track populations over time
    Is = np.array([np.sum(Ii)])
    Ss = np.array([np.sum(Si)])
    ts = np.array([0])
    while (ts[-1] < 100):
        if (np.count_nonzero(Ii) == 0):
            break
        # Set up state vectors to display correctly. NOTE: Only used for display
        # statevec = Si + 2*Ii
        # states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
        # Get infectable nodes
        Iable = (M @ Ii) * Si
        # Calculate next reaction time
        alpha = (kI * np.sum(Iable)) + (kS * np.sum(Ii))
        r1 = np.random.rand()
        tau = 1/alpha * np.log(1/r1)
        # Select a reaction to occur
        r2 = np.random.rand()
        if (r2 < (1/alpha * kI*np.sum(Iable))):
            change = np.random.choice(np.nonzero(Iable)[0])
            Ii[change] = 1
            Si[change] = 0
        else:
            change = np.random.choice(np.nonzero(Ii)[0])
            Ii[change] = 0
            Si[change] = 1
        # Track populations
        ts = np.append(ts, ts[-1] + tau)
        Is = np.append(Is, np.sum(Ii))
        Ss = np.append(Ss, np.sum(Si))
        # Draw the graph
        # draw_lattice(M, states)
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
    draw_lattice(M, states)
    # Return the tracked data
    return ts, Is, Ss


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
    newLattice = make_perc_graph()
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
