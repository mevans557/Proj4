import networkx as nx
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt


def evolve(A, I, N, kI, delta):
    '''
    Runs discrete metapop SIS on a graph for use with a naïve time-step algorithm
    Params:
    A: An adjacency matrix of a graph
    I: State vector storing infected populations on a graph
    N: vector of node populations
    kI: infection rate
    delta: recovery rate (no immunity)
    Returns a state vector storing new infected populations on the graph
    '''
    newI = I + (np.round(kI * (N-I) * (A @ I) - delta*I))
    for i in range(I.size):
        if (newI[i] > N[i]):
            newI[i] = N[i]
    return newI


def Gillespie(A, I, N, kI, delta):
    '''
    Runs discrete metapop SIS on a graph for use with Gillespie iterator
    Params:
    A: An adjacency matrix of a graph
    I: State vector storing infected populations on a graph
    N: vector of node populations
    kI: infection rate
    delta: recovery rate (no immunity)
    Returns a tuple of a state vector storing new infected populations on the graph, and a time elapsed for the reaction
    '''
    alpha_inf_in = kI * I * \
        (N-I)  # ith element is propensity for ith node infects++ internally
    alpha_rec = delta * I  # ith element is propensity for ith node infects-- internally
    # ith element is propensity for ith node infects++ externally
    alpha_inf_bet = kI * (A @ I) * (N-I)
    alpha = sum(alpha_inf_bet + alpha_inf_in + alpha_rec)  # Overall propensity
    r1 = np.random.random()
    tau = 1/alpha * np.log(1/r1)
    r2 = np.random.random()
    inf_cum = np.cumulative_sum((alpha_inf_in + alpha_inf_bet)/alpha)
    rec_cum = np.cumulative_sum(alpha_rec/alpha)
    Inew = I
    # Evolve in case first infection reaction
    if (r2 < inf_cum[0]):
        Inew[0] = Inew[0] + 1
        return Inew, tau
    # Evolve in case one of the other infection reactions
    for i in range(np.size(I)-2):
        if ((r2 >= inf_cum[i]) and (r2 < inf_cum[i+1])):
            Inew[i+1] = Inew[i+1] + 1
            return Inew, tau
    # Evolve in case final infection reaction
    if ((r2 >= inf_cum[-2]) and (r2 < inf_cum[-1])):
        Inew[-1] = Inew[-1] + 1
        return Inew, tau
    # Evolve in case first recovery reaction
    if ((r2 >= inf_cum[-1]) and (r2 < inf_cum[-1] + rec_cum[0])):
        Inew[0] = Inew[0] - 1
        return Inew, tau
    # Evolve in case one of the other recovery reactions
    for i in range(np.size(I)-1):
        if ((r2 >= inf_cum[-1] + rec_cum[i]) and (r2 < inf_cum[-1] + rec_cum[i+1])):
            Inew[i+1] = Inew[i+1] - 1
            return Inew, tau
    # Shouldn't get here
    print("oops, shouldn't be here")


def iterate(G, I, N, kI, delta, stop=100, step=0.1, draw=False):
    '''
    Iterates discrete metapop SIS on a graph using a naïve time-step algorithm
    NOTE: This is quite dodgy and behaves weirdly, particularly on dense
    graphs
    Params:
    G: A networkx graph
    I: State vector of initial infected populations on a graph
    N: vector of node populations
    kI: infection rate
    delta: recovery rate (no immunity)
    stop: optional, time to stop simulating at
    step: optional, time steps to take
    Returns two lists, one storing times, the other storing state vectors of infected population
    '''
    t = 0
    ts = [t]
    Is = [I]
    A = nx.to_scipy_sparse_array(G)
    while (t <= stop):
        if draw:
            draw_graph(G, I, N)
            print(I)
        I = evolve(A, I, N, kI, delta)
        t += step
        ts.append(t)
        Is.append(I)
    return ts, Is


def iterate_Gillespie(G, I, N, kI, delta, stop=100, step=0.1, draw=False):
    '''
    Iterates discrete metapop SIS on a graph using the Gillespie Algorithm
    Params:
    G: A networkx graph
    I: State vector of initial infected populations on a graph
    N: vector of node populations
    kI: infection rate
    delta: recovery rate (no immunity)
    stop: optional, time to stop simulating at
    step: optional, time steps to take
    draw: optional, whether the graph should be drawn (if so, drawn every 50 loops)
    Returns two lists, one storing times, the other storing state vectors of infected population
    '''
    t = 0
    ts = [t]
    Is = [I]
    Isum = [1]
    loop = 0
    A = nx.to_scipy_sparse_array(G)
    while (t <= stop):
        loop += 1
        if (draw and (loop % 50 == 0)):
            draw_graph(G, I, N)
        I, tau = Gillespie(A, I, N, kI, delta)
        t += tau
        ts.append(t)
        Is.append(I)  # Something buggy happens with Is
        Isum.append(np.sum(I))
    return ts, Isum


def draw_graph(G, I, N, colours=plt.cm.viridis):
    p = I/N
    nx.set_node_attributes(G, "green", name="node_color")
    nx.draw_circular(G, node_color=p, node_size=100,
                     vmin=0, vmax=1, cmap=plt.cm.plasma)
    plt.show()


G = nx.cycle_graph(5)
Ii = np.zeros(5)
Ni = np.array([100, 200, 300, 400, 500])
Ii[0] = 1


ts, Isum = iterate_Gillespie(G, Ii, Ni, 0.0002, 0.02, draw=False)

plt.plot(ts, Isum)
plt.show()
