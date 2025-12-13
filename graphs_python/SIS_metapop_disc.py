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


def draw_graph(G, I, N, colours=plt.cm.viridis):
    p = I/N
    nx.set_node_attributes(G, "green", name="node_color")
    nx.draw_circular(G, node_color=p, node_size=100, cmap=plt.cm.rainbow)
    plt.show()


G = nx.complete_graph(5)
Ii = np.zeros(5)
Ni = np.array([100, 200, 300, 400, 500])
Ii[0] = 1


iterate(G, Ii, Ni, 0.002, 0.02, draw=True)
