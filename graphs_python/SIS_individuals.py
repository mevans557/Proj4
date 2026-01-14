import numpy as np
import scipy as sp
import networkx as nx
from matplotlib import pyplot as plt
from multiprocessing import Process, Pipe


def evolve_graph(A, I, k_I, k_S):
    '''
    Takes a single timestep evolution of the SIS model on a graph, where each
    node represents an individual
    Params:
    A: A scipy sparse array, the adjacency matrix of a graph
    I: A numpy array, containing 1 if the node is infected, 0 if not
    k_I: The infection rate
    k_S: The recovery rate (without immunity)
    Returns a tuple of time until reaction, and state after the reaction as tau, Inew
    '''
    r_1 = np.random.random()
    r_2 = np.random.random()
    alpha = sum(k_I * (1 - I) * (A @ I) + k_S * I)
    if (alpha == 0):
        return I, 5
    tau = 1/alpha * np.log(1/r_1)
    Isum = np.cumulative_sum(k_I * (1-I) * (A @ I) / alpha)
    Rsum = np.cumulative_sum(k_S * I/alpha)
    Inew = I
    if (r_2 < sum(k_I * (1-I) * (A @ I)) / alpha):
        for j in range(I.size):
            if (r_2 < Isum[j]):
                Inew[j] = 1
                return Inew, tau
    else:
        r_3 = r_2 - sum(k_I * (1-I) * (A @ I))/alpha
        for j in range(I.size):
            if (r_3 < Rsum[j]):
                Inew[j] = 0
                return Inew, tau
    # Here means r_2 = 1
    Inew[-1] = 0
    return Inew, tau


def iterate_Gillespie(G, I, kI, kS, stop=100, step=0.1, draw=False):
    '''
    Iterates individual SIS on a graph using the Gillespie Algorithm
    Params:
    G: A networkx graph
    I: State vector of initial infected populations on a graph
    kI: infection rate
    delta: recovery rate (no immunity)
    stop: optional, time to stop simulating at
    step: optional, time steps to take
    draw: optional, whether the graph should be drawn (if so, drawn every 50 loops)
    Returns two lists, one storing times, the other storing infected populations
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
            draw_graph(G, I)
        I, tau = evolve_graph(A, I, kI, kS)
        t += tau
        ts.append(t)
        Is.append(I)  # Something buggy happens with Is
        Isum.append(np.sum(I))
    return ts, Isum


def mult_Gillespie(G, I, kI, kS, connect, stop=100, step=0.1, draw=False):
    '''
    Iterates individual SIS on a graph using the Gillespie Algorithm and outputs to a pipe
    Params:
    G: A networkx graph
    I: State vector of initial infected populations on a graph
    kI: infection rate
    delta: recovery rate (no immunity)
    connect: the pipe to output to
    stop: optional, time to stop simulating at
    step: optional, time steps to take
    draw: optional, whether the graph should be drawn (if so, drawn every 50 loops)
    Sends a tuple to the pipe, of two lists, the first storing times, the other infected populations
    '''
    ts, Isum = iterate_Gillespie(G, I, kI, kS, stop, step, draw=draw)
    connect.send((ts, Isum))


def draw_graph(G, I, colours=plt.cm.viridis):
    nx.set_node_attributes(G, "green", name="node_color")
    nx.draw_circular(G, node_color=I, node_size=100,
                     vmin=0, vmax=1, cmap=plt.cm.plasma)
    plt.show()


G = nx.barabasi_albert_graph(500, 7)
I = np.zeros(500)
I[0] = 1
k_I = 0.02
k_S = 0.1

iters = 20

tcollect = []
Icollect = []
conns = []
processes = []

for i in range(iters):
    par_conn, child_conn = Pipe()
    process = Process(target=mult_Gillespie, args=(
        G, I, k_I, k_S, child_conn, 50, 0.1, False))
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


plt.title("SSA of SIS metapopulation model on Graph Infecteds over Time")
plt.xlabel("Time t")
plt.ylabel("Infecteds I(t)")
plt.show()
