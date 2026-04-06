import numpy as np
import scipy as sp
import networkx as nx
import pandas as pd
from matplotlib import pyplot as plt
from multiprocessing import Process, Pipe
from griddify import griddify


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
    gen = np.random.default_rng()
    r_1 = gen.uniform()
    r_2 = gen.uniform()
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


def iterate_Gillespie(G, I, kI, kS, stop=100, step=0.1, drawer=None):
    '''
    Iterates individual SIS on a graph using the Gillespie Algorithm
    Params:
    G: A networkx graph
    I: State vector of initial infected populations on a graph
    kI: infection rate
    delta: recovery rate (no immunity)
    stop: optional, time to stop simulating at
    step: optional, time steps to take (deprecated)
    draw: optional, whether the graph should be drawn (if so, drawn every 50 loops)
    Returns two lists, one storing times, the other storing infected populations
    '''
    t = 0
    ts = [t]
    Is = [I]
    Isum = []
    Isum.append(np.sum(I))
    loop = 0
    A = nx.to_scipy_sparse_array(G)
    while (t <= stop):
        loop += 1
        if ((drawer != None) and (loop % 50 == 0)):
            drawer(G, I)
        I, tau = evolve_graph(A, I, kI, kS)
        t += tau
        ts.append(t)
        Is.append(I)  # Something buggy happens with Is
        Isum.append(np.sum(I))
    return ts, Isum


def mult_Gillespie(G, I, kI, kS, connect, stop=100, step=0.1, drawer=None):
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
    ts, Isum = iterate_Gillespie(G, I, kI, kS, stop, step, drawer=drawer)
    connect.send((ts, Isum))


def draw_graph(G, I, colours=plt.cm.viridis):
    '''
    Draws a given graph, with colouring based on infection
    Params:
    G: A networkx graph
    I: State vector of infected populations on graph
    colours (optional): the colourmap to use
    '''
    nx.set_node_attributes(G, "green", name="node_color")
    nx.draw(G, node_color=I, node_size=100,
            vmin=0, vmax=1, cmap=colours)
    plt.show()


def make_perc_graph(width, length, pr, seed):
    '''
    Creates a networkx lattice graph, on which percolation has been performed.
    Params:
    width: the width of the lattice
    length: the length of the lattice
    pr: the probability of a given edge existing
    seed: seed for random generator, for repeatability
    Returns a networkx graph
    '''
    G = nx.grid_2d_graph(width, length)
    rng = np.random.default_rng(seed)
    remove = [edge for edge in G.edges if (rng.random() > pr)]
    G.remove_edges_from(remove)
    return G


def make_fb_graph():
    facebook = pd.read_csv(
        # Dataset from the SNAP database
        "https://snap.stanford.edu/data/facebook_combined.txt.gz",
        compression="gzip",
        sep=" ",
        names=["start_node", "end_node"],
    )
    G = nx.from_pandas_edgelist(facebook, "start_node", "end_node")
    return G


def draw_perc(gr, I, colours=plt.cm.viridis):
    '''
    Draws a percolation graph, with correct positioning
    Params:
    gr: A networkx lattice graph, on which percolation has been performed.
    '''
    pos = nx.circular_layout(gr)
    for k in pos:
        pos[k] = np.array(k)
    plt.figure(figsize=(7, 7))
    nx.draw(gr, pos=pos, node_color=I, node_size=20,
            vmin=0, vmax=1, cmap=colours)
    plt.show()

if __name__ == '__main__':

    # PERCOLATION STUFF
    # WIDTH = 20
    # LENGTH = 25
    # mid = WIDTH*LENGTH//2 + WIDTH//2
    # G = make_perc_graph(20, 25, 0.55, 20)  # Seeded for repeatability
    # I = np.zeros(500)
    # I[mid] = 1


    # FACEBOOK STUFF
    I = np.zeros(4039)
    G = make_fb_graph()
    for i in range(20):
        I[i] = 1


    # Networkx Builtin Graphs
    # G = nx.watts_strogatz_graph(500, 40, 0.5)
    # I = np.zeros(500)
    # I[0] = 1

    k_I = 0.1
    k_S = 10
    iters = 20
    END = 2

    tcollect = []
    Icollect = []
    conns = []
    processes = []

    for i in range(iters):
        par_conn, child_conn = Pipe()
        process = Process(target=mult_Gillespie, args=(
            G, I, k_I, k_S, child_conn, END, 0.1, None))
        processes.append(process)
        conns.append(par_conn)
        process.start()

    for conn in conns:
        ts, Is = conn.recv()
        tcollect.append(ts)
        Icollect.append(Is)

    for proc in processes:
        proc.join()


    # GRID
    tgrid = np.arange(0, END, END/(I.size*4))
    Igrid = []
    for i in range(iters):
        Igrid.append(griddify(np.array(tcollect[i]), tgrid, np.array(Icollect[i])))

    Igrida = np.array(Igrid)
    avs = np.sum(Igrida, 0)/iters
    sds = np.sqrt(np.sum(Igrida**2, 0)/iters - avs**2)

    plt.plot(tgrid, avs, "g")
    plt.plot(tgrid, avs + sds, "--r")
    plt.plot(tgrid, avs - sds, "--r")
    plt.fill_between(tgrid, avs - sds, avs + sds, color=(1, 0, 0, 0.2))


    # for i in range(iters):
    #     plt.plot(tcollect[i], Icollect[i])


    # plt.title("SSA of SIS Individuals model on Graph Infecteds over Time")
    # plt.xlabel("Time t")
    # plt.ylabel("Infecteds I(t)")
    plt.show()
