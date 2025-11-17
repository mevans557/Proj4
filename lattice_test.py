import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import networkx as nx


SIZE = 100  # Side length of lattice
prob = 0.7  # Probability of edge existence
iters = 10  # Number of iterations to carry out

state_dict = {0: 'H', 1: 'S', 2: 'I'}
colour_dict = {'H': "black", 'S': "green", 'I': "red"}


def draw_graph(M, sts):
    '''
    Draws the graph given in an adjacency matrix
    Params:
    M: An adjacency matrix of a graph
    sts: A dictionary of nodes with states of the graph
    No return types
    '''
    rows, cols, vals = sp.sparse.find(M)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    toadd = [i for i in range(SIZE**2)]
    gr.add_nodes_from(toadd)
    gr.add_edges_from(edges)
    offset = np.array([0, 0])
    pos = nx.circular_layout(gr)
    nx.set_node_attributes(gr, sts, name="state")
    node_colour = [colour_dict[nx.get_node_attributes(gr, "state")[i]]
                   for i in range(SIZE**2)]
    for i in range(SIZE**2):
        pos[i] = offset + np.array([i % SIZE, (i//SIZE)])
    plt.figure(figsize=(7, 7))
    nx.draw(gr, pos, node_size=20, node_color=node_colour)
    plt.show()


def make_perc_graph():
    '''
    Create a simple lattice graph
    Returns a percolation-type lattice graph, using global prob
    '''
    base_row = [0, 1]
    for i in range(2, SIZE):
        base_row.append(0)
    offdi = sp.linalg.toeplitz(base_row)
    I = sp.sparse.eye_array(SIZE)
    A = sp.sparse.kron(offdi, I, format='csr') \
        + sp.sparse.kron(I, offdi, format='csr')
    # plt.matshow(A)
    # plt.show()

    # Randomly remove edges
    rowsout, colsout, vals = sp.sparse.find(
        sp.sparse.triu(A))  # Only upper tri
    edgesout = [edge for edge in zip(rowsout.tolist(), colsout.tolist())
                if prob < np.random.rand()]
    for edge in edgesout:
        A[edge[0], edge[1]] = 0  # FIXTHIS
        A[edge[1], edge[0]] = 0  # Also set lower triangle entry to 0
    A.eliminate_zeros()
    return A


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
        # Set up state vectors to display correctly
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
        # draw_graph(M, states)
        # break loop if no updates occurred
        if (np.array_equal(Ii, oldIi)):
            break
    statevec = Si + 2*Ii
    states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
    # draw_graph(M, states)
    # Return the tracked data
    return ts, Is, Ss


# Initialise SI vectors
suscepts = np.ones(SIZE**2)
infects = np.zeros(SIZE**2)
# Begin with an infected
infects[SIZE**2//2 + SIZE//2] = 1
suscepts[SIZE**2//2 + SIZE//2] = 0


# Track all data over iterations
tcollect = []
Icollect = []
Scollect = []

# Iterate
for i in range(iters):
    newLattice = make_perc_graph()
    tmpt, tmpI, tmpS = naive_iterateSI(newLattice, suscepts, infects)
    tcollect.append(tmpt)
    Icollect.append(tmpI)
    Scollect.append(tmpS)

# Plot data
for i in range(iters):
    plt.plot(tcollect[i], Icollect[i])
    plt.plot(tcollect[i], Scollect[i])

tlog = np.arange(0, 120, 0.1)
Ilog = 9890/(1+np.e**(-(tlog-50)/12))
plt.plot(tlog, Ilog, label="logistic")

plt.title("Graph SI model Infecteds over Time")
plt.xlabel("Time t")
plt.ylabel("Infecteds I(t)")
plt.legend()
plt.show()
