import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import networkx as nx


SIZE = 20  # Side length of lattice
prob = 0.5  # Probability of edge existence

state_dict = {0: 'H', 1: 'S', 2: 'I'}
colour_dict = {'H': "black", 'S': "green", 'I': "red"}


def draw_graph(M, sts):
    '''
    Draws the graph given in an adjacency matrix
    Params:
    M: An adjacency matrix of a graph
    sts: A dictionary of nodes with states of the graph
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


# Create a simple lattice graph
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
rowsout, colsout, vals = sp.sparse.find(sp.sparse.triu(A))  # Only upper tri
edgesout = [edge for edge in zip(rowsout.tolist(), colsout.tolist())
            if prob < np.random.rand()]
for edge in edgesout:
    A[edge[0], edge[1]] = 0  # FIXTHIS
    A[edge[1], edge[0]] = 0  # Also set lower triangle entry to 0
A.eliminate_zeros()


# Initialise SI vectors
suscepts = np.ones(SIZE**2)
infects = np.zeros(SIZE**2)


# Begin with an infected
infects[SIZE**2//2 + SIZE//2] = 1
suscepts[SIZE**2//2 + SIZE//2] = 0


# Naive time-step iteration
while 1 in suscepts:
    # Ensure exit after no change
    oldinfects = np.copy(infects)
    # Set up state vectors to display correctly
    statevec = suscepts + 2*infects
    states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
    # Infect neighbours
    infects = infects + ((A @ infects) * suscepts)
    # These style of expression reset all entries to 1 or 0
    infects[infects > 1] = 1
    infects[infects < 0] = 0
    # Make sure vectors sum to 1
    suscepts = suscepts - infects
    suscepts[suscepts < 0] = 0
    suscepts[suscepts > 1] = 1
    # Draw the graph
    draw_graph(A, states)
    # break loop if no updates occurred
    if (np.array_equal(infects, oldinfects)):
        break
statevec = suscepts + 2*infects
states = {x: state_dict[statevec[x]] for x in range(SIZE**2)}
draw_graph(A, states)
