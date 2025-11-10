import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import networkx as nx


SIZE = 100  # Side length of lattice


def draw_graph(M):
    '''
    Draws the graph given in an adjacency matrix
    Params: M, an adjacency matrix of a graph
    Note: this will not necessarily draw a lattice graph "sensibly"
    '''
    rows, cols, vals = sp.sparse.find(M)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    offset = np.array([0, 0])
    pos = nx.circular_layout(gr)
    for i in range(SIZE**2):
        pos[i] = offset + np.array([i % SIZE, (i//SIZE)])
    nx.draw(gr, pos, node_size=20)
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
prob = 0.5
rowsout, colsout, vals = sp.sparse.find(sp.sparse.triu(A))  # Only upper tri
edgesout = [edge for edge in zip(rowsout.tolist(), colsout.tolist())
            if prob < np.random.rand()]
for edge in edgesout:
    A[edge[0], edge[1]] = 0  # FIXTHIS
    A[edge[1], edge[0]] = 0  # Also set lower triangle entry to 0
A.eliminate_zeros()


draw_graph(A)
