import scipy as sp
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt


# For graph plotting purposes
colour_dict = {'H': "black", 'S': "green", 'I': "red"}


def draw_lattice(M, sts, side_length):
    '''
    Draws the graph given in an adjacency matrix
    Params:
    M: An adjacency matrix of a graph
    sts: A dictionary of nodes with states of the graph
    side_length: the side length of the lattice graph being draw
    No return types
    '''
    rows, cols, vals = sp.sparse.find(M)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    toadd = [i for i in range(side_length**2)]
    gr.add_nodes_from(toadd)
    gr.add_edges_from(edges)
    offset = np.array([0, 0])
    pos = nx.circular_layout(gr)
    nx.set_node_attributes(gr, sts, name="state")
    node_colour = [colour_dict[nx.get_node_attributes(gr, "state")[i]]
                   for i in range(side_length**2)]
    for i in range(side_length**2):
        pos[i] = offset + np.array([i % side_length, (i//side_length)])
    plt.figure(figsize=(7, 7))
    nx.draw(gr, pos, node_size=20, node_color=node_colour)
    plt.show()


def make_perc_graph(side_length, p_exist):
    '''
    Create a simple lattice graph
    Returns the adjacency matrix of a percolation-type square lattice graph
    Params:
    side_length: the side length of the square lattice to create
    p_exist: the probability of an edge to exist
    '''
    base_row = [0, 1]
    for i in range(2, side_length):
        base_row.append(0)
    offdi = sp.linalg.toeplitz(base_row)
    I = sp.sparse.eye_array(side_length)
    A = sp.sparse.kron(offdi, I, format='csr') \
        + sp.sparse.kron(I, offdi, format='csr')
    # plt.matshow(A)
    # plt.show()

    # Randomly remove edges
    rowsout, colsout, vals = sp.sparse.find(
        sp.sparse.triu(A))  # Only upper tri
    edgesout = [edge for edge in zip(rowsout.tolist(), colsout.tolist())
                if p_exist < np.random.rand()]
    for edge in edgesout:
        A[edge[0], edge[1]] = 0
        A[edge[1], edge[0]] = 0
    A.eliminate_zeros()
    return A
