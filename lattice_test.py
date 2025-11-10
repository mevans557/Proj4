import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import networkx as nx


offdi = sp.linalg.toeplitz([0, 1, 0, 0, 0, 0])
I = np.eye(6)

A = np.kron(offdi, I) + np.kron(I, offdi)
# plt.matshow(A)
# plt.show()

rows, cols = np.where(A == 1)
edges = zip(rows.tolist(), cols.tolist())
gr = nx.Graph()
gr.add_edges_from(edges)
nx.draw(gr, node_size = 10)
plt.show()