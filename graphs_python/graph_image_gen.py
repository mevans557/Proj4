import networkx as nx
from matplotlib import pyplot as plt


G = nx.star_graph(4)
plt.figure(figsize=(5, 5))
nx.draw(G)
plt.show()
