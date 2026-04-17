import networkx as nx
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from SIS_individuals import make_perc_graph, make_fb_graph


G = make_fb_graph()
degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
plt.bar(*np.unique(degree_sequence, return_counts=True))
plt.show()