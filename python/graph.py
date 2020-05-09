from IPython.display import SVG, display
import numpy as np
from sknetwork.data import karate_club, painters, movie_actor
from sknetwork.clustering import Louvain, BiLouvain, modularity, bimodularity
from sknetwork.linalg import normalize
from sknetwork.utils import bipartite2undirected, membership_matrix
from sknetwork.visualization import svg_graph, svg_digraph, svg_bigraph
from matplotlib import pyplot as plt


graph = karate_club(metadata=True)

adjacency = graph.adjacency
position = graph.position

svg_gr

print(adjacency)
print(position)

image = svg_graph(adjacency, position)
SVG(image)