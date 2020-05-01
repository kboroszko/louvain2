import networkx as nx
import matplotlib.pyplot as plt

G = nx.Graph()

file = "data/output.dat"
weighted = False
cliques = []
with open(file) as fp:
    line = fp.readline()
    while line and line[0] == "%":
        line = fp.readline()
    [num_rows, num_cols, num_lines] = [int(float(x)) for x in line.split()]
    G.add_nodes_from(range(num_cols))
    for i in range(num_lines):
        line = fp.readline()
        arr = [int(float(x))-1 for x in line.split()]
        # print("parsing line: ", arr)
        if not weighted :
            G.add_edge(arr[0], arr[1])
            G.add_edge(arr[1], arr[0])
    for i in range(num_cols):
        cliques.append(int(fp.readline()))

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n+1)

cmap = get_cmap(max(cliques))
colors = [cmap(x+1) for x in cliques]

nx.draw_networkx(G, node_shape='o', node_color=colors, with_labels=True)

plt.show()