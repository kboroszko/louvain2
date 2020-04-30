import networkx as nx
import matplotlib.pyplot as plt

G = nx.Graph()

file = "data/mycielskian4_bidirectional.mtx"
weighted = False

with open(file) as fp:
    line = fp.readline()
    while line and line[0] == "%" :
        line = fp.readline()
    [num_rows, num_cols, num_lines] = [int(float(x)) for x in line.split()]
    G.add_nodes_from(range(num_cols))
    for i in range(num_lines):
        line = fp.readline()
        arr = [int(float(x))-1 for x in line.split()]
        # print("parsing line: ", arr)
        if not weighted :
            G.add_edge(arr[0], arr[1])
            # G.add_edge(arr[1], arr[0])

def assign_color(x):
    reds = [0,2,4,6]
    greens = [3, 5, 9, 7, 1]
    if x in reds :
        return 'r'
    if x in greens :
        return 'g'
    return 'm'


colors = [assign_color(x) for x in range(G.number_of_nodes())]
# labels = {1:’Start’, 2:’2’, 3:’3’, 4:’4’,
#    5:’5’, 6:’6’, 7:’7’, 8:’End’}
# sizes = [800, 300, 300, 300, 300, 600, 300, 800]
# nx.draw_networkx(G, node_color=colors, node_shape=‘D’,
#      with_labels=True, labels=labels,
#      node_size=sizes)

nx.draw_networkx(G, node_shape='o', node_color=colors, with_labels=True)

plt.show()