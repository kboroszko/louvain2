from pyvis.network import Network
from matplotlib import pyplot as plt
import matplotlib


file = "./../data/example.mtx"

def readData(filename):
    cliques = []
    srcs = []
    dests = []
    vals = []
    with open(file) as fp:
        line = fp.readline()
        while line and line[0] == "%":
            line = fp.readline()
        [num_rows, num_cols, num_lines] = [int(float(x)) for x in line.split()]
        for i in range(num_lines):
            line = fp.readline()
            arr = [float(x)-1 for x in line.split()]
            srcs.append(int(arr[0]))
            dests.append(int(arr[1]))
            if len(arr) > 2 :
                vals.append(arr[2])
            else :
                vals.append(1.)
        for i in range(num_cols):
            line = fp.readline()
            if line :
                cliques.append(int(line))
            else :
                break
    return srcs, dests, vals, cliques, num_cols

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.get_cmap(name, n+1)


def get_hex_color(n, cmap):
    return matplotlib.colors.rgb2hex(cmap(n))

def set_node_size(net, node_id, size):
    for node in net.nodes :
        if node["id"] == node_id :
            node["size"] = size
            break


srcs, dests, vals, cliques, n = readData(file)

print("edgelist=[", end='')
for src, dest, val in zip(srcs, dests, vals) :
    print("({0},{1},{2}), ".format(int(src), int(dest), float(val)), end='')
print(']')

edges = len(vals)
print("found ", n, "nodes and", edges, "edges")

net = Network(height='800px', width='1000px')

for i in range(n):
    net.add_node(str(i), str(i), title=str(i))

for src, dest, val in zip(srcs, dests, vals) :
    print("adding ", src, " -> ", dest, "    v=", val)
    if src == dest :
        # net.add_node(str(src), str(src), size=val)
        set_node_size(net, node_id=str(src), size=val)
    else :
    #     net.add_node(str(src), str(src), title=src)
    #     net.add_node(str(dest), str(dest), title=dest)
        net.add_edge(str(src), str(dest), value=val)

print("got {0} nodes".format(n))

if (len(cliques) > 1) :
    unique_cl = sorted(set(cliques))
    cmap = get_cmap(len(unique_cl))
else :
    unique_cl = range(n)
    cmap = get_cmap(n)
#
if len(cliques) > 1:
    for node in net.nodes:
        c = cliques[int(node["id"])]
        color_idx = unique_cl.index(c)
        node['color'] = get_hex_color(color_idx, cmap)
        node['title'] = "<h1>" + str(c) + "</h1>"


options = """var options = {
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -30000,
      "centralGravity": 0,
      "springLength": 100
    },
    "minVelocity": 0.75
  }
}
"""
net.toggle_physics(True)
net.set_options(options)


net.show("mygraph.html")

# import bs4
#
# # load the file
# with open("mygraph.html") as inf:
#     txt = inf.read()
#     soup = bs4.BeautifulSoup(txt)
#
# # create new link
# div = soup.new_tag("div")
# # insert it into the document
# div.append(soup.new_tag("textarea", id="MyTextarea"))
# div.append(soup.new_tag("br"))
# btn1 = soup.new_tag("button", onclick="clicked()")
# btn1.string = "update"
# div.append(btn1)
# btn2 = soup.new_tag("button", id="prevBtn", onclick="clicked()")
# btn2.string = "<"
# div.append(btn2)
# btn3 = soup.new_tag("button", id="forwBtn", onclick="clicked()")
# btn3.string = ">"
# div.append(btn3)
# div.append(soup.new_tag("br"))
# div.append(soup.new_tag("input", type="range", id="iterSlider", min="0", max="10"))
# div.append(soup.new_tag("br"))
#
# soup.body.append(div)
# # save the file again
# with open("mygraph.html", "w") as outf:
#     outf.write(str(soup))