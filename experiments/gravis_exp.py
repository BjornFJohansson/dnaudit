import gravis as gv
import networkx as nx
from pathlib import Path

g = nx.Graph()
g.add_node('Node1')
g.add_node('Node2')
g.add_edge('Node1', 'Node2')

g.nodes['Node1']["name"] = f"MyNode1"
g.nodes['Node1']["hover"] = '<a href="snippet1.md">snippet1</a>'

g.nodes['Node2']["name"] = f"MyNode2"
g.nodes['Node2']["hover"] = '<a href="snippet2.md">snippet2</a>'

fig = gv.d3(g, node_label_data_source='name')

page = fig.to_html_standalone()

Path("gravis.html").write_text(page)
