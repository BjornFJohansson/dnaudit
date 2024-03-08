#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code prouduces a JavaScript representation of a graph
using the software package gravis (https://github.com/robert-haas/gravis)

Gravis uses JavaScript libraries d3.js, vis.js and 3d-force-graph.js / three.js.

Gravis understands graph objects from graph-tool, igraph, NetworKit, NetworkX, pyntacle and SNAP.

This example uses NetworkX.

To run this script, first install gravis & networkx:

    pip install gravis networkx

Execute this script in this folder and a file called graph.html will be generated.



"""




import gravis as gv
import networkx as nx
from pathlib import Path

g = nx.DiGraph()


nodes = """\
382.fasta
383.fasta
pAG25.fasta
natMX6_cassette_pAG25_PCR.txt
cassette.fasta
PGI_locus.fasta
natMX6_PGI1_RECOMBINATION.txt
PGI_kanMX4_locus.fasta
""".splitlines()

for node in nodes:
    g.add_node(node)
    g.nodes[node]["name"] = node
    g.nodes[node]["hover"] = f'<a href="{node}">{node}</a>'

g.add_edge("382.fasta", "natMX6_cassette_pAG25_PCR.txt")
g.add_edge("383.fasta", "natMX6_cassette_pAG25_PCR.txt")
g.add_edge("pAG25.fasta", "natMX6_cassette_pAG25_PCR.txt")
g.add_edge("natMX6_cassette_pAG25_PCR.txt", "cassette.fasta")
g.add_edge("cassette.fasta", "natMX6_PGI1_RECOMBINATION.txt")
g.add_edge("PGI_locus.fasta", "natMX6_PGI1_RECOMBINATION.txt")
g.add_edge("natMX6_PGI1_RECOMBINATION.txt", "PGI_kanMX4_locus.fasta")

fig = gv.d3(g, node_label_data_source="name", zoom_factor=0.25)
fig = gv.three(g)
page = fig.to_html_standalone()
Path("graph.html").write_text(page)
