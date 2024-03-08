from graphviz import Digraph
dot = Digraph(comment="First ")
dot.attr(rankdir="LR")
dot.node("Snippet1", href="snippet1.md", shape = "note")
dot.node("Snippet2", href="snippet2.txt", shape = "folder")
dot.edge("Snippet1", "Snippet2", href="grp.py")
dot.edge("Snippet1", "Snippet2", href="grp.py")

dot.render(format="svg", cleanup=True, filename="grp")
