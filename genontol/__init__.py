import networkx as nx

if float(nx.__version__)<2.0:
    raise OSError("Networkx >= 2.0 required")
