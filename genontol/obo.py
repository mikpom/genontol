# Copyright 2015-2019 Jan Daniel Rudolph (@jdrudolph)
# Modified version Copyright 2019 Mikhail Pomaznoy (@mikpom)

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import itertools
import networkx as nx
import gzip

def _tokenize(f):
    token = []
    for line in f:
        if line == '\n':
            yield token
            token = []
        else:
            token.append(line)

def _filter_terms(tokens):
    for token in tokens:
        if token[0] == '[Term]\n':
            yield token[1:]

def _parse_terms(terms):
    for term in terms:
        obsolete = False
        node = {}
        parents = {}
        for line in term:
            if line.startswith('id:'):
                id = line[4:-1]
            elif line.startswith('name:'):
                node['name'] = line[6:-1]
            elif line.startswith('namespace:'):
                node['namespace'] = line[11:-1]
            elif line.startswith('is_a:'):
                parents[line[6:16]] = 'is_a'
            elif line.startswith('relationship: part_of'):
                parents[line[22:32]] = 'part_of'
            elif line.startswith('is_obsolete'):
                obsolete = True
                break
        if not obsolete:
            edges = [(p, id, {'rel':r}) for p, r in parents.items()] # will reverse edges later
            yield (id, node), edges
        else:
            continue

_filename = 'db/go-basic.obo'

def read_obo(file):
    """ read ontology from file
    :param file: file path of file handle
    """
    G = nx.DiGraph()

    if isinstance(file, str):
        if file.endswith('.gz'):
            f = gzip.open(file, 'rt')
        else:
            f = open(file)
        we_opened_file = True
    else:
        f = file
        we_opened_file = False

    try:
        tokens = _tokenize(f)
        terms = _filter_terms(tokens)
        entries = _parse_terms(terms)
        nodes, edges = zip(*entries)
        G.add_nodes_from(nodes)
        G.add_edges_from(itertools.chain.from_iterable(edges))
        G.graph['roots'] = {data['name'] : n for n, data in G.node.items()
                            if data['name'] == data['namespace']}
    finally:
        if we_opened_file:
            f.close()

    for root in G.graph['roots'].values():
        for n, depth in nx.shortest_path_length(G, root).items():
            node = G.node[n]
            node['depth'] = min(depth, node.get('depth', float('inf')))
    return G.reverse()
