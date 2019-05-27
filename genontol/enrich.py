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

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from genontol.tools import fdrcorrection

def analyze(O, query, background_attribute, M=None, **kwargs):
    options = {'show' : 'top20'}
    options.update(kwargs)
    _query = set(query)
    terms, nodes = zip(*G.nodes(data=True))
    if not M:
        M = len({x for n in nodes for x in n[background_attribute]}) # all ids used
    N = len(_query)
    ps, xs, ns = calculate_pvalues(nodes, _query, background_attribute,
            M, **options)
    qs, rejs = multiple_testing_correction(ps, **options)
    names, namespaces = zip(*[(n['name'], n['namespace']) for n in nodes])
    df = pd.DataFrame({'name':names, 'namespace':namespaces,
                       'term':terms, 'q':qs, 'rejected':rejs,
                       'p':ps, 'x':xs, 'n':ns, 'M':M, 'N':N})
    return df

def calculate_pvalues(nodes, query, background_attribute, M,
        min_category_size=3, max_category_size=500,
        max_category_depth=5, **kwargs):
    """ calculate pvalues for all categories in the graph

    :param nodes: nodes dictionary from the ontology graph after background was set
    :param query: set of identifiers for which the p value is calculated
    :param background_attribute: node attribute assoc. with the background set
    :param M: background size, total number of genes in the data
    :param min_category_size: categories smaller than this number are ignored
    :param max_category_size: categories larger than this number are ignored
    :param max_category_depth: categories lower in the hierarchy (more specific) will be ignored
    :returns: pvalues, x, n
    """
    N = len(query)
    vals = []
    for node in nodes:
        category = node[background_attribute]
        n = len(category)
        hits = query.intersection(category)
        x = len(hits)
        if ((node.get('depth', 0) > max_category_depth)
            or (n <= min_category_size)
            or (n > max_category_size)):
            vals.append((float('NaN'), x, n))
        else:
            vals.append((hypergeom.sf(x-1, M, n, N), x, n))
    return [np.array(x) for x in zip(*vals)]

def multiple_testing_correction(ps, alpha=0.05,
        method='benjamini-hochberg', **kwargs):
    """ correct pvalues for multiple testing and add corrected `q` value

    :param ps: list of pvalues
    :param alpha: significance level default : 0.05
    :param method: multiple testing correction method [bonferroni|benjamini-hochberg]
    :returns (q, rej): two lists of q-values and rejected nodes
    """
    _p = np.array(ps)
    q = _p.copy()
    rej = _p.copy()
    mask = ~np.isnan(_p)
    p = _p[mask]
    if method == 'bonferroni':
        q[mask] = p * len(p)
        rej[mask] = q[mask] < alpha
    elif method == 'benjamini-hochberg':
        _rej, _q = fdrcorrection(p, alpha)
        rej[mask] = _rej
        q[mask] = _q
    else:
        raise ValueError(method)
    return q, rej
