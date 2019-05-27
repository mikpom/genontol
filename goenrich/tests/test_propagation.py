# Copyright 2015-2019 Jan Daniel Rudolph (@jdrudolph)

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
import unittest
import networkx as nx
import pandas as pd
import goenrich
from goenrich.enrich import propagate


class TestPropagationExample(unittest.TestCase):
    def test_correctness_on_small_example(self): 
        r""" Example graph
            r
          /   \
        c1     c2
          \   /  \
           \ /    \
            c3    c4
        """
        O = nx.DiGraph([('c4', 'c2'), ('c3', 'c1'), ('c3', 'c2'),
            ('c1', 'r'), ('c2', 'r')])
        
        r = set([6])
        c1 = set([])
        c2 = set([4,5])
        c3 = set([1,2,3])
        c4 = set([0])
        x = { 'r' : r, 'c1' : c1, 'c2' : c2, 'c3' : c3, 'c4' : c4 }

        b = 'background'
        propagate(O, x, b)

        self.assertSetEqual(O.node['c3'][b], c3)
        self.assertSetEqual(O.node['c2'][b], c4 | c3 | c2)
        self.assertSetEqual(O.node['c1'][b], c3 | c1)
        self.assertSetEqual(O.node['r'][b], c4 | c3 | c2 | c1 | r)

class TestPropagationReal(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestPropagationReal, self).__init__(*args, **kwargs)
        self.gene2go = goenrich.read.gene2go('db/gene2go.gz')
        self.O = goenrich.obo.ontology('db/go-basic.obo')

    def test_on_gene2go_head(self):
        test = self.gene2go.head(100)
        values = {k: set(v) for k,v in test.groupby('GO_ID')['GeneID']}
        propagate(self.O, values, 'head')

    def test_if_runs_trough_on_real_data(self):
        values = {k: set(v) for k,v in self.gene2go.groupby('GO_ID')['GeneID']}
        propagate(self.O, values, 'real')

if __name__ == '__main__':
    unittest.main()
