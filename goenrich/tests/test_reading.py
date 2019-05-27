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
import pkg_resources
import unittest
import gzip
import goenrich

class TestRead(unittest.TestCase):
    def test_ontology(self):
        G = goenrich.obo.ontology('db/go-basic.obo')

    def test_ontology_from_file_obj(self):
        with open('db/go-basic.obo') as f:
            G = goenrich.obo.ontology(f)
            self.assertFalse(f.closed)

    def test_goa(self):
        background = goenrich.read.goa('db/gene_association.goa_human.gaf.gz')

    def test_goa_from_file_obj(self):
        with gzip.GzipFile('db/gene_association.goa_human.gaf.gz') as f:
            background = goenrich.read.goa(f)
            self.assertFalse(f.closed)

    def test_gene2go(self):
        background = goenrich.read.gene2go('db/gene2go.gz')

    def test_gene2go_from_file_obj(self):
        with gzip.GzipFile('db/gene2go.gz') as f:
            background = goenrich.read.gene2go(f)
            self.assertFalse(f.closed)

    def test_goslim_from_file(self):
        G = goenrich.obo.ontology(pkg_resources.resource_filename(goenrich.__name__, 'tests/test_ontologies/goslim_generic.obo'))
        self.assertEqual(len(G.nodes()), 150)
        self.assertSetEqual(set(G.successors('GO:0009056')), set(['GO:0008150']))
        self.assertSetEqual(set(G.predecessors('GO:0009056')), set(['GO:0034655', 'GO:0006914']))

if __name__ == '__main__':
    unittest.main()
