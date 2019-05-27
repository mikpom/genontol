import unittest
from pkg_resources import resource_filename as pkg_file
from genontol.ontol import GOntology, GOTermNotFoundError, _propagate
import networkx as nx
import pandas as pd

class TestGOntologyClass(unittest.TestCase):
    ont = GOntology.from_obo(pkg_file('genontol.tests', 'data/goslim_generic.obo.gz'))
    annot_df = pd.read_csv(pkg_file('genontol.tests', 'data/gene2go.tsv'), sep='\t')
    go2gene = {k: set(v) for k,v in annot_df.groupby('term')['gene']}

    def test_creation(self):
        self.assertTrue('data/goslim_generic.obo.gz' in self.ont.obofile)
        e = self.ont.G.get_edge_data('GO:0000228', 'GO:0005634')
        self.assertEqual(e['rel'], 'part_of')
        self.assertDictEqual(self.ont.roots, {'biological_process': 'GO:0008150',
                                            'cellular_component': 'GO:0005575',
                                            'molecular_function': 'GO:0003674'})

    def test_get_term(self):
        t = self.ont.get_term('GO:0003729')
        self.assertEqual(t.name, 'mRNA binding')
        with self.assertRaises(GOTermNotFoundError):
            self.ont.get_term('GO:1234567')

    def test_is_child_parent(self):
        self.assertTrue(self.ont.is_child_parent('GO:0065003', 'GO:0022607'))
        self.assertFalse(self.ont.is_child_parent('GO:0065003', 'GO:0000902'))

    def test_get_child_terms(self):
        self.assertListEqual(self.ont.get_child_terms('GO:0009056'),
                             list(map(lambda t: self.ont.get_term(t),
                                      ['GO:0006914', 'GO:0034655'])))
    def test_get_parent_terms(self):
        self.assertListEqual(self.ont.get_parent_terms('GO:0034655'),
                             list(map(lambda t: self.ont.get_term(t),
                                      ['GO:0009056', 'GO:0034641'])))

    def test_search_term_by_name(self):
        self.assertListEqual(self.ont.search_terms_by_name('mitochond'),
                             list(map(lambda t: self.ont.get_term(t),
                                      ['GO:0005739', 'GO:0007005'])))

    def test_getattr(self):
        self.assertEqual(self.ont.get_attr('GO:0140014', 'name'),
                         'mitotic nuclear division')

    def test_all_terms(self):
        self.assertEqual(len(list(self.ont.all_terms())), 149)

    def test_get_relation(self):
        # GO:0006464 - cellular protein modification process
        # GO:0008150 - biological_process
        # GO:0000902 - cell morphogenesis
        # GO:0048856 - anatomical structure development

        # only one is_a path
        self.assertEqual(self.ont.get_relation('GO:0006464', 'GO:0008150'), 'is_a')
        self.assertEqual(self.ont.get_relation('GO:0006464', 'GO:0008150', True),
                         '<GO:0006464 cellular protein modification process> is_a '\
                         +'<GO:0008150 biological_process>')
        
        # only one part_of path
        self.assertEqual(self.ont.get_relation('GO:0000902', 'GO:0048856'), 'part_of')
        
        self.assertEqual(self.ont.get_relation('GO:0006464', 'GO:0008150', True),
                         self.ont.get_relation('GO:0008150', 'GO:0006464', True))
        
        with self.assertRaises(GOTermNotFoundError):
            self.ont.get_relation('GO:0006464', 'GO:0008151')

        # two is_a edges via GO:0065003
        self.assertEqual(self.ont.get_relation('GO:0022618', 'GO:0022607'), 'is_a')

        # is_a directly but part_of x is_a via GO:0048856
        # (anatomical structure development)
        self.assertEqual(self.ont.get_relation('GO:0000902', 'GO:0008150'), 'is_a;part_of')

        # there are five paths from <GO:0005840 ribosome> to <GO:0005575 cellular_component>
        # GO:0005840 is_a GO:0043226 is_a GO:0005575
        # GO:0005840 is_a GO:0032991 is_a GO:0005575
        # GO:0005840 part_of GO:0005737 is_a GO:0005575
        # GO:0005840 part_of GO:0005737 part_of GO:0005622 is_a GO:0005575
        # GO:0005840 part_of GO:0005737 part_of GO:0005622 part_of GO:0005623 is_a GO:0005575
        self.assertEqual(self.ont.get_relation('GO:0005840', 'GO:0005575'), 'is_a;part_of')

        # only one part_of path
        self.assertEqual(self.ont.get_relation('GO:0048856', 'GO:0048856'), 'same_term')
        self.assertEqual(self.ont.get_relation('GO:0048856', 'GO:0048856', True),
                         '<GO:0048856 anatomical structure development> same_term ' \
                         +'<GO:0048856 anatomical structure development>')

    def test_propagation(self):
        self.ont.propagate(self.go2gene, 'prp')
        self.assertEqual(len(self.ont.G.node['GO:0051276']['prp']), 305)
        self.assertEqual(len(self.ont.G.node['GO:0034641']['prp']), 459)

    def test_enrichment_analysis(self):
        self.ont.propagate(self.go2gene, 'bg')

        # GO:0051276 chromosome organization
        t0 = 'GO:0051276'
        genes_w_t0 = self.ont.get_term(t0).node['bg']
        genes_wo_t0 =  set(self.annot_df['gene']).difference(genes_w_t0)
        test_query = list(genes_w_t0)[:100] + list(genes_wo_t0)[:1900]
        
        df = self.ont.get_enrichment(test_query, 'bg')
        self.assertAlmostEqual(0.0063730600891504039, df[df['term']==t0]['p'].values[0])

        df = self.ont.get_enrichment(test_query, 'bg', domain='biological_process')
        self.assertAlmostEqual(0.0063730600891504039, df[df['term']==t0]['p'].values[0])

        self.assertEqual((df['p']>df['q']).sum(), 0)

class testToyExample(unittest.TestCase):
    def test_pval_correctness(self):
        G = nx.DiGraph()
        G.add_node(0, name='r', namespace='a')
        G.add_node(1, name='1', namespace='a')
        G.add_node(2, name='2', namespace='a')
        G.add_edge(0, 1)
        G.add_edge(0, 2)
        O = GOntology.from_graph(G)
        values = {1: set(range(10)), 2: set(range(20))}
        background_attribute = 'bg_attr'
        O.propagate(values, background_attribute)
        query = [1, 2, 3, 4, 5]
        df = O.get_enrichment(query, background_attribute)
        best_pval = float(df.dropna().sort_values('p').head(1).p)
        self.assertAlmostEqual(best_pval, 0.016253869969040255)

if __name__ == '__main__':
    unittest.main()
