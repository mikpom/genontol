from collections import defaultdict
import networkx as nx
from networkx import NetworkXNoPath
from genontol.obo import read_obo
from genontol.enrich import analyze as _analyze
import os

class GOTermNotFoundError(Exception):
    def __init__(self, termid, *args, **kwargs):
        self.termid = termid
        Exception.__init__(self, *args, **kwargs)


def _propagate(O, values, attribute):
    """ Propagate values trough the hierarchy"""
    for n in nx.topological_sort(O):
        current = O.node[n].setdefault(attribute, set())
        current.update(values.get(n, set()))
        for p in O[n]:
            O.node[p].setdefault(attribute, set()).update(current)

class GOterm(object):
    """
    Class implementing GO term.
    """

    def __init__(self, id, name='', definition='', n=None,
                 namespace=None):
        self.id = id
        self.name = name
        self.definition = definition
        self.node = n
        self.namespace = namespace

    def __repr__(self):
        return '<'+self.id + ' ' + self.name+'>'

    def __str__(self):
        return '<GO Term ' + self.id + '>'

    def __eq__(self, other):
        if ((self.id == other.id) and \
            (self.name == other.name) and \
            (self.definition == other.definition)):
            return True
        else:
            return False

class GOntology(object):
    """
    Class representing the Ontology.
    Can be instatiated from an obo file:

    >>> O = GOntology.from_obo('/path/to/obo/file')

    Or directly from networkx graph:

    >>> O = GOntology.from_graph(G)

    Ontology object methods accepting term as an argument accept
    :class:`GOterm` objects as well as term IDs of the form
    GO\:XXXXXXX. Ontology object methods returning terms return
    objects of :class:`GOterm` class.

    """

    def __init__(self, G, obofile=None, roots=None):
        self.G = G
        self.roots = roots
        if obofile:
            self.obofile = os.path.abspath(obofile)
        else:
            self.obofile = None

    @classmethod
    def from_obo(cls, obofile):
        """
        Reads GOntology from obo file privided.
        """
        G = read_obo(obofile)
        return cls(G, obofile, roots=G.graph['roots'])

    @classmethod
    def from_graph(cls, G):
        """
        Creates GOntology instance from a networkx DiGraph.
        """
        return cls(G)

    def __repr__(self):
        if self.obofile:
            return '<GO Ontology f:"{:s}" at {:s}>'.format(self.obofile, hex(id(self)))
        else:
            return '<GO Ontology at {:s}>'.format(hex(id(self)))

    def get_term(self, termid):
        """
        Get :class:`GOterm` object by its id of the form GO\:XXXXXXX.
        """
        if (not isinstance(termid, str)):
            raise TypeError('termid parameter should be a string of the form GO:XXXXXXX')
        if not self.G.has_node(termid):
            raise GOTermNotFoundError(termid, 'Term with id {:s} not found in ontology {:s}'\
                                      .format(termid, repr(self)))
        n = self.G.node[termid]
        term = GOterm(id=termid, name=n['name'], n=n,
                      namespace=n['namespace'])
        return term

    def has_term(self, termid):
        if self.G.has_node(termid):
            return True
        else:
            return False

    def _check_terms(self, *terms):
        for term in terms:
            if (type(term) == str):
                pass
            elif (isinstance(term, GOterm)):
                term = term.id
            if not self.G.has_node(term):
                raise GOTermNotFoundError(term, 'Term with id {:s} not found in ontology {:s}'\
                                      .format(term, repr(self)))
        return

    def all_terms(self):
        for n in self.G.nodes:
            yield n
            
    def get_attr(self, termid, attr):
        """
        Get attribute attr for a term
        """
        return self.G.node[termid][attr]

    def get_relation(self, t1, t2, verbose=False):
        # print('getting relation for pair', t1, t2)
        t1 = getattr(t1, 'id', t1)
        t2 = getattr(t2, 'id', t2)

        rel = None

        if t1==t2:
            rel = 'same_term'
        elif self.is_child_parent(t1, t2):
            pass
        elif self.is_child_parent(t2, t1):
            t1, t2 = t2, t1
        else:
            rel = 'unrelated'

        if not rel:
            pp = nx.all_simple_paths(self.G, t1, t2)
            path_types = set()
            relations = set()
            for p in pp:
                #print("processing path", p)
                path_rels = set()
                for i in range(len(p)-1):
                    d = self.G.edges[p[i], p[i+1]]
                    path_rels.add(d['rel'])
                path_types.add(tuple(sorted(path_rels)))
            for path_type in path_types:
                if path_type == ('is_a', 'part_of'):
                    relations.add('part_of')
                elif path_type == ('is_a',):
                    relations.add('is_a')
                elif path_type == ('part_of',):
                    relations.add('part_of')
            rel = ';'.join(sorted(relations))
        if not verbose:
            return rel
        else:
            return repr(self.get_term(t1))+' '\
                   + rel +' '\
                   +repr(self.get_term(t2))

    def is_child_parent(self, t1, t2):
        """
        Return True if terms t1 and t2 are related by child parent GO
        relashionship. Return False otherwise.

        Parent term is the less specific one and child term is the
        more specific one.  t1 and t2 don't have to be immediate
        child-parent pair but can be separated by a term chain in the
        GO graph.
        """
        self._check_terms(t1, t2)
        t1 = getattr(t1, 'id', t1)
        t2 = getattr(t2, 'id', t2)
        try:
            p = nx.shortest_path(self.G, t1, t2)
            if p:
                return True
        except NetworkXNoPath:
            return False    

    def get_child_terms(self, term):
        """
        Get immediate child terms of a term
        """

        self._check_terms(term)
        term = getattr(term, 'id', term)
        childids = sorted(self.G.predecessors(term))
        return [self.get_term(childid) for childid in childids]

    def get_parent_terms(self, term):
        """
        Get immediate parent terms of a term
        """

        self._check_terms(term)
        term = getattr(term, 'id', term)
        parentids = sorted(self.G.successors(term))
        return [self.get_term(parentid) for parentid in parentids]

    def search_terms_by_name(self, pattern):
        """
        Search the ontology for terms by pattern

        :param pattern: a string to scan for in GO terms' names
        """

        ids_found = filter(lambda tid: pattern in self.G.node[tid]['name'],
                           self.G.nodes())
        terms_found = [self.get_term(i) for i in sorted(ids_found)]
        return terms_found

    def propagate(self, values, attribute):
        """ Propagate values through the ontology graph.

        >>> import genontol
        >>> O = genontol.ontol.GOntology.from_obo('/path/to/go-basic.obo')
        >>> goa_df = genontol.read.goa('/path/to/goa_human.gaf.gz')
        >>> go2prot = {k: set(v) for k,v in goa_df.groupby('go_id')['db_object_id']}
        >>> O.propagate(values, 'bg')

        Uses topological sorting of the vertices. Since degrees are
        usually low performance is almost linear time.

        :param values: mapping of nodes id to set of ids,
                       i.e. dict {node_id:{id1, id2, ...}, ...}
        :param attribute: name of the attribute
        """
        _propagate(self.G, values, attribute)

    def del_attr(self, attr):
        for node in self.G.nodes():
            if attr in self.G.node[node]:
                del self.G.node[node][attr]

    def remove_attr(self, attribute):
        for node in self.G.nodes():
            if attribute in self.G.node[node]:
                del self.G.node[node][attribute]

    def get_enrichment(self, query, background_attr, domain=None, **kwargs):
        """ run enrichment analysis for query

        >>> import genontol
        >>> O = genontol.ontol.GOntology.from_obo('/path/to/go-basic.obo')
        >>> goa_df = genontol.read.goa('/path/to/goa_human.gaf.gz')
        >>> go2prot = {k: set(v) for k,v in goa_df.groupby('go_id')['db_object_id']}
        >>> O.propagate(values, 'bg')
        >>> df = O.get_enrichment(query, 'bg')

        :param query: array like of ids
        :param background_attr: string idicating background attribute to use
        :returns: pandas.DataFrame with results
        """
        if domain:
            items = filter(lambda i: i[1]['namespace']==domain, self.G.nodes.items())
            nodes = [i[0] for i in items]
            subgraph = self.G.subgraph(nodes)
            M = len({x for n in self.G.nodes.values() for x in n[background_attr]})
            return _analyze(subgraph, query, background_attr, M=M, **kwargs)
        else:
            return _analyze(self.G, query, background_attr, **kwargs)
