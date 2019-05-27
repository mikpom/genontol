========
Tutorial
========

Create and use the Ontology
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

  import genontol

  # Create ontology object
  # you can get .obo here http://purl.obolibrary.org/obo/go/go-basic.obo
  O = genontol.ontol.GOntology.from_obo('/path/to/go-basic.obo')

  # getting a term
  t = O.get_term('GO:0003729')
  # <GO:0003729 mRNA binding>

  # getting GO term child terms
  tt = O.get_child_terms(t) #or with Term ID O.get_child_terms('GO:0003729')
  # [<GO:0003730 mRNA 3'-UTR binding>,
  #  <GO:0030350 iron-responsive element binding> ...]

  # search GO term by name
  O.search_terms_by_name('translation factor activity')
  # [<GO:0008135 translation factor activity, RNA binding>,
  #  <GO:0045183 translation factor activity, non-nucleic acid binding>]


Run GO enrichment analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

  # use all geneid associations form HUMAN gaf file as background
  # HUMAN data can be downloaded here:
  # http://geneontology.org/gene-associations/goa_human.gaf.gz
  goa_df = genontol.read.goa('/path/to/goa_human.gaf.gz')
  go2prot = {k: set(v) for k,v in goa_df.groupby('go_id')['db_object_id']}

  # propagate the background through the ontology
  background_attribute = 'human_gaf'
  O.propagate(go2prot, background_attribute)

  # extract large subunit ribosomal proteins as an example query
  query = goa_df[goa_df.db_object_name.apply(lambda n: n.startswith('60S ribosomal'))]\
          .db_object_id.unique()

  # top category should be ribosome related
  df = O.get_enrichment(query, background_attribute)
  
