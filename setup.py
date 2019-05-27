from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='genontol',
      version='0.1',
      description='Gene Ontology analysis with Python',
      long_description=readme(),
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries'],
      keywords= ['GO', 'Gene Ontology', 'Biology', 'Enrichment',
          'Bioinformatics', 'Computational Biology'],
      url='https://github.com/mikpom/genontol',
      author='Mikhail Pomaznoy',
      author_email='pom@mailbox.org',
      license='MIT',
      packages=['genontol'],
      install_requires=[
          'numpy',
          'pandas >= 0.17.0',
          'scipy',
          'networkx >=2.0'],
      test_suite='nose.collector',
      tests_require=['nose', 'setuptools'],
      zip_safe=False)
