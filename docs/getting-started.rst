Getting started
===============

Generate dynophores
-------------------

Download
^^^^^^^^

TBA

Usage
^^^^^

TBA


Visualize dynophores
--------------------

Dynophores can be easily visualized in Jupyter Notebooks with the ``dynophores`` Python package.


Installation
^^^^^^^^^^^^

Install the ``dynophores`` package as described `here <https://dynophores.readthedocs.io/en/latest/installing.html>`_.

Usage
^^^^^

You have different options what to do next:

a. Take a look at a demo notebook showing the dynophore from an MD simulation for ligand-bound kinase CDK2 (PDB ID: 1KE7)::

    dynoviz demo path/to/some/folder

b. Explore your own dynophore data in a new notebook::

    dynoviz create --dyno path/to/dyno/folder --pdb path/to/pdb/file --dcd path/to/dcd/file --workspace path/to/workspace/folder

c. If you already have set up your dynophore notebook (i.e. if you want to revisit a notebook created with option b), run::

    dynoviz open path/to/your/dyno/notebook

   Note: This command is equivalent to::
    
    jupyter lab path/to/your/dyno/notebook