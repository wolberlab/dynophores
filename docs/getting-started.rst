Getting started
===============

The ``dynophores`` package will generate the following dynophore data::

    dynophore_out/
    ├── dynophore.pml     # Point cloud data
    ├── dynophore.cgo     # Point cloud data as PyMol script
    ├── dynophore.json    # Statistics
    └── dynophore.ipynb   # Visualization notebook

Installation
------------

Install the ``dynophores`` package as described `here <https://dynophores.readthedocs.io/en/latest/installing.html>`_.

Usage
-----

a. Take a look at a demo notebook showing the dynophore from an MD simulation for ligand-bound kinase CDK2 (PDB ID: 1KE7)::

    dynophore demo path/to/some/folder

b. If you have an input PDB (topology) and DCD (trajectory) file, run::

    dynophore create -p path/to/pdb/file -d path/to/dcd/file -o path/to/output/folder -n dynophore_name

c. If you have an input PDB (topology) and DCD (trajectory) file *and* you want to specify the chain or ligand to be used for the binding site definiton, run::

    dynophore create -p path/to/pdb/file -c chain_id -3 three_letter_ligand_code -d path/to/dcd/file -o path/to/output/folder -n dynophore_name

d. If you have an input PMZ (topology of binding site) and DCD (trajectory) file, run::

    dynophore create -p path/to/pmz/file -d path/to/dcd/file -o path/to/output/folder -n dynophore_name

e. If you already have a dynophore data folder::

    dynophore visualize -d path/to/dynophore/folder -p path/to/pdb/file -d path/to/dcd/file