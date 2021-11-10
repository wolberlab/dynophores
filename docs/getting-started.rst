Getting started
===============

IO data
-------

Input
~~~~~

The ``dynophores`` package will generate a dynophore based on the following input data:

* ``startframe.(pdb|pmz)``: Topology as PDB file or PMZ file.

    The PMZ file contains the binding site of interest and assigns atoms to either the macromolecule 
    or the ligand side. Interactions will be detected between these two sides.
    Such PMZ files can be generated with the software LigandScout.

    If you do not have such a PMZ file, you can also start from a PDB file. 
    The ``dynophores`` package either infers the binding site automatically or you specify 
    the chain and ligand of interest yourself.

* ``trajectory.dcd``: Trajectory as DCD file.

    Make sure that the trajectory is matching the input PDB or PMZ file.

* ``feature_definitions.xml`` (optional!): Set customized feature definitions.

    Change the default feature definitons by editing `this file <https://github.com/wolberlab/dynophores/tree/master/dynophores/data/custom-chemicalfeature-definitions.xml>`_.

Output
~~~~~~

The ``dynophores`` package will generate the following dynophore data::

    dynophore_out_YYYY-MM-DD_hh-mm-ss/
    ├── dynophore.pml     # Point cloud data
    ├── dynophore.cgo     # Point cloud data as PyMol script
    ├── dynophore.json    # Statistics
    └── dynophore.ipynb   # Visualization notebook

Data from the PML and JSON files are visualized in the interactive Jupyter notebook ``dynophore.ipynb``.

Installation
------------

Install the ``dynophores`` package as described `here <https://dynophores.readthedocs.io/en/latest/installing.html>`_.

Usage
-----

* Take a look at a demo notebook showing the dynophore from an MD simulation for ligand-bound kinase CDK2 (PDB ID: 1KE7)::

    dynophore demo path/to/output/folder

  Generates a dynophore demo notebook at ``path/to/output/folder/dynophore.ipynb``.

* If you have an input PDB (topology) and DCD (trajectory) file, run::

    dynophore create -p path/to/pdb/file -d path/to/dcd/file -o path/to/output/folder -n dynophore_name

  Generates dynophore data at ``path/to/output/folder/dynophore_out_YYYY-MM-DD_hh-mm-ss``, which contains the dynophore notebook.

* If you have an input PDB (topology) and DCD (trajectory) file *and* you want to specify the chain or ligand to be used for the binding site definiton, run::

    dynophore create -p path/to/pdb/file -c chain_id -3 three_letter_ligand_code -d path/to/dcd/file -o path/to/output/folder -n dynophore_name

  Generates dynophore data at ``path/to/output/folder/dynophore_out_YYYY-MM-DD_hh-mm-ss``, which contains the dynophore notebook.

* If you have an input PMZ (topology) and DCD (trajectory) file, run::

    dynophore create -p path/to/pmz/file -d path/to/dcd/file -o path/to/output/folder -n dynophore_name

  Generates dynophore data at ``path/to/output/folder/dynophore_out_YYYY-MM-DD_hh-mm-ss``. 
  To generate the corresponding dynophore notebook, run ``dynophore visualize``.

* If you already have a dynophore data folder, run::

    dynophore visualize -i path/to/dynophore/folder -p path/to/pdb/file -d path/to/dcd/file

  Adds the dynophore notebook to the dynophore data files at  
   ``path/to/dynophore/folder/dynophore.ipynb``.

