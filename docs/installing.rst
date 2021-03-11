Installing
==========

[WIP 🚧] Install from the conda package
---------------------------------------

Eventually, we will have a ``conda`` package, but for now you need to create a new environment manually.

.. note::

    We are assuming you have a working ``conda`` installation in your computer. 
    If this is not the case, please refer to the `official documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation>`_. 
    We recommend using Miniconda.


1. Create a new conda environment::

    conda env create -f https://raw.githubusercontent.com/dominiquesydow/dynophores/master/devtools/conda-envs/test_env.yaml -n dynophores

   Installing via ``conda`` often takes a while. Instead you can also use ``mamba``, a faster reimplementation of the conda package manager in C++.
   Find ``mamba`` installation instructions `here <https://mamba.readthedocs.io/en/latest/getting_started.html#for-new-users>`_.
   If you have ``mamba`` installed successfully, create your conda environment within a couple of seconds::

    mamba env create -f https://raw.githubusercontent.com/dominiquesydow/dynophores/master/devtools/conda-envs/test_env.yaml -n dynophores

2. Activate the new environment and install a few Jupyter extensions we will need for our dynophore notebook::

    conda activate dynophores
    jupyter labextension install @jupyter-widgets/jupyterlab-manager nglview-js-widgets @jupyterlab/toc


   Note: We do not have a ``dynophores`` conda package yet. Until we do, you will need to perform the additional steps 
   described below under "Install from the latest development snapshot" before you can continue.

3. Test that your installation works::

    dynoviz -h

3. Take a look at a demo notebook showing the dynophore from an MD simulation for ligand-bound kinase CDK2 (PDB ID: 1KE7)::

    dynoviz demo

4. Explore your own dynophore data in a new notebook::

    dynoviz create --dyno path/to/dyno/folder --pdb path/to/pdb/file --dcd path/to/dcd/file --workspace path/to/workspace/folder

5. If you already have set up your dynophore notebook (i.e. if you want to revisit a notebook created in step 4), run::

    dynoviz open path/to/your/dyno/notebook

   Note: This command is equivalent to::
    
    jupyter lab path/to/your/dyno/notebook

Install from the latest development snapshot
--------------------------------------------

If you have already created a *conda environment* and it has been activated  (see above steps 1+2), 
the next step is downloading a copy of the current state of the 
`GitHub repository <https://github.com/dominiquesydow/dynophores>`_.

1. Download a zipfile of the repository using `this link <https://github.com/dominiquesydow/dynophores/archive/master.zip>`_.
2. Unzip to your location of choice.
3. Navigate to ``path/to/directory/that/contains/downloaded/folder``.
4. Run ``pip install dynophores-master``. 

   * If the ``pip`` installation failed, try ``python setup.py install`` from within the ``dynophores-master`` folder. 
   * If you cloned instead of downloaded the repository, your folder is called ``dynophores`` not ``dynophores-master``.
5. Start your dynophore notebook as described above in steps 3-5.


.. Unix instructions

.. raw:: html

    <details>
    <summary>Instructions for Linux / MacOS</summary>

.. code-block:: bash

    wget https://github.com/dominiquesydow/dynophores/archive/master.zip -O dynophores.zip
    mkdir -p ~/Documents
    unzip dynophores.zip -d ~/Documents
    cd ~/Documents
    pip install dynophores-master/
    dynoviz -h

.. raw:: html

    </details>

.. Windows instructions

.. raw:: html

    <details>
    <summary>Instructions for Windows (PowerShell)</summary>

.. code-block::

    wget https://github.com/dominiquesydow/dynophores/archive/master.zip -O dynophores.zip
    mkdir ~/Documents/
    Expand-Archive dynophores.zip -d ~/Documents
    cd ~/Documents
    pip install dynophores-master/
    dynoviz -h

.. raw:: html

    </details>