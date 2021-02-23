Installing
==========

Install from the conda package
------------------------------

Eventually, we will have a ``conda`` package, but for now you need to create a new environment manually.

.. note::

    We are assuming you have a working ``conda`` installation in your computer. 
    If this is not the case, please refer to the official documentation 
    `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation>`_. 
    We recommend using Miniconda.


1. Create a new conda environment::

    conda env create -f https://raw.githubusercontent.com/dominiquesydow/dynophores/master/environment.yml -n dynophores

2. Activate the new environment::

    conda activate dynophores

3. Run ``dynophores -h`` to test that it works.
4. [WIP 🚧] Run ``dynophores start .`` to set up a new workspace. This step is not fully functional yet. 
Please refer to "Install from the latest development snapshot" for manual steps.


Install from the latest development snapshot
--------------------------------------------

If you have already created a *conda environment* and it has been activated  (see above), 
the next step is downloading a copy of the current state of the 
`GitHub repository <https://github.com/dominiquesydow/dynophores>`_.

1. Download a zipfile of the repository using `this link <https://github.com/dominiquesydow/dynophores/archive/master.zip>`_.
2. Unzip to your location of choice.
3. Navigate to ``path/to/your/location/dynophores-master/docs``.
4. Start Jupyter Lab.
5. Double click on the lesson you want to start.


.. Unix instructions

.. raw:: html

    <details>
    <summary>Instructions for Linux / MacOS</summary>

.. code-block:: bash

    wget https://github.com/dominiquesydow/dynophores/archive/master.zip -O dynophores.zip
    mkdir -p ~/Documents
    unzip dynophores.zip -d ~/Documents
    cd ~/Documents/dynophores-master/docs
    jupyter lab

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
    cd ~/Documents/dynophores-master/docs
    jupyter lab

.. raw:: html

    </details>