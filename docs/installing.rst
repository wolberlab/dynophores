Installing
==========

.. note::

    We are assuming you have a working ``mamba`` installation in your computer. 
    If this is not the case, please refer to their `official documentation <https://mamba.readthedocs.io/en/latest/installation.html#mamba>`_. 

    If you installed ``mamba`` into an existing ``conda`` installation, also make sure that the ``conda-forge`` channel is configured by running ``conda config --add channels conda-forge``. 


Install from the conda package
------------------------------

1. Create a new conda environment called ``dyno`` with the ``dynophores`` package and all its dependencies installed::

    mamba create -n dyno dynophores

2. Activate the new conda environment::

    conda activate dyno

3. Test that your installation works::

    dynoviz -h


Install from the latest development snapshot
--------------------------------------------

Install the latest development snapshot from the `GitHub repository's master branch <https://github.com/wolberlab/dynophores>`_.


1. Create a new conda environment called ``dyno``::

    mamba env create -f https://raw.githubusercontent.com/wolberlab/dynophores/master/devtools/conda-envs/test_env.yaml -n dyno

2. Activate the new conda environment::

    conda activate dyno

3. Install ``dynophores`` package via pip::

    pip install https://github.com/wolberlab/dynophores/archive/refs/heads/master.zip

4. Test that your installation works::

    dynoviz -h
