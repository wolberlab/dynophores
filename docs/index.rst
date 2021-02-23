.. dynophores documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Dynophores
==========

Dynamic pharmacophore modeling of molecular interactions

.. raw:: html

   <p align="center">
   <img src="_static/dynophore.png" alt="Dynophore 3D visualization" width="400"/>
   <br>
   <font size="1">
   Dynophore 3D visualization.
   </font>
   </p>

Dynophores are an intuitive tool to explore the dynamics of molecular interactions between a 
macromolecule and a ligand throughout an MD simulation. 

A **dynophore** is a collection of so-called **superfeatures**. 
A superfeature is defined as a pharmacophore feature on ligand site 
(defined by a feature type, e.g. HBA, and one or more ligand atom numbers/serials) that occurs 
at least once during an MD simulation. A superfeature can have one or more interaction partner(s) 
on macromolecule-side. 
These interaction partners are called **environmental partners**. 
Superfeatures and their environmental partners are monitored throughout an MD simulation.


This Python library offers visualization schemes for dynophore data generated with the 
``DynophoreApp`` software (TODO: add link to software): 

- Plot statistics on superfeatures and interactions (TODO: link to tutorial)
- Show superfeature clouds in 3D using the NGLviewer (TODO: link to tutorial)
- Show superfeatures mapped onto the ligand in 2D (TODO: link to tutorial)

.. toctree::
   :maxdepth: 1
   :caption: User guide

   installing
   tutorials/dynophore

.. toctree::
   :maxdepth: 1
   :caption: Explore package
   
   tutorials/explore_data
   tutorials/explore_plots
   tutorials/explore_view3d

.. toctree::
   :maxdepth: 1
   :caption: Developers
   
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
