Dynophores
==========

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/wolberlab/dynophores/workflows/CI/badge.svg)](https://github.com/wolberlab/dynophores/actions?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/wolberlab/dynophores/branch/master/graph/badge.svg)](https://codecov.io/gh/wolberlab/dynophores/branch/master)
[![Documentation Status](https://readthedocs.org/projects/dynophores/badge/?version=latest)](https://dynophores.readthedocs.io)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/dynophores.svg)](https://anaconda.org/conda-forge/dynophores)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Dynamic pharmacophore modeling of molecular interactions

![OpenCADD](/docs/_static/dynophore.png)

Dynophores are an intuitive tool to explore the dynamics of molecular interactions between a 
macromolecule and a ligand throughout an MD simulation. 

A **dynophore** is a collection of so-called **superfeatures**. 
A superfeature is defined as a pharmacophore feature on ligand site 
(defined by a feature type, e.g. HBA, and one or more ligand atom numbers/serials) that occurs 
at least once during an MD simulation. A superfeature can have one or more interaction partner(s) 
on macromolecule-side. 
These interaction partners are called **environmental partners**. 
Superfeatures and their environmental partners are monitored throughout an MD simulation.

## Documentation

The documentation is available [here](https://dynophores.readthedocs.io/en/latest/), including [installation instructions](https://dynophores.readthedocs.io/en/latest/installing.html)

## Overview

This Python library offers visualization schemes for dynophore data generated with the 
`DynophoreApp` in a single Jupyter Notebook 
([template dynophore notebook](https://dynophores.readthedocs.io/en/latest/tutorials/dynophore.html)).  

- 3D view of the dynophore's point clouds (one cloud per superfeature) using `nglview` allows easy visual inspection of the dynamic macromolecule-ligand interactions. Point clouds are rendered alongside the topology and (optionally) the trajectory underlying the dynophore. See details in [this tutorial on 3D views](https://dynophores.readthedocs.io/en/latest/tutorials/explore_view3d.html).
- Statistics cover the occurrence of superfeatures and their environmental partners as well as distances between them. See details in [this tutorial on plotting options](https://dynophores.readthedocs.io/en/latest/tutorials/explore_plots.html).
- Dynophore data can be further analyzed conveniently right there in the same notebook by working with the `Dynophore` class. See details in [this tutorial on the dynophore data structure](https://dynophores.readthedocs.io/en/latest/tutorials/explore_data.html).


## License

`dynophores` is a free software and is licensed under the MIT license. 
Copyright (c) 2020, Wolber Lab.


## Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
