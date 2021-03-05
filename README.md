Dynophores
==========

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/dominiquesydow/dynophores/workflows/CI/badge.svg)](https://github.com/dominiquesydow/dynophores/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/dominiquesydow/dynophores/branch/master/graph/badge.svg)](https://codecov.io/gh/dominiquesydow/dynophores/branch/master)
[![Documentation Status](https://readthedocs.org/projects/dynophores/badge/?version=latest)](https://dynophores.readthedocs.io/en/latest/?badge=latest)

> ⚠ This project is work-in-progress. The API is not final.


![OpenCADD](/docs/_static/dynophore.png)

Dynamic pharmacophore modeling of molecular interactions.

Dynophores are an intuitive tool to explore the dynamics of molecular interactions between 
a macromolecule and a ligand throughout an MD simulation: 

A __dynophore__ is a collection of so-called __superfeatures__. 
A superfeature is defined as a pharmacophore feature on ligand site 
(defined by a feature type, e.g. HBA, and one or more ligand atom numbers/serials) that occurs 
at least once during an MD simulation. A superfeature can have one or more interaction partner(s) 
on macromolecule-side. 
These interaction partners are called __environmental partners__. 
Superfeatures and their environmental partners are monitored throughout an MD simulation.

## Documentation

The documentation will be available [here](https://dynophores.readthedocs.io/en/latest/).

## Overview

This Python library offers visualization schemes for dynophore data generated with the `DynophoreApp` in a single Jupyter Notebook.  

- 3D view of the dynophore's point clouds (one cloud per superfeature) using `nglview` allows easy 
visual inspection of the dynamic macromolecule-ligand interactions. Point clouds are rendered 
alongside the topology and (optionally) the trajectory underlying the dynophore.
- Statistics cover the occurrence of superfeatures and their environmental partners as well as distances between them.
- Dynophore data can be further analyzed conveniently right there in the same notebook by working with the `Dynophore` class.


## License

`dynophores` is a free software and is licensed under the MIT license. 
Copyright (c) 2020, Dominique Sydow.


## Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
