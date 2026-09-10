# GR-KG: Massive Scalar Fields and Hyperboloidal Foliations

This repository is a modified copy of the original *GR-EM-KG* code developed by [Dr. Alex Vañó-Viñuales](https://github.com/alexvanov). It was adopted and used during an internship completed by [Javier Riera](https://github.com/javierriera8) during the summer of 2026, at IAC3 (Institute of Applied Computing and Community Code, Universitat de les Illes Balears).

The code simulates the Einstein-Klein-Gordon equations in spherical symmetry using hyperboloidal foliations. It's specifically designed to study the stability of massive scalar fields and mitigate the numerical divergences introduced by the massive term. Note that while the original code implemented Maxwell's equations on electromagnetism, that functionality in not included in this version. 

The contained material is the following:
- **Fortran source code**: Implementation of the evolution equations.
- **Makefile**: Used for compiling the code.
- **Parameter template**: Configuration file for running the simulations.
- **Python Notebook**: A Jupyter Notebook designed to visualize specific variables and analyze convergence tests.

