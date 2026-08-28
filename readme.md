# Introduction

This repository contains code used in the paper *Worst-Case Probability Bounds for Finite-Horizon Safety under Moment Uncertainty*. It provides details on the implementation and a short user guide to reproduce the results of the paper.

**Keywords**: Safety, Moment-sum-of-squares, Worst-case under uncertainty

## Installation

The following tools must be installed:

1. Download [CasADi v3.6.7](https://web.casadi.org/get/) and add it to your MATLAB path.  
2. Download and install [MOSEK v11.1.x](https://www.mosek.com/downloads/) and add the solver to the MATLAB path.  
3. Download [CaΣoS](https://github.com/casos-opti/casos) and add it to your MATLAB path.

Add the root folder (the one that contains the directory `+casos`) to your MATLAB path.

### Folder Structure

```text
.
├── data/            # Saved results
├── experiments/     # Scripts for the experiments
├── figures/         # Figures and .tex files
├── scripts/         # Auxiliary scripts
├── utils/           # Helper functions 
├── LICENSE
└── readme.md
```

### Remarks
In case of problems, questions or remarks, please contact the corresponding authors (see below). 
- Renato Loureiro: renato.loureiro@ifr.uni-stuttgart.de (main)
- Torbjørn Cunis: torbjoern.cunis@ifr.uni-stuttgart.de


### Citation
Please cite the paper as
```
@misc{loureiro2026,
      title={Worst-Case Probability Bounds for Finite-Horizon Safety under Moment Uncertainty}, 
      author={Renato Loureiro and Torbjørn Cunis},
      year={2026},
      eprint={2608.20121},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2608.20121}, 
}
```

The repository material can be cited with 
```
@data{DARUS-6410_2026,
author = {Loureiro, Renato and Cunis, Torbjørn},
publisher = {DaRUS},
title = {{Source code and Numerical Examples for Worst-Case Probability Bounds for Finite-Horizon Safety under Moment Uncertainty}},
year = {2026},
version = {DRAFT VERSION},
doi = {10.18419/DARUS-6410},
url = {https://doi.org/10.18419/DARUS-6410}
}

```
