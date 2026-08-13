# Introduction

This repository contains code used in the paper *Worst-Case Probability Bounds for Finite-Horizon Safety under Moment Uncertainty*. It provides details on the implementation and a short user guide to reproduce the results of the paper.

**Keywords**: Safety, Moment-sum-of-squares, Worst-case under uncertainty

## Installation

The following tools must be installed:

1. Download [CasADi v3.6.7](https://web.casadi.org/get/) and add it to your MATLAB path.  
2. Download and install [MOSEK v11.1.x](https://www.mosek.com/downloads/) and add the solver to the MATLAB path.  
3. Get CaΣoS from **toDO**

Add the root folder (the one that contains the directory `+casos`) to your MATLAB path.

### Folder Structure

```text
.
│   ### Directories ###
│ 
├── fig/                        # Contains the figures used in the paper 
├── results/                    # Saved .mat files
├── utils/                      # Contains auxiliary functions
│ 
│   ### Scripts ###
│   # Example 1
├── ex01_1_meas.m               # ToDO
├── ex01_1_sos_ambiguity.m      # ToDO
├── ex01_1_sos.m                # ToDO
│ 
│   # Example 2
├── ex02_meas.m                 # ToDO
├── ex02_meas_moving.m          # ToDO
├── ex02_sos_ambiguity.m        # ToDO
├── ex02_sos.m                  # ToDO
├── ex02_sos_moving.m           # ToDO
├── ex02_sos_multiple.m         # ToDO
│ 
├── LICENSE
└── readme.md
```

### Remarks
In case of problems, questions or remarks, please contact the corresponding authors (see below). 
- Renato Loureiro: renato.loureiro@ifr.uni-stuttgart.de (main)
- Torbjørn Cunis: torbjoern.cunis@ifr.uni-stuttgart.de


### Citation
Please cite the paper as PREPRINT FOLLOWS SOON
```

```

The repository material can be cited with (SOON)
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