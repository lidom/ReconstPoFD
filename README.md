# ReconstPoFD

This repository contains the R-package `ReconstPoFD`. The only function contained in this package is `reconst_fun()` which allows to reconstruct the (centered) missing parts of a function given the (centered) observed part. 

The underlying statistical method is described in:  
[On the Optimal Reconstruction of Partially Observed Functional Data](https://arxiv.org/abs/1710.10099)  
by Alois Kneip and Dominik Liebl (arXiv:1710.10099)

## Installing the R-package

`library(devtools)`

`install_github("lidom/ReconstPoFD/ReconstPoFD")`

## Replication of Simulation Results

The R-script `MC_simulation.R` allows you to reconstruct the simulation study in:  
[On the Optimal Reconstruction of Partially Observed Functional Data](https://arxiv.org/abs/1710.10099)  
by Alois Kneip and Dominik Liebl (arXiv:1710.10099)

Note: The R-script also contains the AIC-based selection of the truncation parameter K.
