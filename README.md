# ReconstPoFD

This repository contains the R-package `ReconstPoFD`. The only function contained in this package is `reconst_fun()` which allows to reconstruct the (centered) missing parts of a function given the (centered) observed part by using the information from the overall covariance function.

The underlying method is described in:

On the Optimal Reconstruction of Partially Observed Functional Data
by Alois Kneip, Dominik Liebl (arXiv:1710.10099)

## Installing the R-package

`library(devtools)`

`install_github("lidom/ReconstPoFD/ReconstPoFD")`

## Replication of Simulation Results

The R-script `MC_simulation.R` allows you to reconstruct the simulation study in:

On the Optimal Reconstruction of Partially Observed Functional Data
by Alois Kneip, Dominik Liebl (arXiv:1710.10099)

The R-script also contains the AIC-based selection of the truncation parameter K.
