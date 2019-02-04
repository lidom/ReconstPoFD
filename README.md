# R-package: `ReconstPoFD`

This repository contains the R-package `ReconstPoFD`. The underlying statistical methods are described in:  
[On the Optimal Reconstruction of Partially Observed Functional Data](https://arxiv.org/abs/1710.10099)  
by Alois Kneip and Dominik Liebl (arXiv:1710.10099)

## Installing the R-package

`library(devtools)`

`install_github("lidom/ReconstPoFD/ReconstPoFD")`


## Using the Package

`library(ReconstPoFD)`


The main function containing our different reconstruction methods is:
`reconstructKneipLiebl()`

In order to get more information, check out its help file:
`?reconstructKneipLiebl`


The file `Simulation.R` allows to reproduce the simulation study in Section 6 of our manuscript.
