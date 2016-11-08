# MSqRob

Robust statistical inference for quantitative LC-MS proteomics.

## Installation

On Windows, make sure that RTools is installed. Go to: https://cran.r-project.org/bin/windows/Rtools/ to download RTools. A user guide on how to install RTools on Windows can be found at: https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows

First, we need to install the package `devtools`:

~~~~
install.packages("devtools")
library(devtools)
~~~~

Then, we call this to install `MSqRob`:

~~~~
devtools::install_github("ludgergoeminne/MSqRob")
library(MSqRob)
~~~~

## Run the Shiny App

Just enter the following command to run the `MSqRob` Shiny app:

~~~~
shiny::runApp(system.file("App-MSqRob", package="MSqRob"))
~~~~

## References

When using MSqRob, please refer to our published algorithm:

Goeminne, L. J. E., Gevaert, K., and Clement, L. (2016) Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics 15(2), pp 657-668.
