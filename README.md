# MSqRob

Robust statistical inference for quantitative LC-MS proteomics.

The MSqRob package allows a user to do quantitative protein-level statistical inference on LC-MS proteomics data. More specifically, our package makes use of peptide-level input data, thus correcting for unbalancedness and peptide-specific biases. As previously shown (Goeminne et al. (2015)), this approach is both more sensitive and specific than summarizing peptide-level input to protein-level values. Model estimates are stabilized by ridge regression, empirical Bayes variance estimation and downweighing of outliers. Currently, only label-free proteomics data types are supported. The authors kindly ask to make a reference to (Goeminne et al. (2016)) when making use of this package in any kind of publication or presentation.

## Installation

On Windows, make sure that RTools is installed. Go to: https://cran.r-project.org/bin/windows/Rtools/ to download RTools. A user guide on how to install RTools on Windows can be found at: https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows. Errors in MSqRob on Windows related to unable to zip the results Excel file might be related to errors in configuring RTools.

To install MSqRob directly from GitHub, we first need to install the package `devtools`:

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

Ludger Goeminne, Andrea Argentini, Lennart Martens and Lieven Clement (2015). Summarization vs Peptide-Based Models in Label-Free Quantitative Proteomics: Performance, Pitfalls, and Data Analysis Guidelines. Journal of Proteome Research 15(10), 3550-3562.

When using MSqRob, please refer to our published algorithm:

Goeminne, L. J. E., Gevaert, K., and Clement, L. (2016) Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics 15(2), pp 657-668.
