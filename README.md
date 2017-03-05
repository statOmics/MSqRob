# MSqRob

Robust statistical inference for quantitative LC-MS proteomics.

The MSqRob package allows a user to do quantitative protein-level statistical inference on LC-MS proteomics data. More specifically, our package makes use of peptide-level input data, thus correcting for unbalancedness and peptide-specific biases. As previously shown (Goeminne et al. (2015)), this approach is both more sensitive and specific than summarizing peptide-level input to protein-level values. Model estimates are stabilized by ridge regression, empirical Bayes variance estimation and downweighing of outliers. Currently, only label-free proteomics data types are supported.

The MSqRob Shiny App allows for an easy-to-use graphical user interface that requires no programming skills. Interactive plots are outputted and results are automatically exported to Excel.

The authors kindly ask to make a reference to (Goeminne et al. (2016)) when making use of this package in any kind of publication or presentation.

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

## Run the MSqRob Shiny App (graphical user interface)

Just enter the following command in RStudio to run the `MSqRob` Shiny app:

~~~~
shiny::runApp(system.file("App-MSqRob", package="MSqRob"))
~~~~

## Getting started with MSqRob in R

To get started with MSqRob, we suggest to take a look at the MSqRob vignette at https://github.com/statOmics/MSqRob/blob/master/vignettes/MSqRob.Rmd

If you want to try a clean, simple analysis, make sure to check out our two examples at:
https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/analysis_CPTAC.Rmd
and
https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/Francisella/analysis_Francisella.Rmd

These files can be easily adjusted to cope with the analysis of your own experiment.

## Run MSqRob on the terminal on a Linux system

In order to run the first example in the terminal on a Linux system, it is first needed to tell the system where to find the Pandoc installation. This folder can easily be found by typing the following command in RStudio:

~~~~
Sys.getenv("RSTUDIO_PANDOC")
~~~~

In our case, this folder is `/Applications/RStudio.app/Contents/MacOS/pandoc`.
We thus execute the following in the Terminal:
~~~~
export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc
~~~~

Then we knit out Rmarkdown file (which we saved at `/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/analysis_CPTAC.Rmd`) by running the following command in the Terminal:

~~~~
Rscript -e "require( 'rmarkdown' ); render('/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/analysis_CPTAC.Rmd', 'html_document')"
~~~~

## Run MSqRob via a bash script

An example bash script can be found at:

https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/bash_CPTAC.sh

Note that in this script, the path `/Applications/RStudio.app/Contents/MacOS/pandoc` might need to be changed to the path where your Pandoc installation is saved (see previous paragraph). Likewise `/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/analysis_CPTAC.Rmd` needs to be changed to the folder where your Rmarkdown file is saved. To execute the script, the following two lines need to be passed to the Terminal (where `/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh` needs to be the file path to where your bash script is saved). The first line of code makes the script executable, the second line executes the script.

~~~~
chmod +x /Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh
bash /Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh
~~~~

## Contact

We are happy to help you with any problems you might encounter when using MSqRob.
For any issues or questions, please e-mail us at: ludger.goeminne@vib-ugent.be

## References

Ludger Goeminne, Andrea Argentini, Lennart Martens and Lieven Clement (2015). Summarization vs Peptide-Based Models in Label-Free Quantitative Proteomics: Performance, Pitfalls, and Data Analysis Guidelines. Journal of Proteome Research 15(10), 3550-3562.

When using MSqRob, please refer to our published algorithm:

Goeminne, L. J. E., Gevaert, K., and Clement, L. (2016) Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics 15(2), pp 657-668.
