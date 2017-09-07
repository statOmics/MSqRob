# MSqRob

Robust statistical inference for quantitative LC-MS proteomics.

The MSqRob package allows a user to do quantitative protein-level statistical inference on LC-MS proteomics data. More specifically, our package makes use of peptide-level input data, thus correcting for unbalancedness and peptide-specific biases. As previously shown (Goeminne et al. (2015)), this approach is both more sensitive and specific than summarizing peptide-level input to protein-level values. Model estimates are stabilized by ridge regression, empirical Bayes variance estimation and downweighing of outliers. Currently, only label-free proteomics data types are supported.

The MSqRob Shiny App allows for an easy-to-use graphical user interface that requires no programming skills. Interactive plots are outputted and results are automatically exported to Excel.

The authors kindly ask to make a reference to Goeminne et al. (2016) and Goeminne et al. (2017) when making use of this package in any kind of publication or presentation.

## 1. MSqRob as a standalone for Windows

When you are a Windows user who only wants to use the functionality of the MSqRob Shiny GUI and who doesn't foresee to use R in the near future, you can go for the MSqRob standalone version. Note that this distribution is around 500 MB in size, as it contains the R framework internally, so if you have R already installed on your computer, you might consider to install MSqRob within the R framework (see below).

### 1.1. Installation

Download the MSqRob GUI through this link: https://github.com/statOmics/MSqRob/releases/download/MSqRob_standalone_Win_0.7.0/MSqRob.zip.  

### 1.2. Run the MSqRob Shiny App (graphical user interface)

To start MSqRob, doubleclick the file "run_MSqRob.vbs" in the "dist" folder.

## 2. MSqRob within the R framework

### 2.1. Installation

If you're not a Windows user, and/or you are already using R and/or RStudio, and/or you want to use MSqRob on the R command line, follow the following instructions to install MSqRob.

First, make sure that R is installed on your computer. Optionally, you can also install RStudio as an integrated development invironment (IDE) for R. Instructions on how to install R and RStudio can be found here: http://web.cs.ucla.edu/~gulzar/rstudio/.

On Windows, make sure that RTools is installed. Go to: https://cran.r-project.org/bin/windows/Rtools/ to download RTools. A user guide on how to install RTools on Windows can be found at: https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows. Errors in MSqRob on Windows related to unable to zip the results Excel file might be related to errors in configuring RTools.

To install MSqRob directly from GitHub, we first need to install the package `devtools`. Just enter the following commands in RStudio:

~~~~
install.packages("devtools")
library(devtools)
~~~~

Then, install `Bioconductor`:

~~~~
source("https://bioconductor.org/biocLite.R")
biocLite()
~~~~

Finally, we call this to install the latest version of `MSqRob` (0.7.0):

~~~~

devtools::install_github("statOmics/MSqRob@MSqRob0.7.0")

library(MSqRob)
~~~~

If you encounter any issues in installing MSqRob, please let us know at https://github.com/statOmics/MSqRob/issues/1!

### 2.2. Run the MSqRob Shiny App (graphical user interface)

Enter the following command in RStudio to run the `MSqRob` Shiny app:

~~~~
shiny::runApp(system.file("App-MSqRob", package="MSqRob"))
~~~~

## 3. Tutorials for the MSqRob Shiny App

A tutorial on using the MSqRob Shiny App can be found in Goeminne et al. (2017). For the *Francisella* example, one finds the peptides.txt file, the experimental annotation file (label-free_Francisella_annotation.xlsx) and the proteinGroups.txt file under  https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/Francisella/. For the CPTAC example, these files can be found under https://github.com/statOmics/MSqRobData/tree/master/inst/extdata/CPTAC (use label-free_CPTAC_annotation.xlsx as the experimental annotation file).

## 4. Getting started with MSqRob on the command line in R

To get started with MSqRob, we suggest to take a look at the MSqRob vignette at https://github.com/statOmics/MSqRob/blob/master/vignettes/MSqRob.Rmd

If you want to try a clean, simple analysis, make sure to check out our two examples at:
https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/analysis_CPTAC.Rmd
and
https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/Francisella/analysis_Francisella.Rmd

These files can be easily adjusted to cope with the analysis of your own experiment.

## 5. Run MSqRob on the terminal on a Linux system

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

## 6. Run MSqRob via a bash script

An example bash script can be found at:

https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/bash_CPTAC.sh

Note that in this script, the path `/Applications/RStudio.app/Contents/MacOS/pandoc` might need to be changed to the path where your Pandoc installation is saved (see previous paragraph). Likewise `/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/analysis_CPTAC.Rmd` needs to be changed to the folder where your Rmarkdown file is saved. To execute the script, the following two lines need to be passed to the Terminal (where `/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh` needs to be the file path to where your bash script is saved). The first line of code makes the script executable, the second line executes the script.

~~~~
chmod +x /Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh
bash /Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh
~~~~

## 7. Troubleshooting

**If you experience any kind of inconvencience or something is unclear to you, please do not hesitate to contact us at: [ludger.goeminne@vib-ugent.be](mailto:ludger.goeminne@vib-ugent.be).** For us it is important to get user-feedback so that we can continue improving MSqRob's user-friendliness. Therefore, we are happy to help you with any problems you might encounter when installing or using MSqRob. So for any issues or questions, even if you think that they might be trivial at first sight, please do send us a mail.

## 8. References

Ludger Goeminne, Andrea Argentini, Lennart Martens and Lieven Clement (2015). Summarization vs Peptide-Based Models in Label-Free Quantitative Proteomics: Performance, Pitfalls, and Data Analysis Guidelines. Journal of Proteome Research 15(10), 3550-3562.

**When using MSqRob, please refer to our published algorithm and our tutorial article on how to use the graphical user interface:**

Ludger Goeminne, Kris Gevaert, and Lieven Clement (2016). Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics 15(2), pp 657-668.

Ludger Goeminne, Kris Gevaert and Lieven Clement (2017). Experimental design and data-analysis in label-free quantitative LC/MS proteomics: A tutorial with MSqRob. Journal of Proteomics (in press).
