# MSqRob

Robust statistical inference for quantitative LC-MS proteomics.

## IMPORTANT THE PACKAGE IS NO LONGER MAINTAINED ==> [msqrob2](https://www.bioconductor.org/packages/release/bioc/html/msqrob2.html)
All functionalities are ported to our novel bioconductor tool [msqrob2](https://www.bioconductor.org/packages/release/bioc/html/msqrob2.html)

The MSqRob package allows a user to do quantitative protein-level statistical inference on LC-MS proteomics data. More specifically, our package makes use of peptide-level input data, thus correcting for unbalancedness and peptide-specific biases. As previously shown (Goeminne et al. (2015)), this approach is both more sensitive and specific than summarizing peptide-level input to protein-level values. Model estimates are stabilized by ridge regression, empirical Bayes variance estimation and downweighing of outliers. Currently, only label-free proteomics data types are supported.

The MSqRob Shiny App allows for an easy-to-use graphical user interface that requires no programming skills. Interactive plots are outputted and results are automatically exported to Excel.

**The authors kindly ask to make a reference to Goeminne et al. (2016) and Goeminne et al. (2017) when making use of this package in any kind of publication or presentation.**

## 1. MSqRob as a standalone application for Windows

When you are a Windows user who only wants to use the functionality of the MSqRob Shiny App and who doesn't foresee to use R in the near future, you can go for the MSqRob standalone version. Note that this distribution is currently 596 MB in size, as it contains the R framework internally, so if you have R already installed on your computer, you might consider to install MSqRob within the R framework (see below).

### 1.1. Install MSqRob as a standalone application (Windows only!)

Download the MSqRob GUI through this link: https://github.com/statOmics/MSqRob/releases/download/MSqRob_standalone_Win_0.7.7/MSqRob.zip.  

### 1.2. Run the standalone MSqRob Shiny App (graphical user interface)

To start MSqRob, first unzip the zip file and then doubleclick the file "run_MSqRob.vbs" in the "dist" folder.

## 2. MSqRob within the R framework

### 2.1. Installation

If you're not a Windows user, and/or you are already using R and/or RStudio, and/or you also want to have access to MSqRob on the R command line, follow the following instructions to install MSqRob.

First, make sure that R is installed on your computer. Optionally, you can also install RStudio as an integrated development environment (IDE) for R. Instructions on how to install R and RStudio can be found here: http://web.cs.ucla.edu/~gulzar/rstudio/.

On Windows, make sure that RTools is installed. Go to: https://cran.r-project.org/bin/windows/Rtools/ to download RTools. A user guide on how to install RTools on Windows can be found at: https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows. Errors in MSqRob on Windows related to unable to zip the results Excel file might be related to errors in configuring RTools.

To install MSqRob directly from GitHub, we first need to install the package `devtools`. Just enter the following commands in RStudio:

~~~~
install.packages("devtools")
library(devtools)
~~~~

Then, install `Bioconductor`:

~~~~
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()
~~~~

Finally, we call this to install the latest version of `MSqRob` (0.7.7):

~~~~
devtools::install_github("statOmics/MSqRob@MSqRob0.7.7")
library(MSqRob)
~~~~

If you encounter any issues in installing MSqRob, please let us know at https://github.com/statOmics/MSqRob/issues/1!

### 2.2. Run the MSqRob Shiny App (graphical user interface)

Enter the following command in RStudio to run the `MSqRob` Shiny app:

~~~~
shiny::runApp(system.file("App-MSqRob", package="MSqRob"))
~~~~

## 3. Input and output

MSqRob needs only two inputs:

1. A file with peptide-specific intensities.
2. An experimental annotation file (Excel or tab-delimited). A template for this file can be generated within the Shiny App.

In the Shiny App, the file with peptide-specific intensities can be in one of the following formats: MaxQuant (peptides.txt file), moFF (Argentini et al. (2016)) peptide-level output, mzTab format or Progenesis csv. Optionally for MaxQuant, a proteinGroups.txt file is needed as well if proteins that are only identified by peptides carrying a chemical modification need to be removed in the preprocessing stage.

When using the R script, any type of peptide-level input will do, as long as you can convert it to an MSnbase object (Gatto and Lilley (2012)).

## 4. Tutorials for the MSqRob Shiny App

A tutorial on using the MSqRob Shiny App can be found in Goeminne et al. (2017). For the *Francisella* example, one finds the peptides.txt file, the experimental annotation file (label-free_Francisella_annotation.xlsx) and the proteinGroups.txt file under  https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/Francisella/. For the CPTAC example, these files can be found under https://github.com/statOmics/MSqRobData/tree/master/inst/extdata/CPTAC (use label-free_CPTAC_annotation.xlsx as the experimental annotation file).

## 5. Getting started with MSqRob on the command line in R

To get started with MSqRob, we suggest to take a look at the MSqRob vignette at https://github.com/statOmics/MSqRob/blob/master/vignettes/MSqRob.Rmd

If you want to try a clean, simple analysis, make sure to check out our two examples at:
https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/analysis_CPTAC.Rmd
and
https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/Francisella/analysis_Francisella.Rmd

These files can be easily adjusted to cope with the analysis of your own experiment.

## 6. Run MSqRob on the terminal on a Linux system

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

## 7. Run MSqRob via a bash script

An example bash script can be found at:

https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/bash_CPTAC.sh

Note that in this script, the path `/Applications/RStudio.app/Contents/MacOS/pandoc` might need to be changed to the path where your Pandoc installation is saved (see previous paragraph). Likewise `/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/analysis_CPTAC.Rmd` needs to be changed to the folder where your Rmarkdown file is saved. To execute the script, the following two lines need to be passed to the Terminal (where `/Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh` needs to be the file path to where your bash script is saved). The first line of code makes the script executable, the second line executes the script.

~~~~
chmod +x /Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh
bash /Users/lgoeminn/MSqRobData/inst/extdata/CPTAC/bash_CPTAC.sh
~~~~

## 8. Troubleshooting

**If you experience any kind of inconvencience or something is unclear to you, please do not hesitate to contact us at: [ludger.goeminne@epfl.ch](mailto:ludger.goeminne@epfl.ch).** For us it is important to get user-feedback so that we can continue improving MSqRob's user-friendliness. Therefore, we are happy to help you with any problems you might encounter when installing or using MSqRob. So for any issues or questions, even if you think that they might be trivial at first sight, please do send us a mail.

## 9. Extensions of MSqRob

If you are interested in increased power to detect proteins with many missing values, check out our new hurdle method at: https://github.com/statOmics/MSqRobHurdlePaper/ (Goeminne et al. (2020), R script only).

If you are interested in generating protein-level summaries, check out our new modular MSqRob framework for protein-level summarization at https://github.com/statOmics/MSqRobSum (Sticker et al. (2020), R script only).

## 10. References

Ludger Goeminne, Andrea Argentini, Lennart Martens and Lieven Clement (2015). Summarization vs Peptide-Based Models in Label-Free Quantitative Proteomics: Performance, Pitfalls, and Data Analysis Guidelines. Journal of Proteome Research 15(10), pp. 3550-3562.

Ludger Goeminne, Adriaan Sticker, Lennart Martens, Kris Gevaert and Lieven Clement (2020). MSqRob Takes the Missing Hurdle: Uniting Intensity- and Count-Based Proteomics. Analytical Chemistry 92(9), pp. 6278–6287.

Adriaan Sticker, Ludger Goeminne, Lennart Martens and Lieven Clement (2020). Robust Summarization and Inference in Proteome-wide Label-free Quantification. Molecular & Cellular Proteomics 19(7), pp. 1209-1219.

Andrea Argentini, Ludger Goeminne, Kenneth Verheggen, Niels Hulstaert, An Staes, Lieven Clement and Lennart Martens (2016). moFF: a robust and automated approach to extract peptide ion intensities. Nature methods 13, pp. 964–966.

Laurent Gatto and Kathryn Lilley (2012). MSnbase - an R/Bioconductor Package for Isobaric Tagged Mass Spectrometry Data Visualization, Processing and Quantitation. Bioinformatics 28(2), pp. 288–9.

**When using MSqRob, please refer to our published algorithm and our tutorial article on how to use the graphical user interface:**

Ludger Goeminne, Kris Gevaert, and Lieven Clement (2016). Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics 15(2), pp. 657-668.

Ludger Goeminne, Kris Gevaert and Lieven Clement (2018). Experimental design and data-analysis in label-free quantitative LC/MS proteomics: A tutorial with MSqRob. Journal of proteomics, 171(Supplement C), pp. 23-36.
