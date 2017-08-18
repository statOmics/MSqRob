# Change Log
All notable changes to this project will be documented in this file.

## 0.6.5 - 2017-08-18

### Added

 - The function "prot.p.adjust" now gets the option to calculate the values for more than one multiple testing procedure at once. An extra option "fdrtool" is included, which allows FDR correction as specified in Strimmer, K. (2008a & 2008b). An extra element "threshold_excess1" is added to the function. This allows to ignore the excess p-values close to 1 when determining the shape of the null distribution. These excess 1 values are indeed often observed when fitting ridge regression models with sparse data. "threshold_excess1" defaults to 1e-90. In order to get the exact same results from previous MSqRob versions, set "threshold_excess1" to 0 (then a correction will never be performed, no matter how strong the evidence for enrichment in p-values close to 1).
 
 Strimmer, K. (2008a). A unified approach to false discovery rate estimation. BMC Bioinformatics 9: 303. Available from http://www.biomedcentral.com/1471-2105/9/303/.

Strimmer, K. (2008b). fdrtool: a versatile R package for estimating local and tail area- based false discovery rates. Bioinformatics 24: 1461-1462. Available from http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/12/1461.

### Fixed

 - Fixed a bug in which adding no fixed effects to the fit.model function would cause an error.
 - Fixed a bug in which the number of residual degrees of freedom for unfitted empty models appeared incorrect.
 - Various minour stability issues.
 
### Changed

 - Ordinary least squares models now also remove factors from their formulae for which only one level is present, instead of returning an empty model.
 - Added a robust version of M estimation.

## 0.6.4 - 2017-05-28

### Added

- Added question marks to the MSqRob GUI that link the user to a page with extra info on how to use the input terms.
- Added possibility to zoom in on density plots.

### Changed

- Added support for calculating degrees of freedom in 4 different ways: "residual", "Satterthwaite", "exp_between" and "custom". Removed the separate "Satterthwaite" argument from all functions that test contrasts. Instead, "Satterthwaite" is now an option in the new argument "type_dfs".

- Improved preprocessing_wide and preprocessing_long functions by adding default settings and updating their descriptions. Fixed some minor issues that could lead to errors.

### Fixed

- Added an error message when trying to calculate Satterthwaite degrees of freedom with ANOVA, as this is not yet implemented.
- xlim of density plots is now calculated based on the density object and not on the data itself, giving a better overview of the densities.
- Default size of MDS plot is adjusted so that labels are always in the plot window.
- Fixed a bug in which the ANOVA option would wrongly give NA numerator degrees of freedom and p values when the rank of the contrast matrix would be 1. Also, the estimate of the average effect would sometimes be estimated wrongly when the contrast matrix was of less than full rank.

## 0.6.3 - 2017-04-14

- .adjustNames now also works for models with only one random effect.

### Added
- This CHANGELOG file.

### Fixed
- .adjustNames now also works for models with only one random effect.
