# Change Log
All notable changes to this project will be documented in this file.

## 0.7.0 - 2017-05-28 [in progress]

### Changed

 - The function "prot.p.adjust" now gets the option to calculate the values for more than one multiple testing procedure at once. An extra option "fdrtool" is included, which allows FDR correction as specified in Strimmer, K. (2008a & 2008b). Extra elements "mitigate_excess1" and "threshold_excess1" are added to the function. These allow to ignore the excess p values close to 1 when determining the shape of the null distribution. These excess 1 values are indeed often observed when fitting ridge regression models with sparse data.
 
 Strimmer, K. (2008a). A unified approach to false discovery rate estimation. BMC Bioinformatics 9: 303. Available from http://www.biomedcentral.com/1471-2105/9/303/.

Strimmer, K. (2008b). fdrtool: a versatile R package for estimating local and tail area- based false discovery rates. Bioinformatics 24: 1461-1462. Available from http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/12/1461.

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
