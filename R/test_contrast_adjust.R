#' Test a contrast and perform FDR correction
#'
#' @description A shortcut version for contrast testing and multiple testing correction using the \code{\link{test.protLMcontrast}}, \code{\link{prot.p.adjust}} and \code{\link{prot.signif}} functions sequentially.
#' @param protLM An object of class \code{\link[=protLM-class]{protLM}}.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param level	The significance level at which the q-value needs to be controlled. Defaults to 5\%.
#' @param method Correction method. Can be abbreviated. Defaults to "fdr". To get all available methods, type \code{p.adjust.methods}. For more information on these methods, see the \code{\link{p.adjust}} function.
#' @param add.annotations A logical indicating whether the \code{annotations} slot of the \code{\link[=protLM-class]{protLM}} object should be added as extra columns to each matrix in the returned list of matrices. Defaults to \code{TRUE}.
#' @param simplify A logical indicating wheter, if there is only one contrast, a matrix should be returned instead of a list containing one matrix. Defaults to \code{TRUE}.
#' @param lfc The minimum (log2) fold-change that is considered scientifically meaningful. Defaults to \code{0}. Ignored when \code{anova = TRUE}.
#' @param anova A logical indicating whether the contrasts should be tested simultaneously in one F-test (\code{anova = TRUE}) or as separate t-tests (\code{anova = FALSE}). Defaults to \code{FALSE}.
#' @param anova.na.ignore A logical indicating whether contrasts that cannot be fitted due to missing observations should be ignored when calculating an F value (\code{anova.na.ignore = TRUE}) or NA should be returned when at least one contrast cannot be estimated (\code{anova.na.ignore = FALSE}). Defaults to \code{TRUE}. Ignored when \code{anova = FALSE}.
#' @param type_dfs Either one of \code{"residual"}, \code{"between-within"}, \code{"Satterthwaite"}, \code{"exp_between"} or \code{"custom"}. This argument indicates how the degrees of freedom should be calculated. Defaults to \code{"residual"}.
#' More information is given under 'Details'.
#' @param custom_dfs Only used if \code{type_dfs="custom"}. A list of length equal to the number of models in \code{protLM} containing vectors of lenght \code{ncol(L)} (if \code{anova = FALSE}) representing the degrees of freedom that should be used for each contrast in \code{L} and each model in the \code{protLM} object. If \code{anova = TRUE}, each element of the list should contain a single numeric value representing the degrees of freedom to be used for the anova test. Defaults to \code{NULL}.
#' @param exp_unit Only used if \code{type_dfs="exp_between"}. The effect that in all models corresponds to the experimental unit.
#' @param pars_between Only used if \code{type_dfs="exp_between"}. Character vector indicating all parameters in the models that are between-treatment effects. If left to default (\code{NULL}), all parameters in the models will be asumed to be between-treatment effects (this is not adviced as the result will mostly be too conservative).
#' @param lmerModFun Only used when \code{satterthwaite=TRUE}. \code{lmerModFun} indicates which deviance function should be used when calculating the Satterthwaite approximation for the degrees of freedom. The default (\code{NULL}) uses the lme4 \code{\link[=lme4::mkLmerDevfun]{mkLmerDevfun}} function to generate the deviance function. This parameter should only rarely, if ever, be changed.
#' @param gradMethod Only used when \code{satterthwaite=TRUE}. One of "Richardson", "simple", or "complex" indicating the method to use for the gradient calculation by numerical approximation during the calculation of the Satterthwaite approximation for the degrees of freedom. Defaults to "simple".
#' @param printProgress A logical indicating whether the R should print a message before calculating the contrasts for each accession. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message_extract Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the extraction of beta, vcov, df and sigma, or \code{NULL} to hide the current message (if any).
#' @param message_test Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the testing of the contrasts, or \code{NULL} to hide the current message (if any).
#' @details Calculating degrees of freedom (and hence p values) for mixed models with unbalanced designs is an unresolved issue in the field (see for example here https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html and here https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q2/000904.html).
#' We offer different approximations and leave it up to the user to select his/her preferred approach.
#' \code{"residual"} calculates approximative degrees of freedom by subtracting the trace of the hat matrix from the number of observations. It is the default setting, but this approach might be somewhat too liberal.
#' \code{"Satterthwaite"} calculates approximative degrees of freedom using Satterthwaite's approximation (Satterthwaite, 1946). This approximative approach is used in many applications but is rather slow to calculate and might lead to some missing values due difficulties in calculating the Hessian.
#' \code{"exp_between"} calculates approximative degrees of freedom by defining on which level the treatments were executed and substracting all degrees of freedom lost due to between-treatement effects (\code{pars_between}) from the number of experimental units \code{exp_unit}. This allows to mimick the behaviour of \code{type_dfs="between-within"} for more complex designs.
#' \code{"custom"} Allows the user to provide his/her own degrees of freedom for each contrast and each protein. Custom degrees of freedom should be entered in the \code{custom_dfs} field.
#' @return A list of matrices, with each matrix in the list corresponding to a contrast in \code{L}. Each row of the matrix corresponds to an accession in the \code{\link[=protLM-class]{protLM}} object.
#' The \code{estimate} column contains the size estimate of the contrast, the \code{se} column contains the estimated standard error on the contrast, the \code{Tval} column contains the T-value corresponding to the contrast, the \code{pval} column holds the p-value corresponding to the contrast and the \code{qval} column holds the corrected p-values.
#' Each matrix is sorted from smalles to largest \code{pval} with \code{NA} values at the bottom of the matrices.
#' If \code{simplify=TRUE} and the \code{\link[=protLM-class]{protLM}} object contains only one element, the matrix is not present in a list.
#' @references F. E. Satterthwaite.
#' An Approximate Distribution of Estimates of Variance Components.
#' Biometrics Bulletin, 1946.
#' @include testProtLMContrast.R
#' @include updateProgress.R
#' @include prot_p_adjust.R
#' @include prot_signif.R
#' @export
test.contrast_adjust <- function(protLM, L, level = 0.05, method = "fdr", add.annotations = TRUE, simplify = TRUE, lfc = 0, anova = FALSE, anova.na.ignore = TRUE, type_dfs = "residual", custom_dfs = NULL, exp_unit = NULL, pars_between = NULL, lmerModFun = NULL, gradMethod = "simple", printProgress = FALSE, shiny = FALSE, message_extract = NULL, message_test = NULL)
{
  contrasts <- test.protLMcontrast(protLM, L, add.annotations=add.annotations, simplify=simplify, lfc=lfc, anova=anova, anova.na.ignore=anova.na.ignore, type_dfs=type_dfs, custom_dfs=custom_dfs, exp_unit=exp_unit, pars_between=pars_between, lmerModFun = lmerModFun, gradMethod = gradMethod, printProgress=printProgress, shiny=shiny, message_extract=message_extract, message_test=message_test)
  contrasts <- prot.p.adjust(contrasts, method=method)
  contrasts <- prot.signif(contrasts, level=level)
  return(contrasts)
}
