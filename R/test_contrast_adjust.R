#' Test a contrast and perform FDR correction
#'
#' @description A shortcut version for contrast testing and multiple testing correction using the \code{\link{test.protLMcontrast}}, \code{\link{prot.p.adjust}} and \code{\link{prot.signif}} functions sequentially.
#' @param protLM An object of class \code{\link[=protLM-class]{protLM}}.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param level	The significance level at which the q-value needs to be controlled. Defaults to 5\%.
#' @param method Correction method. Can be abbreviated. Defaults to "fdr". To get all available methods, type \code{p.adjust.methods}. For more information on these methods, see the \code{\link{p.adjust}} function.
#' @param add.annotations A logical indicating whether the \code{annotations} slot of the \code{\link[=protLM-class]{protLM}} object should be added as extra columns to each matrix in the returned list of matrices. Defaults to \code{TRUE}.
#' @param squeezeVar A logical indicating whether residual standard deviations should be squeezed together by computing empirical Bayes posterior means. For more information, see limma's \code{\link{squeezeVar}} function. Defaults to \code{TRUE}.
#' @param par_squeeze Character vector indicating which model parameters need to be squeezed. When squeezing random effects, provide their names. Fixed effects are present in shrinkage groups, e.g. ridgeGroup.1. If you want them to be squeezed as well, provide the names of the shrinkage groups that need to be squeezed. The default \code{NULL} indicates that no parameters will be squeezed.
#' @param min_df A numeric value indicating the minimal number of residual degrees of freedom that is required for a protein to be taken into account for the squeezing of the variances. Only used when \code{squeezeVar} is set to \code{TRUE}. Defaults to \code{1}.
#' @param robust_var A logical indicating wheter the squeezing of the variances should be robust to the presence of outlier sample variances. Defaults to \code{TRUE}.
#' @param simplify A logical indicating wheter, if there is only one contrast, a matrix should be returned instead of a list containing one matrix. Defaults to \code{TRUE}.
#' @param lfc The minimum (log2) fold-change that is considered scientifically meaningful. Defaults to \code{0}. Ignored when \code{anova = TRUE}.
#' @param anova A logical indicating whether the contrasts should be tested simultaneously in one F-test (\code{anova = TRUE}) or as separate t-tests (\code{anova = FALSE}). Defaults to \code{FALSE}.
#' @param anova.na.ignore A logical indicating whether contrasts that cannot be fitted due to missing observations should be ignored when calculating an F value (\code{anova.na.ignore = TRUE}) or NA should be returned when at least one contrast cannot be estimated (\code{anova.na.ignore = FALSE}). Defaults to \code{TRUE}. Ignored when \code{anova = FALSE}.
#' @param exp_unit The effect that in all models corresponds to the experimental unit. Only needed when one would like to calculate a more conservative way of estimating the degrees of freedom.
#' The default way of estimating the degrees of freedom (\code{exp_unit=NULL}) subtracts the total number of observations by the trace of the Hat matrix. However, often, observations are not completely independent. A more conservative way (\code{df_exp}) is defining on which level the treatments were executed and substracting all degrees of freedom lost due to between-treatement effects (\code{pars_df}) from the number of treatments.
#' @param pars_df Only used if exp_unit is not \code{NULL}. Character vector indicating all parameters in the models that are between-treatment effects in order to calculate a more conservative degrees of freedom (\code{df_exp}). If left to default (\code{NULL}), all parameters in the models will be asumed to be between-treatment effects (this is not adviced as the result will mostly be too conservative).
#' @param satterthwaite A logical indicating whether the Satterthwaite approximation for the degrees of freedom should be used instead of the classical trace of the Hat matrix. Defaults to \code{FALSE}.
#' @param ... Additional arguments to be passed to the squeezeVarRob function internally.
#' @return A list of matrices, with each matrix in the list corresponding to a contrast in L. Each row of the matrix corresponds to an accession in the \code{\link[=protLM-class]{protLM}} object.
#' The \code{estimate} column contains the size estimate of the contrast, the \code{se} column contains the estimated standard error on the contrast, the \code{Tval} column contains the T-value corresponding to the contrast, the \code{pval} column holds the p-value corresponding to the contrast and the \code{qval} column holds the corrected p-values.
#' Each matrix is sorted from smalles to largest \code{pval} with \code{NA} values at the bottom of the matrices.
#' If \code{simplify=TRUE} and the \code{\link[=protLM-class]{protLM}} object contains only one element, the matrix is not present in a list.
#' @include testProtLMContrast.R
#' @include prot_p_adjust.R
#' @include prot_signif.R
#' @export
test.contrast_adjust <- function(protLM, L, level=0.05, method="fdr", add.annotations=TRUE, squeezeVar=TRUE, par_squeeze=NULL, min_df=1, custom_dfs=NULL, robust_var=TRUE, simplify=TRUE, lfc=0, anova=FALSE, anova.na.ignore=TRUE, exp_unit=NULL, pars_df=NULL, satterthwaite=FALSE, ...)
{
  contrasts <- test.protLMcontrast(protLM, L, add.annotations=add.annotations, squeezeVar=squeezeVar, par_squeeze=par_squeeze, min_df=min_df, custom_dfs=custom_dfs, robust_var=robust_var, simplify=simplify, lfc=lfc, anova=anova, anova.na.ignore=anova.na.ignore, exp_unit=exp_unit, pars_df=pars_df, satterthwaite=satterthwaite, ...)
  contrasts <- prot.p.adjust(contrasts, method=method)
  contrasts <- prot.signif(contrasts, level=level)
  return(contrasts)
}
