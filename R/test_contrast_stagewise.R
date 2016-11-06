#' Perform stage-wise contrast testing
#'
#' @description A shortcut version for contrast testing and multiple testing correction using first \code{\link{test.contrast_adjust}} to perform an ANOVA as a first stage.
#' In the second stage, the functions \code{\link{test.protLMcontrast}}, \code{\link{prot.p.adjust_protwise}} and \code{\link{prot.signif}} are used sequentially only on the proteins retained in the first stage.
#' This approach is more powerful than the classic \code{\link{test.contrast_adjust}} when testing for contrasts where related proteins are expected to be differentially abundant. This approach is also better for identifying interactions.
#' @param protLM An object of class \code{\link[=protLM-class]{protLM}}.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param level	The significance level at which the q-value needs to be controlled. Defaults to 5\%.
#' @param method_stage1 Correction method to be used in the first stage ANOVA. Can be abbreviated. Defaults to "fdr". To get all available methods, type \code{p.adjust.methods}. For more information on these methods, see the \code{\link{p.adjust}} function.
#' @param add.annotations A logical indicating whether the \code{annotations} slot of the \code{\link[=protLM-class]{protLM}} object should be added as extra columns to each matrix in the returned list of matrices. Defaults to \code{TRUE}.
#' @param squeezeVar A logical indicating whether residual standard deviations should be squeezed together by computing empirical Bayes posterior means. For more information, see limma's \code{\link{squeezeVar}} function. Defaults to \code{TRUE}.
#' @param par_squeeze Character vector indicating which model parameters need to be squeezed. When squeezing random effects, provide their names. Fixed effects are present in shrinkage groups, e.g. ridgeGroup.1. If you want them to be squeezed as well, provide the names of the shrinkage groups that need to be squeezed. The default \code{NULL} indicates that no parameters will be squeezed.
#' @param min_df A numeric value indicating the minimal number of residual degrees of freedom that is required for a protein to be taken into account for the squeezing of the variances. Only used when \code{squeezeVar} is set to \code{TRUE}. Defaults to \code{1}.
#' @param robust_var A logical indicating wheter the squeezing of the variances should be robust to the presence of outlier sample variances. Defaults to \code{TRUE}.
#' @param simplify A logical indicating wheter, if there is only one contrast, a matrix should be returned instead of a list containing one matrix. Defaults to \code{TRUE}.
#' @param lfc The minimum (log2) fold-change that is considered scientifically meaningful. Defaults to \code{0}. Ignored when \code{anova = TRUE}.
#' @param exp_unit The effect that in all models corresponds to the experimental unit. Only needed when one would like to calculate a more conservative way of estimating the degrees of freedom.
#' The default way of estimating the degrees of freedom (\code{exp_unit=NULL}) subtracts the total number of observations by the trace of the Hat matrix. However, often, observations are not completely independent. A more conservative way (\code{df_exp}) is defining on which level the treatments were executed and substracting all degrees of freedom lost due to between-treatement effects (\code{pars_df}) from the number of treatments.
#' @param pars_df Only used if exp_unit is not \code{NULL}. Character vector indicating all parameters in the models that are between-treatment effects in order to calculate a more conservative degrees of freedom (\code{df_exp}). If left to default (\code{NULL}), all parameters in the models will be asumed to be between-treatment effects (this is not adviced as the result will mostly be too conservative).
#' @param satterthwaite A logical indicating whether the Satterthwaite approximation for the degrees of freedom should be used instead of the classical trace of the Hat matrix. Defaults to \code{FALSE}.
#' @param printProgress A logical indicating whether the R should print a message before calculating the contrasts for each accession. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @param ... Additional arguments to be passed to the squeezeVarRob function internally.
#' @return A list of matrices, with each matrix in the list corresponding to a contrast in L. Each row of the matrix corresponds to a protein in the \code{\link[=protLM-class]{protLM}} object.
#' The \code{estimate} column contains the size estimate of the contrast, the \code{se} column contains the estimated standard error on the contrast, the \code{Tval} column contains the T-value corresponding to the contrast, the \code{pval} column holds the p-value corresponding to the contrast and the \code{qval} column holds the corrected p-values.
#' Each matrix is sorted from smalles to largest \code{pval} with \code{NA} values at the bottom of the matrices.
#' If \code{simplify=TRUE} and the \code{\link[=protLM-class]{protLM}} object contains only one element, the matrix is not present in a list.
#' @include testProtLMContrast.R
#' @include updateProgress.R
#' @include prot_p_adjust.R
#' @include prot_signif.R
#' @include test_contrast_adjust.R
#' @include prot_p_adjust_protwise.R
#' @export
test.contrast_stagewise <- function(protLM, L, add.annotations=TRUE, level=0.05, method_stage1="fdr", squeezeVar=TRUE, par_squeeze=NULL, min_df=1, custom_dfs=NULL, robust_var=TRUE, simplify=TRUE, lfc=0, exp_unit=NULL, pars_df=NULL, satterthwaite=FALSE, printProgress=FALSE, shiny=FALSE, message=NULL, ...)
{

  contrast_S1 <- test.contrast_adjust(protLM, L, level=level, method=method_stage1, add.annotations=add.annotations, squeezeVar=squeezeVar, par_squeeze=par_squeeze, min_df=min_df, custom_dfs=custom_dfs, robust_var=robust_var, simplify=TRUE, lfc=lfc, anova=TRUE, anova.na.ignore=TRUE, exp_unit=exp_unit, pars_df=pars_df, satterthwaite=satterthwaite, printProgress=printProgress, shiny=shiny, message=message)

  #Extract the ones that are significant in stage 1
  sign_setS1 <- subset(contrast_S1, contrast_S1[,"signif"]==1)
  significantS1 <- rownames(sign_setS1)

  contrasts <- test.protLMcontrast(protLM, L, add.annotations=add.annotations, squeezeVar=squeezeVar, par_squeeze=par_squeeze, min_df=min_df, custom_dfs=custom_dfs, robust_var=robust_var, simplify=FALSE, lfc=lfc, exp_unit=exp_unit, pars_df=pars_df, satterthwaite=satterthwaite, printProgress=printProgress, shiny=shiny, message=message)

  n_ann <- ncol(getAnnotations(protLM))

  contrasts <- lapply(contrasts, function(x){
    #Order each matrix according to the first stage
    x <- x[match(rownames(contrast_S1), rownames(x)),]
    #Add the pval, qval and signif columns of the first stage at the appropriate place in the data frame
    x <- data.frame(x[,1:(4+n_ann)],pvalS1=contrast_S1[,"pval"],qvalS1=contrast_S1[,"qval"],signifS1=contrast_S1[,"signif"],pval=x[,(5+n_ann)])
    return(x)})

  contrasts <- prot.p.adjust_protwise(contrasts, L, stage2=TRUE, significant_stage1=significantS1)
  contrasts <- prot.signif(contrasts, level=level)

  #Add the results of the ANOVA
  contrasts$ANOVA <- contrast_S1

  if(isTRUE(simplify) & length(contrasts)==1)
  {
    contrasts <- contrasts[[1]]
  }

  return(contrasts)
}
