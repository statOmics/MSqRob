#' Adjust P-values for Multiple Comparisons
#'
#' Given data frame or a list of matrices with a column containing p-values, adds a column containing p-values adjusted using one of several methods.
#' @param coefmatlist Either a data frame containing a column with p-values or a list of such matrices.
#' @param method	Correction method. Can be abbreviated. Defaults to "fdr". To get all available methods, type \code{p.adjust.methods}. For more information on these methods, see the \code{\link{p.adjust}} function.
#' @param pcol A character string or numeric index indicating the column containing the p-values. Defaults to "pval".
#' @param qname A character string indicating the name that will be given to the column containing the q-values. Defaults to "qval".
#' @return Depending on the input, either a data frame or a list of data frames with an extra numeric column containing corrected p-values.
#' @examples #Create a contrast matrix L for investigating whether all 10 pairwise contrasts between conditions conc6A, conc6B, conc6D and conc6E differ from zero:
#' L <- makeContrast(contrasts=c("conc6B-conc6A","conc6C-conc6A","conc6D-conc6A","conc6E-conc6A","conc6C-conc6B","conc6D-conc6B","conc6E-conc6B","conc6D-conc6C","conc6E-conc6C","conc6E-conc6D"),
#' levels=c("conc6A","conc6B","conc6C","conc6D","conc6E"))
#' #Load the protLM object protmodel:
#' protmodel <- data(modelRRCPTAC, package="MSqRob")
#' #Test the contrast L:
#' protvalues <- test.protLMcontrast(protmodel,L)
#' #Adjust the p-values by Benjamini-Hochberg FDR:
#' protvalues <- prot.p.adjust(protvalues)
#' @export
prot.p.adjust <- function(coefmatlist, method = "fdr", pcol="pval", qname="qval")
{

  nolist <- FALSE

  if(!is.list(coefmatlist) || is.data.frame(coefmatlist)){
    coefmatlist <- list(coefmatlist)
    nolist <- TRUE}

    coefmatlist <- lapply(coefmatlist, function(x) {
      y <- cbind.data.frame(x, p.adjust(x[,pcol], method = method))
      colnames(y)[ncol(y)] <- qname
      return(y)
    })

    if(isTRUE(nolist)){coefmatlist <- coefmatlist[[1]]}

  #else{stop("coefmatlist should either be a matrix or a list of matrices.")}

  return(coefmatlist)
}
