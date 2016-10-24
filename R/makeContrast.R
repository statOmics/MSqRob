#' Create a contrast matrix
#'
#' @description This function returns a contrast matrix corresponding to specified contrasts of a set of parameters.
#' In order to view all possible predictor levels of a specific model, use attr(mymodel,"MSqRob_levels"), e.g. attr(getModels(modelRRCPTAC[1]),"MSqRob_levels").
#' The matrix specifies which comparisons between the coefficients are to be extracted from the fit.
#' The output from this function is usually used as input to the \code{\link{test.protLMcontrast}} function.
#' @param contrasts A character vector specifying containing each contrast of interest.
#' @param levels A character vector containing the parameter levels to be used in 1 or more contrasts.
#' @return A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @examples #Create a contrast matrix L for investigating whether all 10 pairwise contrasts between conditions conc6A, conc6B, conc6D and conc6E differ from zero:
#' L <- makeContrast(contrasts=c("conc6B-conc6A","conc6C-conc6A","conc6D-conc6A","conc6E-conc6A","conc6C-conc6B","conc6D-conc6B","conc6E-conc6B","conc6D-conc6C","conc6E-conc6C","conc6E-conc6D"),
#' levels=c("conc6A","conc6B","conc6C","conc6D","conc6E"))
#' #Create a contrast matrix L for investigating whether the size of the effects conc6A, conc6B, conc6D and conc6E differs from zero:
#' L <- makeContrast(contrasts=c("conc6A","conc6B","conc6C","conc6D","conc6E"),
#' levels=c("conc6A","conc6B","conc6C","conc6D","conc6E"))
#' #Create a contrast matrix L for investigating whether the size of conc6A minus 2 times conc6B differs from zero:
#' L <- makeContrast(contrasts=c("conc6A-2*conc6B"),
#' levels=c("conc6A","conc6B","conc6C","conc6D","conc6E"))
#' @export
makeContrast <- function(contrasts, levels)
{
  n <- length(levels)
  if (n < 1)
    stop("No levels to construct contrasts from.")
  if (is.factor(levels))
    levels <- levels(levels)
  if (!is.character(levels))
    levels <- colnames(levels)

  indicator <- function(i,n) {
         out <- rep(0,n)
         out[i] <- 1
         out
    }

  if (!is.null(contrasts)) {

    #In order to remove invalid level values which are not syntactically valid variable names in R and, for example, do not begin with a letter.
    from <- levels
    to <- make.names(levels)

    gsub2 <- function(pattern, replacement, x, ...) {
      for(i in 1:length(pattern))
        x <- gsub(pattern[i], replacement[i], x, ...)
      x
    }

    e <- gsub2(from, to, contrasts)

    levelsenv <- new.env()
    for (i in 1:n) assign(to[i], indicator(i, n), pos = levelsenv)

    ne <- length(contrasts)
    L <- matrix(0, nrow=length(levels), ncol=length(contrasts), dimnames=list(Levels = levels, Contrasts = contrasts))

    for (j in 1:ne) {
      ej <- parse(text = e[j])
      L[, j] <- eval(ej, envir = levelsenv)
    }

  }

  return(L)
}
