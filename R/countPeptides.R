#' Count the number of peptides in each group
#'
#' @description Given a \code{\link[=protdata-class]{protdata}} object, a dataframe pData containing experiment information and a predictor of interest,
#' returns a data frame giving the amount of peptides observed for each level of the predictor.
#' @param protdata A \code{\link[=protdata-class]{protdata}} object
#' @param pData A dataframe that contains experimenter-supplied variables describing phenotypes and/or experimental factors that might influence the quantitative value of interest for each mass spec run. The \code{colnames} should be all the different effects that influence the outcome. The first column should hold all the different run names. The corresponding levels of the other effects should be in the subsequent columns.
#' @param predictor A predictor for which you would like to summarize the number of identified peptides per level for each protein. This predictor should correspond to one element of the \code{colnames} of \code{pData}.
#' @include protdata.R
#' @export
countPeptides <- function(protdata, pData, predictor){

data <- getData(protdata)

pred_lvls <- levels(pData[, predictor])

countframe <- as.data.frame(matrix(ncol=length(pred_lvls), nrow=0, dimnames=list(NULL, pred_lvls)))

for(i in 1:length(protdata))
{
  for(j in 1:length(pred_lvls)){
    countframe[i,pred_lvls[j]] <- sum(data[[i]][,predictor] %in% pred_lvls[j])
  }
}

rownames(countframe) <- getAccessions(protdata)

return(countframe)
}

#pData <- pData(peptides)
#protdata <- proteins
#predictor <- "treat"


