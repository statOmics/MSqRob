#' Adjust P-values for Multiple Comparisons per Protein
#'
#' Given matrix or a list of matrices with a column containing p-values and a contrast matrix, adds a column containing Holm-adjusted p-values correcting for the number of independent contrasts.
#' This function can also be used to calculate the corrected p-values in the second stage when performing stage-wise testing.
#' @param coefmatlist Either a matrix containing a column with p-values or a list of such matrices.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param method	Correction method. Can be abbreviated. Defaults to "Holm", which is currently the only available option.
#' @param pcol A character string or numeric index indicating the column containing the p-values. Defaults to "pval".
#' @param qname A character string indicating the name that will be given to the column containing the q-values. Defaults to "qval".
#' @param stage2 Set this to \code{TRUE} when using this function in the 2nd stage of a stage-wise test.
#' @param significant_stage1 Only used when \code{stage2 = TRUE}. This should be a character vector of protein identifiers that were retained in the first stage during stage-wise testing.
#' @return Depending on the input, either a matrix or a list of matrices with an extra numeric column containing corrected p-values.
#' @examples #....
#' @export
prot.p.adjust_protwise <- function(coefmatlist, L, method = "Holm", pcol="pval", qname="qval", stage2=FALSE, significant_stage1=NULL)
{

if(is.list(coefmatlist)){

  p_frame <- data.frame(accession=character(0), stringsAsFactors = TRUE)
  for(i in 1:length(coefmatlist))
    {
    p_column <- data.frame(accession=rownames(coefmatlist[[i]]), pval=coefmatlist[[i]][,pcol], stringsAsFactors = TRUE)
    p_frame <- suppressWarnings(merge(p_frame,p_column, all=TRUE, by="accession"))
    }

  q_frame <- p_frame
  q_frame[,-1] <- NA #Initialisation: set all values for q_frame to NA
  #q_matrix <- matrix(NA, nrow=nrow(p_frame), ncol=ncol(p_frame)-1) #1 column less because of the accession number in p_frame

  for(j in 1:nrow(p_frame))
  {
  p <- as.numeric(p_frame[j,-1])

  #Remove NA values from L
  Lnona <- L
  Lnona[,is.na(p)] <- 0

  o <- order(p)
  ro <- order(o)
  #method="qrLINPACK" is needed to give rank 0 to a matrix containing only zeros
  ranks <- unlist(lapply((as.list(1:ncol(Lnona))), function(x){return(Matrix::rankMatrix(Lnona[,o[x:ncol(Lnona)], drop=FALSE], method="qrLINPACK"))}))
  q_frame[j,-1] <- pmin(1, cummax(ranks * p[o]))[ro]
  }

  for(i in 1:length(coefmatlist))
  {

    q_column <- q_frame[,i+1]
    names(q_column) <- q_frame$accession
    #If there would be different lengths:
    q_column <- q_column[q_frame$accession %in% rownames(coefmatlist[[i]])]

    if(isTRUE(stage2)){
    q_column <- adjust_stage2(q_column, significant_stage1)
    }

    q_column <- q_column[order(match(names(q_column),rownames(coefmatlist[[i]])))]

    coefmatlist[[i]] <- cbind(coefmatlist[[i]], q_column)
    colnames(coefmatlist[[i]])[ncol(coefmatlist[[i]])] <- qname
  }

  #If there is only one contrast, the q-value is equal to the p-value
  }else{

    q_column <- coefmatlist[,pcol]

    if(isTRUE(stage2)){
      q_column <- adjust_stage2(q_column, significant_stage1)
    }

    q_column <- q_column[order(match(names(q_column),rownames(coefmatlist[[i]])))]

    coefmatlist <- cbind(coefmatlist, q_column)
    colnames(coefmatlist)[ncol(coefmatlist)] <- qname
  }
  return(coefmatlist)
}

adjust_stage2 <- function(q_column, significant_stage1){

  q_retained <- q_column[names(q_column) %in% significant_stage1] #q_column[significant_stage1]
  #Remove everything that did not pass stage 1
  q_column[!(names(q_column) %in% significant_stage1)] <- NA

  #Only execute when the number of proteins in stage 1 is given.
  if(!is.null(significant_stage1)){
    q_column <- q_column*sum(!is.na(q_column))/sum(!is.na(q_retained)) #Adjust p values to FDR significance level of stage 1
  }
  return(q_column)
}



