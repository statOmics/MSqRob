#' Adjust P-values for Multiple Comparisons
#'
#' Given data frame or a list of matrices with a column containing p-values, adds one or more columns containing p-values adjusted using one or more methods.
#' @param coefmatlist Either a data frame containing a column with p-values or a list of such matrices. A column with T-values is also necessary when \code{method = "fdrtool"}.
#' @param method	Correction method. Can be abbreviated. Defaults to \code{"fdr"}. To get all available methods, type \code{c("fdrtool", p.adjust.methods)}. For more information on the "fdrtool" method, check out the package \code{fdrtool} by Bernd Klaus and Korbinian Strimmer. For more information on the other methods, see the \code{\link{p.adjust}} function.
#' @param pcol A character string or numeric index indicating the column containing the p-values. Defaults to \code{"pval"}.
#' @param qname A character string indicating the name that will be given to the column containing the q-values. Defaults to \code{"qval"}.
#' @param Tcol Only needed when \code{method = "fdrtool"}. A character string or numeric index indicating the column containing the T-values. Defaults to \code{NULL}.
#' @param threshold_excess1 A numeric value between 0 and 1 indicating the threshold at which correcting for an excess in p-values close to 1 should be performed. This is sometimes needed when a null distribution needs to be fitted. Note that the method still accounts for the total number of p-values when performing this correction. Only if there is, according to a binomial test, a significant enrichment of p-values above an automatically determined cut-off at the \code{threshold_excess1} significance level, a correction for enrichment in  p-values close to 1 will be performed. If you never want to correct for an excess in p-values close to 1, set \code{threshold_excess1} to \code{0}. If you always want to correct for an excess in p-values close to 1, set \code{threshold_excess1} to \code{1}. Defaults to \code{1e-90}.
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
#' @include fdrtool_subset.R
#' @export
prot.p.adjust <- function(coefmatlist, method = "fdr", pcol="pval", qname="qval", Tcol=NULL, threshold_excess1=1e-90)
{

  method <- match.arg(method, choices=c("fdrtool", p.adjust.methods), several.ok = TRUE)

  if(is.null(Tcol) & "fdrtool" %in% method){stop("Method \"fdrtool\" requires that you provide a column with T values in the \"Tcol\" argument.")}

  nolist <- FALSE

  if(!is.list(coefmatlist) || is.data.frame(coefmatlist)){
    coefmatlist <- list(coefmatlist)
    nolist <- TRUE}

    coefmatlist <- lapply(coefmatlist, function(x) {

      pvals <- x[,pcol]

      if(is.null(Tcol) || !(Tcol %in% colnames(x))){signs <- NULL
      } else{signs <- sign(x[,Tcol])}

      ###Order p values from small to large###
      ord <- order(pvals)

      pvals <- pvals[ord]
      signs <- signs[ord]
      x <- x[ord,]
      ########################################

      #Adjust p values and add them as an extra column
      qvals <- get_qvals(pvals = pvals, signs=signs, method = method, threshold_excess1 = threshold_excess1)
      y <- cbind.data.frame(x, qvals)
      if(length(method)==1){colnames(y)[ncol(y)] <- qname}
      attr(y, "MSqRob_fdrtool") <- attr(qvals, "MSqRob_fdrtool")

      ###Put data frame back in original order###
      y <- y[order(ord),]
      ###########################################

      return(y)
    })

    if(isTRUE(nolist)){coefmatlist <- coefmatlist[[1]]}

  #else{stop("coefmatlist should either be a matrix or a list of matrices.")}

  return(coefmatlist)
}



#' Adjust P-values for Multiple Comparisons
#'
#' Given data frame or a vector of p-values, creates a data frame with column(s) containing p-values adjusted using one or more methods.
#' @param pvals A numeric vector of p-values (possibly with \link{NA}s).
#' @param signs A numeric vector containing only -1, 1, 0 and \link{NA}. It indicates the sign of the corresponding T statistics.
#' @param method	Correction method. Can be abbreviated. Defaults to \code{"fdr"}. To get all available methods, type \code{c("fdrtool", p.adjust.methods)}. For more information on the "fdrtool" method, check out the package \code{fdrtool} by Bernd Klaus and Korbinian Strimmer. For more information on the other methods, see the \code{\link{p.adjust}} function.
#' @param threshold_excess1 A numeric value between 0 and 1 indicating the threshold at which correcting for an excess in p-values close to 1 should be performed. This is sometimes needed when a null distribution needs to be fitted. Note that the method still accounts for the total number of p-values when performing this correction. Only if there is, according to a binomial test, a significant enrichment of p-values above an automatically determined cut-off at the \code{threshold_excess1} significance level, a correction for enrichment in  p-values close to 1 will be performed. If you never want to correct for an excess in p-values close to 1, set \code{threshold_excess1} to \code{0}. If you always want to correct for an excess in p-values close to 1, set \code{threshold_excess1} to \code{1}. Defaults to \code{1e-90}.
#' @return A data frame containing either one column with corrected p-values, or when \code{method = "fdrtool"}, a data frame with three columns: p-values adjusted for scale and proportion of null values, q-values (FDR) and local fdr.
#' @examples .....
#' @export
get_qvals <- function(pvals, signs, method, threshold_excess1){

  n_tot <- length(pvals) #length of pvals with NA

  no_NA_indices <- !(is.na(pvals)) #Must be like this, otherwise pval is changed in the second rule

  #If all p values are NA
  if(sum(no_NA_indices)==0){

    pvals_red <- numeric(0)
    signs_red <- numeric(0)
    excess1 <- numeric(0)
    n2 <- 0
    n_NA <- length(pvals)

  #If not, proceed with normal calculations
  } else{

  pvals <- pvals[no_NA_indices]
  signs <- signs[no_NA_indices]

  if(!is.null(signs) & suppressWarnings(any(is.na(signs)))){stop("Some signs contain NA values for valid p-values. Please check your input. If those p-values should not be taken into account, please set them to NA.")}

  n2 <- length(pvals) #length of pvals without NA
  n_NA <- n_tot-n2 #number of NA's in pvals

  pvals_red <- pvals #if there is no excess in 1's => reduced p-values==p-values
  signs_red <- signs

    #Test if there is an enrichment

    cutoffs <- cutOffPval(pvals)
    enrich_pval <- (1-pbinom((sum(pvals>cutoffs["cutoff_l"], na.rm=TRUE)), size=n2, prob=(1-cutoffs["cutoff_l"])))

    if(enrich_pval<threshold_excess1){
      #p-waarden die naar FDR gaan aanpassen
      #geef melding dat ze aangepast zijn
      ###bv.: all p values larger than... are set to 1
      pvals_red <- pvals[pvals<cutoffs["cutoff_u"]] #-qnorm(p[p<cutoffs["cutoff_u"]]/2)*signs[p<cutoffs["cutoff_u"]]
      signs_red <- signs_red[pvals<cutoffs["cutoff_u"]]
    }

    excess1 <- rep(1,n2-length(pvals_red)) #A vector with the excess 1's

  }

  #Loop to calculate all q values, loop is only needed in the rare case that more than one method is specified.
  qvals <- data.frame(matrix(ncol = 0, nrow = n_tot))
  for(i in 1:length(method)){
    qvals2 <- calculate_qvals(pvals, pvals_red, signs, signs_red, method[i], excess1, n2, n_NA)
    qvals <- cbind.data.frame(qvals,qvals2)
    if(is.null(attr(qvals, "MSqRob_fdrtool"))){attr(qvals, "MSqRob_fdrtool") <- attr(qvals2, "MSqRob_fdrtool")}
  }

  return(qvals) #data frame

}
#####################

calculate_qvals <- function(pvals, pvals_red, signs, signs_red, method, excess1, n2, n_NA){

  if(length(method)!=1){stop("Use only one method at a time with this function!")}

  if(method=="fdrtool"){

    #!!!transformeren naar z values!!! -> anders geen shift!!!
    x <- as.vector(-qnorm(pvals/2)*signs)
    x_red <- as.vector(-qnorm(pvals_red/2)*signs_red)

    MSqRob_fdrtool <- fdrtool_subset(x,x_red)
    #geef melding dat ze aangepast zijn
    ###bv.: all p values larger than... are set to 1

    #excess1 must not be added here because it's already included in fdrtool_subset:
    qvals <- cbind.data.frame(pval_adj = c(MSqRob_fdrtool$pval, rep(NA, n_NA)), lfdr = c(MSqRob_fdrtool$lfdr, rep(NA, n_NA)), qval = c(MSqRob_fdrtool$qval, rep(NA, n_NA)))

    #Pass on this result as an attribute in order to be able to plot later on!
    attr(qvals, "MSqRob_fdrtool") <- MSqRob_fdrtool

  } else{

    #To prevent an annoying error when all p values are NA: p.adjust works for numeric(0), but not if n is set to 0...
    if(n2==0){p_adjusted <- numeric(0)} else{
    p_adjusted <- p.adjust(pvals_red, method = method, n = n2)}

    qvals <- data.frame(c(p_adjusted, excess1, rep(NA, n_NA)))
    colnames(qvals) <- method
    ###all(c(p.adjust(p_red, n=length(p)),rep(1,length(p)-length(p_red)))==p.adjust(p)) geeft TRUE!
    #-> only in the exceptional case when the removed p values have no FDR=1, which is our intention!
    #because they must be always set to 1 if enrich_pval<threshold_excess1!
  }

  return(qvals)

}

#' Find cut off for p-values enriched in values close to 1
#'
#' Contrast tests based on mixed models often produce p-values containing an excess of values close to 1.
#' This function will, given a vector of p-values find a cut-off above which the values are excessively enriched.
#' The cut-off is determined based on the difference coefficient and 2nd order difference coefficient of sequential p-values.
#' Therefore, it will (almost) always find a cut-off, even if there is no real enrichment in p-values.
#' One should thus always assess whether there is truly an enrichment when using this function.
#' @param p A vector of p-values.
#' @return A named vector with two values: "cutoff_l", the maximal p-value below the cut-off and "cutoff_u", the minimal p-value that is above the cut-off.
#' @examples .......
#' @export
cutOffPval <- function(p){

  #Old: order p values from small to big with NA at the end
  #New: we only want to find the cut off => remove NA values
  p <- na.omit(p)
  #We don't return p, so we don't need to put it back into the right order!!!
  #names(p) <- 1:length(p)
  p <- p[order(p)]

  ###Set p values that are too close to 1 to 1, find cut-off

  #Log-transform p values
  data <- log(1-p)

  derivMatrix <- cbind((1:(length(p)-1)),diff(data), c(diff(diff(data)),NA))
  derivMatrix <- na.omit(derivMatrix[order(derivMatrix[,2]),])

  found <- FALSE
  j <- 1
  while(!isTRUE(found) && j<=nrow(derivMatrix)){
    i <- derivMatrix[j,1]
    #i cannot be 1, otherwise the second part is logical(0): the || prevents getting a logical(0)!!!

    index1 <- which(derivMatrix[,1]==i-1)
    index2 <- which(derivMatrix[,1]==i)

    if(i==1 || !(derivMatrix[index1,3]<0 & derivMatrix[index2,3]>0)){j <- j+1
    }else{found <- TRUE
    }
  }

  if(isTRUE(found)){
    #remove names! Needed for names of "cutoffs", otherwise they are appended!
    p <- unname(p)
    result <- c(cutoff_l=p[i], cutoff_u=p[i+1])
  } else{result <- c(cutoff_l=1, cutoff_u=1)}

  return(result)

}
