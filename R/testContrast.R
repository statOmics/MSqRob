#' Test a contrast
#'
#' @description This function can test a contrast based on parameter estimates \code{beta}, a variance covariance matrix \code{vcov}, a number of degrees of freedom \code{df}, a residual error term \code{sigma} and a contrast matrix \code{L}.
#' @param beta A matrix with a single column containing parameter estimates.
#' @param vcov A matrix of the estimated covariances between the parameter estimates.
#' @param df A vector of numeric values indicating the degrees of freedom corresponding to each contrast or a single numeric value if the degrees of freedom for each contrast are equal.
#' @param sigma A numeric value indicating a residual standard deviation.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param lfc The minimum log2-fold-change that is considered scientifically meaningful. Defaults to \code{0}.
#' @return A matrix with a row for each contrast in \code{L} and with four columns: (1) \code{estimate}, containing the size estimates of the contrasts; (2) \code{se}, containing the estimated standard error on the contrast; (3) \code{Tval}, containing the T-value corresponding to the contrast; (4) \code{pval}, containing the p-value corresponding to the contrast.
#' @examples #Load the protLM object protmodel:
#' protmodel <- data(modelRR, package="MSqRob")
#' betaVcovDfList <- getbetaVcovDfList(protmodel)
#' contrasts <- test.contrast(betaVcovDfList)
#' contrasts
#' @references David Ruppert, M.P. Want and R.J. Carroll.
#' Semiparametric Regression.
#'  Cambridge University Press, 2003.
#' @include protdata.R
#' @include protLM.R
#' @include getBetaVcovDf.R
#' @include getBetaVcovDfList.R
#' @export
test.contrast=function(beta, vcov, df, sigma, L, lfc=0)
{
  estimate <- as.double(crossprod(L,beta)) #=as.double(t(L)%*%beta)=colSums(L*beta)
  se <- Matrix::diag(t(L)%*%(vcov*sigma^2)%*%L)^.5

  lfc <- abs(lfc)
  aest <- abs(estimate)

  Tval <- setNames(rep(0, length(estimate)),names(se))
  tstat.right <- (aest - lfc)/se
  tstat.left <- (aest + lfc)/se
  pval <- pt(tstat.right, df = df, lower.tail = FALSE) +
    pt(tstat.left, df = df, lower.tail = FALSE)
  tstat.right <- pmax(tstat.right, 0)
  fc.up <- (estimate >= lfc)
  fc.up[is.na(fc.up)] <- FALSE
  fc.down <- (estimate < -lfc)
  fc.down[is.na(fc.down)] <- FALSE
  Tval[fc.up] <- tstat.right[fc.up]
  Tval[fc.down] <- -tstat.right[fc.down]
  Tval[is.na(estimate)] <- NA

  #Tval <- estimate/se
  #pval <- 2*pt(-abs(Tval),df)

  returnmatrix <- cbind(estimate, se, df, Tval, pval)
  rownames(returnmatrix) <- colnames(L)
  return(returnmatrix)
}

