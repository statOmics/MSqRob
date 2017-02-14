#' ANOVA for multiple contrasts
#'
#' @description This function can test an ANOVA for multiple contrasts based on parameter estimates \code{beta}, a variance covariance matrix \code{vcov}, a number of degrees of freedom \code{df}, a residual error term \code{sigma} and a contrast matrix \code{L}.
#' @param beta A matrix with a single column containing parameter estimates.
#' @param vcov A matrix of the estimated covariances between the parameter estimates.
#' @param df A numeric value indicating the residual degrees of freedom.
#' @param sigma A numeric value indicating a residual standard deviation.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param anova.na.ignore A logical indicating whether contrasts that cannot be fitted due to missing observations should be ignored when calculating an F value (\code{anova.na.ignore = TRUE}) or NA should be returned when at least one contrast cannot be estimated (\code{anova.na.ignore = FALSE}). Defaults to \code{TRUE}.
#' @return A matrix with a single row and with three columns: (1) \code{AveExpr}, containing the average estimate of the contrasts; (2) \code{Fval}, containing the F-value corresponding to the contrasts; (3) \code{pval}, containing the p-value corresponding to the contrasts.
#' @examples ....
#' @references ....
#' @include protdata.R
#' @include protLM.R
#' @include getBetaVcovDf.R
#' @include getBetaVcovDfList.R
#' @export
test.ANOVA=function(beta, vcov, df, sigma, L, anova.na.ignore=TRUE)
{

  #If anova.na.ignore==TRUE, remove NA contrast columns from L
  #If all columns in L are NA, keep NA
  if(isTRUE(anova.na.ignore) &  !all(is.na(L))){
  L <- L[,!is.na(L[1,]), drop=FALSE]
  }

  #If there are still any NA values in L => means that anova.na.ignore=FALSE or all(is.na(L)) => everything should be NA
  if(any(is.na(L))){L[L==L] <- NA}

  #Remove columns that make L less than full rank.
  #It doesn't matter which columns you remove (as long as the rank stays the same), the result will always be the same.
  if(!any(is.na(L))){
  qrc <- qr(L)
  ncontrasts <- qrc$rank
  if (ncontrasts == 0)
    stop("contrasts are all zero")
  coef <- 1:ncontrasts
  if (ncontrasts < ncol(L))
    L <- L[, qrc$pivot[coef]]
  }

  #Estimates
  estimate <- as.double(crossprod(L,beta))

  #AveExpr
  AveExpr <- mean(estimate) #, na.rm=anova.na.ignore

  #F value
  Fval <- as.double(estimate%*%solve(t(L)%*%(vcov*sigma^2)%*%L)%*%estimate)

  #df numerator:
  df_num <- tryCatch(Matrix::rankMatrix(L, method="qrLINPACK"), error=function(e){return(as.numeric(NA))})

  #df denominator: equal to t-test
  df_den <- df

  #p value
  pval <- pf(Fval, df_num, df_den, lower.tail=FALSE)

  returnmatrix <- cbind(AveExpr, df_num, df_den, Fval, pval)

  return(returnmatrix)
}

