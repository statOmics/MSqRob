fitFDist_LG <- function(x, df1, covariate=NULL, min_df=1)
  #  Moment estimation of the parameters of a scaled F-distribution
  #	The first degrees of freedom are given
  #	Gordon Smyth and Belinda Phipson
  #	8 Sept 2002.  Last revised 27 Oct 2012.
{
  #	Check covariate
  if(!is.null(covariate)) {
    if(length(covariate) != length(x)) stop("covariate and x must be of same length")
    if(any(is.na(covariate))) stop("NA covariate values not allowed")
    isfin <- is.finite(covariate)
    if(!all(isfin)) {
      if(!any(isfin))
        covariate <- sign(covariate)
      else {
        r <- range(covariate[isfin])
        covariate[covariate == -Inf] <- r[1]-1
        covariate[covariate == Inf] <- r[2]+1
      }
    }
    splinedf <- min(4,length(unique(covariate)))
    if(splinedf < 2) covariate <- NULL
  }
  #	Remove missing or infinite values and zero degrees of freedom
  ok <- is.finite(x) & is.finite(df1) & (x > -1e-15) & (df1 > min_df) #1e-15
  notallok <- !all(ok)
  if(notallok) {
    x <- x[ok]
    df1 <- df1[ok]
    if(!is.null(covariate)) {
      covariate2 <- covariate[!ok]
      covariate <- covariate[ok]
    }
  }
  n <- length(x)
  if(n==0) return(list(scale=NA,df2=NA))

  #	Avoid exactly zero values
  x <- pmax(x,0)
  m <- median(x)
  if(m==0) {
    warning("More than half of residual variances are exactly zero: eBayes unreliable")
    m <- 1
  } else {
    if(any(x==0)) warning("Zero sample variances detected, have been offset",call.=FALSE)
  }
  x <- pmax(x, 1e-5 * m)

  #	Better to work on with log(F)
  z <- log(x)
  e <- z-digamma(df1/2)+log(df1/2)

  if(is.null(covariate)) {
    emean <- mean(e)
    evar <- sum((e-emean)^2)/(n-1) # equals evar <- var(e) if emean=mean(e)

  } else {
    require(splines)
    design <- try(ns(covariate,df=splinedf,intercept=TRUE),silent=TRUE)
    if(is(design,"try-error")) stop("Problem with covariate")
    fit <- lm.fit(design,e)
    if(notallok) {
      design2 <- predict(design,newx=covariate2)
      emean <- rep.int(0,n+length(covariate2))
      emean[ok] <- fit$fitted
      emean[!ok] <- design2 %*% fit$coefficients
    } else {
      emean <- fit$fitted
    }
    evar <- mean(fit$residuals[-(1:fit$rank)]^2)
  }

  #MSqRob added: avoid NaN in evar
  if(n==1) evar <- 0

  evar <- evar - mean(trigamma(df1/2))
  if(evar > 0) {
    df2 <- 2*limma::trigammaInverse(evar)
    s20 <- exp(emean+digamma(df2/2)-log(df2/2))
  } else {
    df2 <- Inf
    s20 <- exp(emean)
  }
  list(scale=s20,df2=df2)
}


#Just adapted to use fitFDist_LG internally instead of fitFDist
fitFDistRobustly_LG <- function(x, df1, covariate=NULL, winsor.tail.p=c(0.05,0.1), trace=FALSE, min_df=1)
  #	Robust estimation of the parameters of a scaled F-distribution,
  #	given the first degrees of freedom, using first and second
  #	moments of Winsorized z-values
  #	Gordon Smyth and Belinda Phipson
  #	Created 7 Oct 2012.  Last revised 20 November 2015.
{
  #	Check x
  n <- length(x)

  #	Eliminate cases of no useful data
  if(n<2) return(list(scale=NA,df2=NA,df2.shrunk=NA))
  if(n==2) return(fitFDist_LG(x=x,df1=df1,covariate=covariate,min_df=min_df))

  #	Check df1
  if(all(length(df1) != c(1,n))) stop("x and df1 are different lengths")

  #	Check covariate
  if(!is.null(covariate)) {
    if(length(covariate) != n) stop("x and covariate are different lengths")
    if(!all(is.finite(covariate))) stop("covariate contains NA or infinite values")
  }


  #	Treat zero df1 values as non-informative cases
  #	Similarly for missing values or x or df1
  ok <- !is.na(x) & is.finite(df1) & (df1 > min_df) #1e-6
  notallok <- !all(ok)
  if(notallok) {
    df2.shrunk <- x
    x <- x[ok]
    if(length(df1)>1) df1 <- df1[ok]
    if(!is.null(covariate)) {
      covariate2 <- covariate[!ok]
      covariate <- covariate[ok]
    }
    fit <- Recall(x=x,df1=df1,covariate=covariate,winsor.tail.p=winsor.tail.p,trace=trace)
    df2.shrunk[ok] <- fit$df2.shrunk
    df2.shrunk[!ok] <- fit$df2
    if(is.null(covariate))
      scale <- fit$scale
    else {
      scale <- x
      scale[ok] <- fit$scale
      scale[!ok] <- exp(approx(covariate,log(fit$scale),xout=covariate2,rule=2)$y)
    }
    return(list(scale=scale,df2=fit$df2,df2.shrunk=df2.shrunk))
  }

  #	Avoid zero or negative x values
  m <- median(x)
  if(m<=0)	stop("x values are mostly <= 0")
  i <- (x < m*1e-12)
  if(any(i)) {
    warning("small x values have been offset away from zero")
    x[i] <- m*1e-12
  }

  #	Store non-robust estimates
  NonRobust <- fitFDist_LG(x=x,df1=df1,covariate=covariate)

  #	Check winsor.tail.p
  prob <- winsor.tail.p <- rep(winsor.tail.p,length=2)
  prob[2] <- 1-winsor.tail.p[2]
  if(all(winsor.tail.p < 1/n)) {
    NonRobust$df2.shrunk <- rep.int(NonRobust$df2,n)
    return(NonRobust)
  }

  #	Transform x to constant df1
  if(length(df1)>1) {
    df1max <- max(df1)
    i <- (df1 < (df1max - 1e-14))
    if(any(i)) {
      if(is.null(covariate)) s <- NonRobust$scale else s <- NonRobust$scale[i]
      f <- x[i]/s
      df2 <- NonRobust$df2
      pupper <- pf(f,df1=df1[i],df2=df2,lower.tail=FALSE,log.p=TRUE)
      plower <- pf(f,df1=df1[i],df2=df2,lower.tail=TRUE,log.p=TRUE)
      up <- pupper<plower
      if(any(up)) f[up] <- qf(pupper[up],df1=df1max,df2=df2,lower.tail=FALSE,log.p=TRUE)
      if(any(!up)) f[!up] <- qf(plower[!up],df1=df1max,df2=df2,lower.tail=TRUE,log.p=TRUE)
      x[i] <- f*s
      df1 <- df1max
    } else {
      df1 <- df1[1]
    }
  }

  #	Better to work with log(F)
  z <- log(x)

  #	Demean or Detrend
  if(is.null(covariate)) {
    ztrend <- mean(z,trim=winsor.tail.p[2])
    zresid <- z-ztrend
  } else {
    lo <- loessFit(z,covariate,span=0.4)
    ztrend <- lo$fitted
    zresid <- lo$residual
  }

  #	Moments of Winsorized residuals
  zrq <- quantile(zresid,prob=prob)
  zwins <- pmin(pmax(zresid,zrq[1]),zrq[2])

  # !!! Changed to make it more robust (mainly variance was a problem!)

  #zwmean <- mean(zwins) #median(zwins)==median(zresid)
  zwmean <- median(zwins)
  #zwvar <- mean((zwins-zwmean)^2)*n/(n-1) #mad(zwins)==mad(zresid)
  zwvar <- mad(zwins)

  if(trace) cat("Variance of Winsorized Fisher-z",zwvar,"\n")

  #	Theoretical Winsorized moments
  if(!requireNamespace("statmod",quietly=TRUE)) stop("statmod package required but is not installed")
  g <- statmod::gauss.quad.prob(128,dist="uniform")
  linkfun <- function(x) x/(1+x)
  linkinv <- function(x) x/(1-x)
  winsorizedMoments <- function(df1=df1,df2=df2,winsor.tail.p=winsor.tail.p) {
    fq <- qf(p=c(winsor.tail.p[1],1-winsor.tail.p[2]),df1=df1,df2=df2)
    zq <- log(fq)
    q <- linkfun(fq)
    nodes <- q[1] + (q[2]-q[1]) * g$nodes
    fnodes <- linkinv(nodes)
    znodes <- log(fnodes)
    f <- df(fnodes,df1=df1,df2=df2)/(1-nodes)^2
    q21 <- q[2]-q[1]
    m <- q21*sum(g$weights*f*znodes) + sum(zq * winsor.tail.p)
    v <- q21*sum(g$weights*f*(znodes-m)^2) + sum((zq-m)^2 * winsor.tail.p)
    list(mean=m,var=v)
  }

  #	Try df2==Inf
  mom <- winsorizedMoments(df1=df1,df2=Inf,winsor.tail.p=winsor.tail.p)
  funvalInf <- log(zwvar/mom$var)
  if(funvalInf <= 0) {
    df2 <- Inf
    #		Correct trend for bias
    ztrendcorrected <- ztrend+zwmean-mom$mean
    s20 <- exp(ztrendcorrected)
    #		Posterior df for outliers
    Fstat <- exp(z-ztrendcorrected)
    TailP <- pchisq(Fstat*df1,df=df1,lower.tail=FALSE)
    r <- rank(Fstat)
    EmpiricalTailProb <- (n-r+0.5)/n
    ProbNotOutlier <- pmin(TailP/EmpiricalTailProb,1)
    df.pooled <- n*df1
    df2.shrunk <- rep.int(df2,n)
    O <- ProbNotOutlier < 1
    if(any(O)) {
      df2.shrunk[O] <- ProbNotOutlier[O]*df.pooled
      o <- order(TailP)
      df2.shrunk[o] <- cummax(df2.shrunk[o])
    }
    return(list(scale=s20,df2=df2,df2.shrunk=df2.shrunk))
  }

  #	Estimate df2 by matching variance of zwins
  #	Use beta distribution Gaussian quadrature to find mean and variance
  #	of Winsorized F-distribution
  fun <- function(x) {
    df2 <- linkinv(x)
    mom <- winsorizedMoments(df1=df1,df2=df2,winsor.tail.p=winsor.tail.p)
    if(trace) cat("df2=",df2,", Working Var=",mom$var,"\n")
    log(zwvar/mom$var)
  }

  #	Use non-robust estimate as lower bound for df2
  if(NonRobust$df2==Inf) {
    NonRobust$df2.shrunk <- rep.int(NonRobust$df2,n)
    return(NonRobust)
  }
  rbx <- linkfun(NonRobust$df2)
  funvalLow <- fun(rbx)
  if(funvalLow >= 0) {
    df2 <- NonRobust$df2
  } else {
    u <- uniroot(fun,interval=c(rbx,1),tol=1e-8,f.lower=funvalLow,f.upper=funvalInf)
    df2 <- linkinv(u$root)
  }

  #	Correct ztrend for bias
  mom <- winsorizedMoments(df1=df1,df2=df2,winsor.tail.p=winsor.tail.p)
  ztrendcorrected <- ztrend+zwmean-mom$mean
  s20 <- exp(ztrendcorrected)

  #	Posterior df for outliers
  zresid <- z-ztrendcorrected
  Fstat <- exp(zresid)
  TailP <- pf(Fstat,df1=df1,df2=df2,lower.tail=FALSE)
  r <- rank(Fstat)
  EmpiricalTailProb <- (n-r+0.5)/n
  ProbNotOutlier <- pmin(TailP/EmpiricalTailProb,1)
  ProbOutlier <- 1-ProbNotOutlier
  if(any(ProbNotOutlier<1)) {
    o <- order(TailP)
    #		Old calculation for df2.outlier
    VarOutlier <- max(zresid)^2
    VarOutlier <- VarOutlier-trigamma(df1/2)
    if(trace) cat("VarOutlier",VarOutlier,"\n")
    if(VarOutlier > 0) {
      df2.outlier.old <- 2*limma::trigammaInverse(VarOutlier)
      if(trace) cat("df2.outlier.old",df2.outlier.old,"\n")
      if(df2.outlier.old < df2) {
        df2.shrunk.old <- ProbNotOutlier*df2+ProbOutlier*df2.outlier.old
        #				Make df2.shrunk.old monotonic in TailP
        df2.ordered <- df2.shrunk.old[o]
        df2.ordered[1] <- min(df2.ordered[1],NonRobust$df2)
        m <- cumsum(df2.ordered)
        m <- m/(1:n)
        imin <- which.min(m)
        df2.ordered[1:imin] <- m[imin]
        df2.shrunk.old[o] <- cummax(df2.ordered)
      }
    }

    #		New calculation for df2.outlier
    #		Find df2.outlier to make maxFstat the median of the distribution
    #		Exploit fact that log(TailP) is nearly linearly with positive 2nd deriv as a function of df2
    #		Note that minTailP and NewTailP are always less than 0.5
    df2.outlier <- log(0.5)/log(min(TailP))*df2
    #		Iterate for accuracy
    NewTailP <- pf(max(Fstat),df1=df1,df2=df2.outlier,lower.tail=FALSE)
    df2.outlier <- log(0.5)/log(NewTailP)*df2.outlier
    df2.shrunk <- ProbNotOutlier*df2+ProbOutlier*df2.outlier

    #		Force df2.shrunk to be monotonic in TailP
    o <- order(TailP)
    df2.ordered <- df2.shrunk[o]
    m <- cumsum(df2.ordered)
    m <- m/(1:n)
    imin <- which.min(m)
    df2.ordered[1:imin] <- m[imin]
    df2.shrunk[o] <- cummax(df2.ordered)

    #		Use isoreg() instead. This gives similar results.
    df2.shrunk.iso <- rep.int(df2,n)
    o <- o[1:(n/2)]
    df2.shrunk.iso[o] <- ProbNotOutlier[o]*df2+ProbOutlier[o]*df2.outlier
    df2.shrunk.iso[o] <- isoreg(TailP[o],df2.shrunk.iso[o])$yf

  } else {
    df2.outlier <- df2.outlier2 <- df2
    df2.shrunk2 <- df2.shrunk <- rep.int(df2,n)
  }

  list(scale=s20,df2=df2,tail.p.value=TailP,prob.outlier=ProbOutlier,df2.outlier=df2.outlier,df2.outlier.old=df2.outlier.old,df2.shrunk=df2.shrunk,df2.shrunk.old=df2.shrunk.old,df2.shrunk.iso=df2.shrunk.iso)

}

#' Robustly Squeeze Sample Variances
#'
#' @description This function squeezes a set of sample variances together by computing empirical Bayes posterior means in a way that is robust against the presence of very small non-integer degrees of freedom values.
#' @param var A numeric vector of independent sample variances.
#' @param df A numeric vector of degrees of freedom for the sample variances.
#' @param covariate If \code{non-NULL}, \code{var.prior} will depend on this numeric covariate. Otherwise, \code{var.prior} is constant.
#' @param robust A logical indicating wheter the estimation of \code{df.prior} and \code{var.prior} should be robustified against outlier sample variances. Defaults to \code{FALSE}.
#' @param winsor.tail.p A numeric vector of length 1 or 2, giving left and right tail proportions of \code{x} to Winsorize. Used only when \code{robust=TRUE}.
#' @param k A numeric value indicating that the calculation of the robust squeezed variances should Winsorize at \code{k} standard deviations.
#' @param min_df A numeric value indicating the minimal degrees of freedom that will be taken into account for calculating the prior degrees of freedom and prior variance.
#' @return A list with components:
#' \code{var.post} A numeric vector of posterior variances.
#' \code{var.prior} The location of prior distribution. A vector if \code{covariate} is non-\code{NULL}, otherwise a scalar.
#' \code{df.prior} The degrees of freedom of prior distribution. A vector if \code{robust=TRUE}, otherwise a scalar.
#' @export
squeezeVarRob <- function(var, df, covariate=NULL, robust=FALSE, winsor.tail.p=c(0.05,0.1), min_df=1)
  #  Empirical Bayes posterior variances
  #	Adapted from Gordon Smyth
{
  n <- length(var)
  if(n == 0) stop("var is empty")
  #If there is only one variance that is not NA:
  if(sum(!is.na(var)) == 1){
    df[!is.na(df)] <- 0
    return(list(var.post=var,var.prior=var,df.prior=df))
  }
  if(length(df)==1) {
    df <- rep.int(df,n)
  } else {
    if(length(df) != n) stop("lengths differ")
  }

  #	Estimate prior var and df
  if(robust)
    fit <- fitFDistRobustly_LG(var, df1=df, covariate=covariate, winsor.tail.p=winsor.tail.p, min_df=min_df)
  else
    fit <- fitFDist_LG(var, df1=df, covariate=covariate, min_df=min_df)

  #	Prior var will be vector if robust=TRUE, otherwise scalar
  var.prior <- fit$scale

  #	Prior df will be vector if covariate is non-NULL, otherwise scalar
  df.prior <- fit$df2.shrunk
  if(is.null(df.prior)) df.prior <- fit$df2

  #	Check estimated prior df
  if(is.null(df.prior) || any(is.na(df.prior))) stop("Could not estimate prior df")

  #	Squeeze the posterior variances
  df.total <- df + df.prior
  var[df==0] <- 0 # guard against missing or infinite values
  Infdf <- df.prior==Inf
  if(any(Infdf)) {
    var.post <- rep(var.prior,length.out=n)
    i <- which(!Infdf)
    if(length(i)) {
      if(is.null(covariate))
        s02 <- var.prior
      else
        s02 <- var.prior[i]
      var.post[i] <- (df[i]*var[i] + df.prior[i]*s02) / df.total[i]
    }
  } else {
    var.post <- (df*var + df.prior*var.prior) / df.total
  }

  list(df.prior=df.prior,var.prior=var.prior,var.post=var.post)
}
