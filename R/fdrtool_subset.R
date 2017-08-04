#' Estimate (Local) False Discovery Rates For Diverse Test Statistics
#'
#' This function is largely based on the \code{\link[fdrtool]{fdrtool}} function from the \pkg{fdrtool} package by Korbinian Strimmer (\url{http://strimmerlab.org}).
#' \pkg{fdrtool} takes a vector of z-scores (or of correlations, p-values, or t-statistics), and estimates for each case both the tail area-based Fdr as well as the density-based fdr (=q-value resp. local false discovery rate). The parameters of the null distribution are estimated adaptively from the data (except for the case of p-values where this is not necessary).
#' Our only modification is adding the possibility to provide a reduced input (e.g. an input of which the excess in p values equal to 1 are removed) for estimating the null distribution while still correcting for the total number of input values.
#' @param x	A vector of the observed test statistics.
#' @param x_red A subset vector of x that will be given to the fdrtool algorithm.
#' @param statistic	One of "normal" (default), "correlation", "pvalue". This species the null model.
#' @param verbose	A logical indicating if status messages should be printed out. Defaults to \code{TRUE}.
#' @param cutoff.method	One of "fndr" (default), "pct0", "locfdr".
#' @param pct0	The fraction of data x_red that will be used for fitting null model - only used if cutoff.method="pct0".
#' @details See package \pkg{fdrtool}.
#' @return A list with the following components:
#'   \itemize{
#'   \item \code{x} The input vector of the observed test statistics.
#'   \item \code{x_red} The subset of the input vector of the observed test statistics that was given to the fdrtool algorithm.
#'   \item \code{pval}	A vector with adjusted p-values for each element of \code{x_red}.
#'   \item \code{pval_full} A vector with adjusted p-values for each element of \code{x}.
#'   \item \code{qval}	A vector with q-values (Fdr) for each element of \code{x_red}.
#'   \item \code{qval_full} A vector with q-values (Fdr) for each element of \code{x}.
#'   \item \code{lfdr}	A vector with local fdr values for each element of \code{x_red}.
#'   \item \code{lfdr_full} A vector with local fdr values for each element of \code{x}.
#'   \item \code{statistic} The specified type of null model.
#'   \item \code{param} A vector containing the estimated parameters (\code{eta0} (the null proportion of \code{x_red}) and the free parameter of the null model).
#'   \item \code{x0} The cut off value of \code{x_red}, indicating the separation between the null model and the mixture model.
#'   \item \code{f.pval}, \code{F.pval}, \code{fdr.pval}, \code{Fdr.pval}, \code{f0}, \code{F0}, \code{get.pval}, \code{fdr}, \code{Fdr} Functions needed to make diagnostic plots.
#' }
#' @examples .......
#' @references Strimmer, K. (2008a). A unified approach to false discovery rate estimation. BMC Bioinformatics 9: 303. Available from \url{http://www.biomedcentral.com/1471-2105/9/303/}.
#'
#' Strimmer, K. (2008b). fdrtool: a versatile R package for estimating local and tail area- based false discovery rates. Bioinformatics 24: 1461-1462. Available from \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/12/1461}.
#' @export
fdrtool_subset <- function (x, x_red=x, statistic = c("normal", "correlation", "pvalue"),
          verbose = TRUE, cutoff.method = c("fndr", "pct0", "locfdr"), pct0 = 0.75)
  {
  statistic = match.arg(statistic)
  cutoff.method = match.arg(cutoff.method)
  if (is.vector(x) == FALSE | is.vector(x_red) == FALSE)
    stop("input test statistics must be given as a vector!")
  #Added: x_red must be a subset of x
  tablx <- table(x)
  tablxred <- table(x_red)
  tablx_in_xred <- tablx[names(tablxred)]
  #all(x_red %in% x) cannot be used because of possible duplicates, e.g. c(x_red,x_red[1]) %in% x will still give TRUE
  if (any(tablx_in_xred<tablxred) || !all(x_red %in% x))
    stop("x_red should either be equal to x or a subset of x.")
  ###
  if (length(x_red) < 200)
    warning("There may be too few input test statistics for reliable FDR calculations!")

  if (length(x_red)==0){
    cf.out <- rep(NA,6)
    names(cf.out) <- c("cutoff","N.cens","eta0","eta0.SE","sd","sd.SE")
    x0 <- cf.out["cutoff"]
    result = list(x=numeric(0), x_red=numeric(0), pval = numeric(0), qval = numeric(0), lfdr = numeric(0), statistic = "normal",
                  param = cf.out, x0 = x0, f.pval = NULL, F.pval = NULL, fdr.pval=NULL, Fdr.pval=NULL, f0=NULL, F0=NULL, get.pval=NULL, fdr=NULL, Fdr=NULL)
  } else{

  if (statistic == "pvalue") {
    if (max(x) > 1 | min(x) < 0)
      stop("input p-values must all be in the range 0 to 1!")
  }
  if (verbose)
    cat("Step 1... determine cutoff point\n")
  if (cutoff.method == "pct0") {
    if (statistic == "pvalue")
      x0 = quantile(x_red, probs = 1 - pct0)
    else x0 = quantile(abs(x_red), probs = pct0)

  }  else if (cutoff.method == "locfdr" & (statistic == "normal" |
                                        statistic == "correlation")) {
    if (statistic == "normal")
      z = x_red
    if (statistic == "correlation")
      z = atanh(x_red)
    iqr = as.double(diff(quantile(z, probs = c(0.25, 0.75))))
    sdhat = iqr/(2 * qnorm(0.75))
    N = length(x)
    b = ifelse(N > 5e+05, 1, 4.3 * exp(-0.26 * log(N, 10)))
    z0 = b * sdhat
    if (statistic == "normal")
      x0 = z0
    if (statistic == "correlation")
      x0 = tanh(z0)

  }  else {
    if (cutoff.method == "locfdr")
      warning("cutoff.method=\"locfdr\" only available for normal and correlation statistic.")
    x0 = fdrtool::fndr.cutoff(x_red, statistic) #-> uitgevoerd
  }
  if (verbose)
    cat("Step 2... estimate parameters of null distribution and eta0\n")

  ###Alle x0 moet met aangepaste x'n
  ###Dit moet met aangepaste x'n:###
  cf.out <- fdrtool::censored.fit(x = x_red, cutoff = x0, statistic = statistic) #-> uitgevoerd
  ##############################

  if (statistic == "pvalue")
    scale.param = NULL
  else scale.param <- cf.out[1, 5]
  eta0 = cf.out[1, 3]
  if (verbose)
    cat("Step 3... compute p-values and estimate empirical PDF/CDF\n")
  nm = fdrtool:::get.nullmodel(statistic)

  ###Dit moet met aangepaste pval berekend worden:###
  pval_full = nm$get.pval(x, scale.param) #scale.param is the sd that will be used, set to 1 to get normal p values
  pval = nm$get.pval(x_red, scale.param)
  ee <- fdrtool:::ecdf.pval(pval, eta0 = eta0) #empirical cdf
  ##############################

  g.pval <- fdrtool::grenander(ee)
  f.pval = approxfun(g.pval$x.knots, g.pval$f.knots, method = "constant",
                     rule = 2)
  f0.pval = function(x) return(ifelse(x > 1 | x < 0, 0, rep(1,
                                                            length(x))))
  F.pval = approxfun(g.pval$x.knots, g.pval$F.knots, method = "linear",
                     yleft = 0, yright = g.pval$F.knots[length(g.pval$F.knots)])
  F0.pval = function(x) return(ifelse(x > 1, 1, ifelse(x <
                                                         0, 0, x)))
  fdr.pval = function(p, eta0) { #eta0 added as argument
    p[p == .Machine$double.eps] = 0
    pmin(eta0/f.pval(p), 1)
  }
  Fdr.pval = function(p, eta0) pmin(eta0 * p/F.pval(p), 1) #eta0 added as argument
  if (verbose)
    cat("Step 4... compute q-values and local fdr\n")

  eta0_full <- (eta0*length(x_red)+(length(x)-length(x_red)))/length(x)

  qval_full <- Fdr.pval(pval_full, eta0=eta0_full)
  lfdr_full <- fdr.pval(pval_full, eta0=eta0_full)

  ###Should be equal!!!
  #For some reason, they are not entirely equal if you go through the functions manually, but they are if you execute the functions
  #plot(fdrtool::fdrtool(x_red, plot=FALSE)$lfdr)
  #plot(fdr.pval(pval, eta0=eta0))

  ###Exactly eta0.SE (remains as is because completely calculated based on x_red and only adjusted to account for difference in lenght):
  #Adjustments would make it smaller, would not be correct
  # m <- 1-nm$get.pval(x0, scale.param)
  # N.cens <- cf.out[1,2]
  # N <- length(x_red)
  # th <- N.cens/N
  #
  # sqrt(th * (1 - th)/(N * m * m))
  ###############

  ###Extra functions only needed for plotting###
  if (statistic == "pvalue") {
    f0 <- function(zeta) return(nm$f0(zeta, scale.param))
    F0 <- function(zeta) return(nm$F0(zeta, scale.param))
    get.pval <- function(zeta) return(nm$get.pval(1 -
                                                    zeta, scale.param))
    x0 = 1 - x0
  }  else {
    f0 <- function(zeta) return(2 * nm$f0(zeta, scale.param))
    F0 <- function(zeta) return(2 * nm$F0(zeta, scale.param) - 1)
    get.pval <- function(zeta) return(nm$get.pval(zeta,
                                                  scale.param))
  }

  fdr = function(zeta) fdr.pval(get.pval(zeta), eta0_full)
  Fdr = function(zeta) Fdr.pval(get.pval(zeta), eta0_full)
  ########

  result = list(x=x, x_red=x_red, pval = pval_full, qval = qval_full, lfdr = lfdr_full, statistic = statistic,
                param = cf.out, x0 = x0, f.pval = f.pval, F.pval = F.pval, fdr.pval=fdr.pval, Fdr.pval=Fdr.pval, f0=f0, F0=F0, get.pval=get.pval, fdr=fdr, Fdr=Fdr)
  if (verbose)
    cat("\n")
  }

  return(result)
}


#' Make diagnostic plots for fdrtool results
#'
#' This function allows to make the three diagnostic plots from the \pkg{fdrtool} package by Korbinian Strimmer (\url{http://strimmerlab.org}) in separate windows.
#' The implementation is identical to that in the \pkg{fdrtool} package, except that it allows to plot the functions separately.
#' The first plot is also different in showing a histogram for both the vector of the observed test statistics \code{x} in white
#' and a histogram for the subset of this vector that was used in the \code{\link{fdrtool_subset}} function in gray.
#' @param result The result of a call to the \code{\link{fdrtool_subset}} function.
#' @param make_plots A vector of three logical values. Each value indicates whether a certain plot should be made. Defaults to \code{c(TRUE, TRUE, TRUE)} indicating that all three plots should be made.
#' @param color.figure A logical indicating whether a color figure or a black and white figure is produced. Defaults to \code{TRUE}, i.e. to color figure.
#' @param verbose A logical indicating if status messages should be printed out. Defaults to \code{TRUE}.
#' @return Diagnostic plots.
#' @examples .......
#' @references Strimmer, K. (2008a). A unified approach to false discovery rate estimation. BMC Bioinformatics 9: 303. Available from \url{http://www.biomedcentral.com/1471-2105/9/303/}.
#'
#' Strimmer, K. (2008b). fdrtool: a versatile R package for estimating local and tail area- based false discovery rates. Bioinformatics 24: 1461-1462. Available from \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/12/1461}.
#' @export
plot_fdrtool <- function(result, make_plots = c(TRUE, TRUE, TRUE), color.figure = TRUE, verbose = TRUE) {

  ###Added###
  x = result$x
  x_red = result$x_red
  x0 = result$x0
  param = result$param
  statistic = result$statistic
  nm = fdrtool:::get.nullmodel(statistic)

  eta0 <- param[1, 3]
  eta0_full <- (eta0*length(x_red)+(length(x)-length(x_red)))/length(x)

  if (statistic == "pvalue")
    scale.param = NULL
  else scale.param <- param[1, 5]

  if (verbose) cat("Prepare for plotting\n")

  ax = abs(x)
  ax_red = abs(x_red)
  if (statistic == "pvalue") {ax = 1 - ax
  ax_red = 1 - ax_red
  }
  xxx = seq(0, max(ax[!is.infinite(ax)]), length.out = 500)
  ll = fdrtool:::pvt.plotlabels(statistic, scale.param, eta0_full)

  if (color.figure)
    cols = c(2, 4)
  else cols = c(1, 1)

  f.pval = result$f.pval #approxfun
  F.pval = result$F.pval #approxfun
  fdr.pval = result$fdr.pval #needs p and eta0_full
  Fdr.pval = result$Fdr.pval #needs p and eta0_full
  f0 = result$f0 #needs zeta, nm and scale.param
  F0 = result$F0 #needs zeta, nm and scale.param
  get.pval = result$get.pval #needs zeta and scale.param
  fdr = result$fdr #needs fdr.pval, get.pval, zeta and eta0_full
  Fdr = result$Fdr #needs Fdr.pval, get.pval, zeta and eta0_full

  ####################


  #Plot 1

  if(isTRUE(make_plots[1])){
  f = function(zeta) eta0_full * (f0(zeta))/fdr(zeta)
  fA = function(zeta) (f(zeta) - eta0_full * f0(zeta))/(1 -
                                                          eta0_full)

  hist(ax, freq = FALSE, bre = 50, main = ll$main, xlab = ll$xlab,
       cex.main = 1.8)
  hist(ax_red, freq = FALSE, bre = 50, col="gray", add=TRUE)
  lines(xxx, eta0_full * f0(xxx), col = cols[1], lwd = 2, lty = 3)
  lines(xxx, (1 - eta0_full) * fA(xxx), col = cols[2], lwd = 2)
  if (statistic == "pvalue")
    pos1 = "topleft"
  else pos1 = "topright"
  legend(pos1, c("Mixture", "Null Component", "Alternative Component"),
         lwd = c(1, 2, 2), col = c(1, cols), lty = c(1, 3,
                                                     1), bty = "n", cex = 1.5)
  }


  #Plot 2

  if(isTRUE(make_plots[2])){
  F = function(zeta) 1 - eta0_full * get.pval(zeta)/Fdr(zeta)
  FA = function(zeta) (F(zeta) - eta0_full * F0(zeta))/(1 -
                                                          eta0_full)

  plot(xxx, F(xxx), lwd = 1, type = "l", ylim = c(0, 1),
       main = "Density (first row) and Distribution Function (second row)",
       xlab = ll$xlab, ylab = "CDF", cex.main = 1.5)
  lines(xxx, eta0_full * F0(xxx), col = cols[1], lwd = 2, lty = 3)
  lines(xxx, (1 - eta0_full) * FA(xxx), col = cols[2], lwd = 2)
  }


  #Plot 3

  if(isTRUE(make_plots[3])){
  plot(xxx, Fdr(xxx), type = "l", lwd = 2, ylim = c(0,
                                                    1), main = "(Local) False Discovery Rate", ylab = "Fdr and fdr",
       xlab = ll$xlab, lty = 3, cex.main = 1.5)
  lines(xxx, fdr(xxx), lwd = 2)
  if (eta0_full > 0.98)
    pos2 = "bottomleft"
  else pos2 = "topright"
  legend(pos2, c("fdr (density-based)", "Fdr (tail area-based)"),
         lwd = c(2, 2), lty = c(1, 3), bty = "n", cex = 1.5)
  }

  rm(ax)
}






