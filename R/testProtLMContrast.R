#' Test a contrast
#'
#' @description This function can test a contrast based on a protLM object \code{protLM} and a contrast matrix \code{L}.
#' @param protLM An object of class \code{\link[=protLM-class]{protLM}}.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param add.annotations A logical indicating whether the \code{annotations} slot of the \code{\link[=protLM-class]{protLM}} object should be added as extra columns to each matrix in the returned list of matrices. Defaults to \code{TRUE}.
#' @param custom_dfs An optional vector of length equal to the number of models in \code{protLM} containing the degrees of freedom that should be used for each model in the \code{protLM} object. By default (\code{custom_dfs=NULL}), \code{MSqRob} uses the trace of the Hat matrix to determine the degrees of freedom or the Satterhtwiate approximation (if \code{satterthwaite=TRUE}). \code{custom_dfs} allows you to overrule this and use your own custom approximation for the degrees of freedom.
#' @param simplify A logical indicating wheter, if there is only one contrast, a matrix should be returned instead of a list containing one matrix. Defaults to \code{TRUE}.
#' @param lfc The minimum (log2) fold-change that is considered scientifically meaningful. Defaults to \code{0}. Ignored when \code{anova = TRUE}.
#' @param anova A logical indicating whether the contrasts should be tested simultaneously in one F-test (\code{anova = TRUE}) or as separate t-tests (\code{anova = FALSE}). Defaults to \code{FALSE}.
#' @param anova.na.ignore A logical indicating whether contrasts that cannot be fitted due to missing observations should be ignored when calculating an F value (\code{anova.na.ignore = TRUE}) or NA should be returned when at least one contrast cannot be estimated (\code{anova.na.ignore = FALSE}). Defaults to \code{TRUE}. Ignored when \code{anova = FALSE}.
#' @param exp_unit The effect that in all models corresponds to the experimental unit. Only needed when one would like to calculate a more conservative way of estimating the degrees of freedom.
#' The default way of estimating the degrees of freedom (\code{exp_unit=NULL}) subtracts the total number of observations by the trace of the Hat matrix. However, often, observations are not completely independent. A more conservative way (\code{df_exp}) is defining on which level the treatments were executed and substracting all degrees of freedom lost due to between-treatement effects (\code{pars_df}) from the number of treatments.
#' @param pars_df Only used if exp_unit is not \code{NULL}. Character vector indicating all parameters in the models that are between-treatment effects in order to calculate a more conservative degrees of freedom (\code{df_exp}). If left to default (\code{NULL}), all parameters in the models will be asumed to be between-treatment effects (this is not adviced as the result will mostly be too conservative).
#' @param satterthwaite A logical indicating whether the Satterthwaite approximation for the degrees of freedom should be used instead of the classical trace of the Hat matrix. Defaults to \code{FALSE}.
#' @param lmerModFun Only used when \code{satterthwaite=TRUE}. \code{lmerModFun} indicates which deviance function should be used when calculating the Satterthwaite approximation for the degrees of freedom. The default (\code{NULL}) uses the lme4 \code{\link[=lme4::mkLmerDevfun]{mkLmerDevfun}} function to generate the deviance function. This parameter should only rarely, if ever, be changed.
#' @param gradMethod Only used when \code{satterthwaite=TRUE}. One of "Richardson", "simple", or "complex" indicating the method to use for the gradient calculation by numerical approximation during the calculation of the Satterthwaite approximation for the degrees of freedom. Defaults to "simple".
#' @param printProgress A logical indicating whether the R should print a message before calculating the contrasts for each accession. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message_extract Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the extraction of beta, vcov, df and sigma, or \code{NULL} to hide the current message (if any).
#' @param message_test Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the testing of the contrasts, or \code{NULL} to hide the current message (if any).
#' @return A list of data frames, with each data frame in the list corresponding to a contrast in L. Each row of the data frame corresponds to a protein in the \code{\link[=protLM-class]{protLM}} object.
#' The \code{estimate} column contains the size estimate of the contrast, the \code{se} column contains the estimated standard error on the contrast, the \code{Tval} column contains the T-value corresponding to the contrast and the \code{pval} column holds the p-value corresponding to the contrast.
#' If \code{simplify=TRUE} and the \code{\link[=protLM-class]{protLM}} object contains only one element, the data frame is not present in a list.
#' @include protdata.R
#' @include updateProgress.R
#' @include protLM.R
#' @include getBetaVcovDf.R
#' @include getBetaVcovDfList.R
#' @include testContrast.R
#' @include squeezePars.R
#' @export
test.protLMcontrast <- function(protLM, L, add.annotations=TRUE, custom_dfs=NULL, simplify=TRUE, lfc=0, anova=FALSE, anova.na.ignore=TRUE, exp_unit=NULL, pars_df=NULL, satterthwaite=FALSE, lmerModFun=NULL, gradMethod="simple", printProgress=FALSE, shiny=FALSE, message_extract=NULL, message_test=NULL)
{

  if(is.null(rownames(L))) stop("L should be a matrix with row names corresponding to model predictors.")

  #if(!all(rownames(L) %in% rownames(betaVcovDf$beta))) stop(paste0("\"",rownames(L)[!(rownames(L) %in% rownames(betaVcovDf$beta))],"\" does not correspond to any model predictor."))

  betaVcovDfList <- getBetaVcovDfList(protLM, exp_unit=exp_unit, pars_df=pars_df, printProgress=printProgress, shiny=shiny, message=message_extract)

  #Progress bar for testing contrasts
  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message_test, value = 0)
  }

  betas <- lapply(betaVcovDfList, function(x) {return(x$beta)})
  vcovs <- lapply(betaVcovDfList, function(x) {return(x$vcov)})
  sigmas <- sapply(betaVcovDfList, function(x) {return(x$sigma)})
  dfs <- sapply(betaVcovDfList, function(x) {return(x$df)})
  #Must be lapply, because if all NULL, makes list of nulls with sapply!
  dfs_exp <- lapply(betaVcovDfList, function(x) {return(x$df_exp)})

  #List of models:
  models <- getModels(protLM, simplify=FALSE)
  #List of proteins:
  proteins <- getAccessions(protLM)
  #Types of models: either "lmerMod" or "lm":
  classes <- sapply(models, "class")
  #Matrix of annotations
  annotation_matrix <- getAnnotations(protLM)

  if(!isTRUE(add.annotations)){annotation_matrix <- matrix(nrow=nrow(annotation_matrix),ncol=0)}

  if(!is.null(custom_dfs)){dfs <- custom_dfs}

  #Initialize coefmatlist

  #If it is NOT an ANOVA
  if(!isTRUE(anova)){
  coefmatlist <- rep(list(cbind.data.frame(annotation_matrix,matrix(rep(NA,length(proteins)*5), nrow=length(proteins), ncol=5, dimnames=list(proteins,c("estimate","se","df","Tval","pval"))))), ncol(L))
  names(coefmatlist) <- colnames(L)

  #If it is an ANOVA
  } else{  coefmatlist <- list(cbind.data.frame(annotation_matrix,matrix(rep(NA,length(proteins)*5), nrow=length(proteins), ncol=5, dimnames=list(proteins,c("AveExpr", "df_num", "df_den", "Fval", "pval"))))) }

  for(i in 1:length(models))
  {

    updateProgress(progress=progress, detail=paste0("Testing contrast(s) for protein ",i," of ",length(models),"."), n=length(models), shiny=shiny, print=isTRUE(printProgress))

    beta <- betas[[i]]
    vcov <- vcovs[[i]]
    sigma <- sigmas[i]
    df <- dfs[i]
    df_exp <- dfs_exp[[i]]
    model <- models[[i]]

   if(!(classes[i] %in% c("lm","lmerMod"))){

     if(!isTRUE(anova)){
      contrasts <- matrix(rep(NA,5), nrow=ncol(L), ncol=5, dimnames=list(rep(proteins[i],ncol(L)),c("estimate","se","df","Tval","pval")))
     } else{
      contrasts <- matrix(rep(NA,5), nrow=ncol(L), ncol=5, dimnames=list(rep(proteins[i],ncol(L)),c("AveExpr", "df_num", "df_den", "Fval", "pval")))
     }
      warning(paste0("Unrecognized class for model ",i,". Models should be either of type \"lm\" or \"lmerMod\". Returning NA."))

    #if(classes[i] %in% c("lm","lmerMod")):
      } else {

      if(!all(is.na(beta))){beta <- na.omit(beta)}
      betanames <- rownames(beta)

      xlevels <- getXlevels(model,classes[i])

      L2 <- makeContrastMatrix(L, betanames, xlevels)

    #If satterthwaite is TRUE, we use the Satterthwaite approximation for df
    if(isTRUE(satterthwaite)){
    df <- dfSatterthwaite(lmerModFun, model, vcov, sigma, L2, gradMethod)

    #Change df to df_exp if you use the more conservative dfs based on the experimental units
    #Note: you should use the normal dfs for the squeezing of the variances!!!
    } else if(!is.null(df_exp)){df <- df_exp}


    if(!isTRUE(anova)){
    contrasts <- test.contrast(beta, vcov, df, sigma, L2, lfc=lfc)
    #if(isTRUE(anova)):
    } else{
    contrasts <- test.ANOVA(beta, vcov, df, sigma, L2, anova.na.ignore=anova.na.ignore)
    }


    rownames(contrasts) <- rep(proteins[i],nrow(contrasts))

    }

    #Add each contrast to its corresponding matrix in the list
    for(j in 1:length(coefmatlist))
    {
      coefmatlist[[j]][i,(ncol(annotation_matrix)+1):(ncol(annotation_matrix)+5)] <- contrasts[j,,drop=FALSE]
    }

  }

  #If there is only 1 contrast and simplify==TRUE, then return a matrix instead of a list with a matrix
  if(simplify & length(coefmatlist)==1)
  {
    coefmatlist <- coefmatlist[[1]]
  }

  return(coefmatlist)
}


makeContrastMatrix=function(L, betanames, xlevels){
a <- matrix(!(rownames(L) %in% xlevels), nrow=nrow(L), ncol=ncol(L))
b <- L!=0
#If a parameter is not estimated in the model, but is present with a non-zero value in a contrast, this value should be NA:
L[a&b] <- NA

L2 <- matrix(0, nrow=length(betanames), ncol=ncol(L), dimnames=list(betanames,colnames(L)))
welke <- rownames(L)[rownames(L) %in% rownames(L2)]
L2[welke,] <- L[welke,]
#If any value in a contrast cannot be estimated, the whole contrast value should be NA:
L2[,is.na(colSums(L))] <- NA
return(L2)
}

getXlevels=function(model, class){
  if(class=="lm"){

    xlevels <- vector()
    for(k in 1:length(model$xlevels))
    {
      xlevels <- c(xlevels,paste0(names(model$xlevels[k]),model$xlevels[[k]]))
    }
    if("(Intercept)" %in% names(coef(model)))
    {xlevels <- c("(Intercept)",xlevels)}

    #!!!Uitgebreid testen met gevallen waarin bv. A ontbreekt, A en B ontbreken, A en C, C en D,...

    #Een keuze die we maken: standaard gooit vcov er de NA-parameters uit, wij doen dit ook voor de betas
    #Alternatief: van de moment dat er een NA is, alles op NA zetten; nadeel: sommige modellen werken dan niet terwijl als je 1 observatie zou weglaten, ze wel zouden werken...
    #lmer does this automatically with its unshrunken fixed effects

  } else if(class=="lmerMod"){
    xlevels <- unlist(attr(model,"MSqRob_levels"))
  }
}

dfSatterthwaite=function(lmerModFun, model, vcov, sigma, L2, gradMethod){

  #This function is only used for df calculation for inference => we use the vcov that is calculated back to the original scale!

  #extract useful elements from lmer model
  #CalculateHessian
  if (is.null(lmerModFun)) lmerModFun=devFun(model)
  environment(lmerModFun)$pars=""
  sigmaTheta=c(sigma,model@pp$theta)
  h=numDeriv::hessian(function(x) devianceForHessian(x,lmerModFun),x=sigmaTheta)

  #Calculate variance covariance for variance parameters, we use the Schnabel/Eskow Cholesky factorization when the hessian is non-invertable
  vcovTheta=tryCatch(2*chol2inv(sechol(h)), error=function(e)-matrix(data=NA,nrow=ncol(h),ncol=ncol(h)))

  #Calculate variance covariance for BLUP
  #vcov in these functions is equal to our vcov*sigma^2
  vcov=vcov*sigma^2

  #Calculate variance of contrast
  vcovContrast=t(L2)%*%vcov%*%L2

  #Create gew: the sqrt of the weigths
  if(is.null(model@frame$"(weights)")){gew <- rep(1,nobs(model))} else{
    gew <- sqrt(model@frame$"(weights)")}


  df <- sapply(1:ncol(L2),function(ic,L2,vcovContrast)
  {
    Lh=as.matrix(L2[,ic])
    #CalculateGradient
    g=numDeriv::grad(function(x) {

      return(as.matrix(t(Lh)%*%getVcovBetaBforGrad(x,lmerModFun, R_fix=attr(model,"MSqRob_R_fix"), Zt_indices=attr(model, "MSqRob_Zt_indices"), unshr_pos=attr(model, "MSqRob_unshr_pos"), gew=gew)%*%Lh))
    } , x=sigmaTheta,method=gradMethod)
    #Calculate denominator for sattertwaite approximation of degrees of freedom
    denom <- t(g) %*% vcovTheta %*% g
    #Calculate degrees of freedom
    return(as.double(2*(vcovContrast[ic,ic])^2/denom))
  },L=L2,vcovContrast=vcovContrast)

  return(df)

}


# devFun=function(lmerMod)
# {
#   return(update(lmerMod,devFunOnly=TRUE))
# }
#outputFun= function(theta) .Call(MSqRob_lmer_Deviance,pp$ptr(),resp$ptr(),as.double(theta)))
#die moet nog in environment draaien met lower, pp, resp en MSqRob_lmer_Deviance
#return(outputFun)
devFun=function(lmerMod)
{
  reTrmsHlp=list(Zt=lmerMod@pp$Zt,theta= lmerMod@pp$theta,Lambdat= lmerMod@pp$Lambdat, Lind=lmerMod@pp$Lind,Gp=lmerMod@Gp,lmerMod@lower,flist=lmerMod@flist, cnms=lmerMod@cnms)
  return(lme4::mkLmerDevfun(lmerMod@frame, lmerMod@pp$X, reTrmsHlp, REML = TRUE, start = NULL))
}




getVcovBetaBforGrad<-function(sigmaTheta,devFunUpdate,Ginvoffset = 1e-18, R_fix, Zt_indices, unshr_pos, gew)
{
  sigma <- sigmaTheta[1]
  theta <- sigmaTheta[-1]
  devUpdate=devFunUpdate(theta)
  rho=environment(devFunUpdate)
  p=ncol(rho$pp$X)
  q=nrow(rho$pp$Zt)
  Ct=rbind2(t(rho$pp$X),rho$pp$Zt)
  #B, the inverse of the estimated variance-covariance matrix of the random effects. We add a small offset Ginvoffset to the diagonal of B to prevent near-singularity:
  lambdat <- rho$pp$Lambdat
  B=Matrix::solve(Matrix::crossprod(lambdat)+Matrix::Diagonal(q,Ginvoffset))

  #Sqrt-weighted design matrix Ct_w (sparse matrix):
  Ct_w <- Matrix::tcrossprod(Ct,Matrix::Diagonal(length(gew))*gew)

  vcovInv <- Matrix::tcrossprod(Ct_w)
  vcovInv[((p+1):(q+p)),((p+1):(q+p))]=vcovInv[((p+1):(q+p)),((p+1):(q+p))]+B

  if(!is.null(R_fix)){

    R_fix_large <- Matrix::Matrix(diag(nrow(vcovInv)), sparse=TRUE)
    R_fix_large[c(1:p,Zt_indices+p),c(1:p,Zt_indices+p)] <- R_fix
    # R_fix_large[1:p,1:p] <- R_fix[unshr_pos,unshr_pos]
    # R_fix_large[Zt_indices+p,Zt_indices+p] <- R_fix[-unshr_pos,-unshr_pos]
    # R_fix_large[1:p,Zt_indices+p] <- R_fix[unshr_pos,-unshr_pos]
    # R_fix_large[Zt_indices+p,1:p] <- R_fix[-unshr_pos,unshr_pos]

    colnames(R_fix_large) <- colnames(vcovInv)
    rownames(R_fix_large) <- rownames(vcovInv)

    #Put the shrunken fixed effects back to their original scale:
    vcovInv <-Matrix::t(R_fix_large)%*%vcovInv%*%R_fix_large
    #Ct[Zt_indices+p,] <- R_fix%*%Ct[Zt_indices+p,]
    #Ct_w[Zt_indices+p,] <- R_fix%*%Ct_w[Zt_indices+p,]
  }

  return(chol2inv(sechol(as.matrix(vcovInv)))*sigma^2)
}



devianceForHessian <- function(sigmaTheta,devFunUpdate) {

  sigma <- sigmaTheta[1]
  theta <- sigmaTheta[-1]
  devtmp=devFunUpdate(theta)
  rho=environment(devFunUpdate)

  if(!is.null(rho$getCovBlocks)){
    beta <- rho$pp$beta(1)
    Lambda.ts <- rho$getCovBlocks(rho$pp$Lambdat, rho$ranefStructure)
    exponentialTerms <- rho$calculatePriorExponentialTerms(rho$priors,
                                                           beta, Lambda.ts)
    polynomialTerm <- rho$calculatePriorPolynomialTerm(rho$priors$covPriors,
                                                       Lambda.ts)
    return(rho$lmmObjective(rho$pp, rho$resp, sigma, exponentialTerms, polynomialTerm,
                            rho$blmerControl))
  }
  else{
    prss=rho$pp$sqrL(1)+rho$resp$wrss()
    n <- ncol(rho$pp$Zt)
    dev <- n * log(2 * pi * sigma^2) + rho$pp$ldL2() + prss/sigma^2 +
      rho$pp$ldRX2()* (rho$resp$REML>0) -
      rho$resp$REML * log(2 * pi * sigma^2)-sum(log(rho$resp$weights))

    return(as.vector(dev))
  }
}




