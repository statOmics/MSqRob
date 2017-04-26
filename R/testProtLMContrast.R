#' Test a contrast
#'
#' @description This function can test a contrast based on a protLM object \code{protLM} and a contrast matrix \code{L}.
#' @param protLM An object of class \code{\link[=protLM-class]{protLM}}.
#' @param L A contrast matrix with the parameter levels as rows and a column for each contrast.
#' @param add.annotations A logical indicating whether the \code{annotations} slot of the \code{\link[=protLM-class]{protLM}} object should be added as extra columns to each matrix in the returned list of matrices. Defaults to \code{TRUE}.
#' @param simplify A logical indicating wheter, if there is only one contrast, a matrix should be returned instead of a list containing one matrix. Defaults to \code{TRUE}.
#' @param lfc The minimum (log2) fold-change that is considered scientifically meaningful. Defaults to \code{0}. Ignored when \code{anova = TRUE}.
#' @param anova A logical indicating whether the contrasts should be tested simultaneously in one F-test (\code{anova = TRUE}) or as separate t-tests (\code{anova = FALSE}). Defaults to \code{FALSE}.
#' @param anova.na.ignore A logical indicating whether contrasts that cannot be fitted due to missing observations should be ignored when calculating an F value (\code{anova.na.ignore = TRUE}) or NA should be returned when at least one contrast cannot be estimated (\code{anova.na.ignore = FALSE}). Defaults to \code{TRUE}. Ignored when \code{anova = FALSE}.
#' @param type_dfs Either one of \code{"residual"}, \code{"between-within"}, \code{"Satterthwaite"}, \code{"exp_between"} or \code{"custom"}. This argument indicates how the degrees of freedom should be calculated. Defaults to \code{"residual"}.
#' More information is given under 'Details'.
#' @param custom_dfs Only used if \code{type_dfs="custom"}. A list of length equal to the number of models in \code{protLM} containing vectors of lenght \code{ncol(L)} (if \code{anova = FALSE}) representing the degrees of freedom that should be used for each contrast in \code{L} and each model in the \code{protLM} object. If \code{anova = TRUE}, each element of the list should contain a single numeric value representing the degrees of freedom to be used for the anova test. Defaults to \code{NULL}.
#' @param exp_unit Only used if \code{type_dfs="exp_between"}. The effect that in all models corresponds to the experimental unit.
#' @param pars_between Only used if \code{type_dfs="exp_between"}. Character vector indicating all parameters in the models that are between-treatment effects. If left to default (\code{NULL}), all parameters in the models will be asumed to be between-treatment effects (this is not adviced as the result will mostly be too conservative).
#' @param lmerModFun Only used when \code{satterthwaite=TRUE}. \code{lmerModFun} indicates which deviance function should be used when calculating the Satterthwaite approximation for the degrees of freedom. The default (\code{NULL}) uses the lme4 \code{\link[=lme4::mkLmerDevfun]{mkLmerDevfun}} function to generate the deviance function. This parameter should only rarely, if ever, be changed.
#' @param gradMethod Only used when \code{satterthwaite=TRUE}. One of "Richardson", "simple", or "complex" indicating the method to use for the gradient calculation by numerical approximation during the calculation of the Satterthwaite approximation for the degrees of freedom. Defaults to "simple".
#' @param printProgress A logical indicating whether the R should print a message before calculating the contrasts for each accession. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message_extract Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the extraction of beta, vcov, df and sigma, or \code{NULL} to hide the current message (if any).
#' @param message_test Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the testing of the contrasts, or \code{NULL} to hide the current message (if any).
#' @details Calculating degrees of freedom (and hence p values) for mixed models with unbalanced designs is an unresolved issue in the field (see for example here https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html and here https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q2/000904.html).
#' We offer different approximations and leave it up to the user to select his/her preferred approach.
#' \code{"residual"} calculates approximative degrees of freedom by subtracting the trace of the hat matrix from the number of observations. It is the default setting, but this approach might be somewhat too liberal.
#' \code{"Satterthwaite"} calculates approximative degrees of freedom using Satterthwaite's approximation (Satterthwaite, 1946). This approximative approach is used in many applications but is rather slow to calculate and might lead to some missing values due difficulties in calculating the Hessian.
#' \code{"exp_between"} calculates approximative degrees of freedom by defining on which level the treatments were executed and substracting all degrees of freedom lost due to between-treatement effects (\code{pars_between}) from the number of experimental units \code{exp_unit}. This allows to mimick the behaviour of \code{type_dfs="between-within"} for more complex designs.
#' \code{"custom"} Allows the user to provide his/her own degrees of freedom for each contrast and each protein. Custom degrees of freedom should be entered in the \code{custom_dfs} field.
#' @return A list of data frames, with each data frame in the list corresponding to a contrast in \code{L}. Each row of the data frame corresponds to a protein in the \code{\link[=protLM-class]{protLM}} object.
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
test.protLMcontrast <- function(protLM, L, add.annotations=TRUE, simplify=TRUE, lfc=0, anova=FALSE, anova.na.ignore=TRUE, type_dfs="residual", custom_dfs=NULL, exp_unit=NULL, pars_between=NULL, lmerModFun=NULL, gradMethod="simple", printProgress=FALSE, shiny=FALSE, message_extract=NULL, message_test=NULL)
{

  if(is.null(rownames(L))) stop("L should be a matrix with row names corresponding to model predictors.")

  #if(!all(rownames(L) %in% rownames(betaVcovDf$beta))) stop(paste0("\"",rownames(L)[!(rownames(L) %in% rownames(betaVcovDf$beta))],"\" does not correspond to any model predictor."))

  betaVcovDfList <- getBetaVcovDfList(protLM, exp_unit=exp_unit, pars_between=pars_between, printProgress=printProgress, shiny=shiny, message=message_extract)

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
  dfs <- lapply(betaVcovDfList, function(x) {return(x$df)})
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
    df <- dfs[[i]] #default: type_dfs=="residual" Other options: see below
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
    if(type_dfs=="Satterthwaite"){
      df <- dfSatterthwaite(lmerModFun, model, vcov, sigma, L2, gradMethod, anova=anova)
    #Change df to df_exp if you use the more conservative dfs based on the experimental units
    #Note: you should use the normal dfs for the squeezing of the variances!!!
    } else if(type_dfs=="exp_between"){
      df <- df_exp
    } else if(type_dfs=="custom"){
      df <- custom_dfs[[i]]
    }

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

  return(xlevels)
}

########################################################################
##########Satterthwaite degrees of freedom calculation##################
########################################################################

dfSatterthwaite=function(lmerModFun, model, vcov, sigma, L2, gradMethod, anova=FALSE){

  if(isTRUE(anova)){stop("Satterthwaire degrees of freedom not yet implemented for ANOVA.")}

  #This function is only used for df calculation for inference => we use the vcov that is calculated back to the original scale!

  #extract useful elements from lmer model
  #CalculateHessian
  if (is.null(lmerModFun)) lmerModFun=devFun(model)
  environment(lmerModFun)$pars=""
  sigmaTheta=c(sigma,model@pp$theta) #sigmaTheta=c(sigma,model@pp$theta[model@pp$theta!=0])
  h=numDeriv::hessian(function(x) devianceForHessian(x,lmerModFun),x=sigmaTheta)

  #Calculate variance covariance for variance parameters, we use the Schnabel/Eskow Cholesky factorization when the hessian is non-invertable
  vcovTheta=tryCatch(2*solve(h), error=function(e)-matrix(data=NA,nrow=ncol(h),ncol=ncol(h)))

  #old:
  #chol2inv(sechol(h))

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

  ## IMPORTANT! Deep-copy the (only) reference-class slots! ...
  lmerMod3 <- lmerMod
  lmerMod3@pp <- lmerMod@pp$copy()
  lmerMod3@resp <- lmerMod@resp$copy()

  # lmerMod2@cnms <- lmerMod2@cnms[which(lmerMod2@theta!=0)]
  # lmerMod2@flist <- lmerMod2@flist[which(lmerMod2@theta!=0)]
  # lmerMod2@theta <- lmerMod2@theta[which(lmerMod2@theta!=0)]
  # # object2@pp$setTheta(lmerMod2@theta[which(lmerMod2@theta!=0)])
  #
  # dd <- as.function(lmerMod2)
  # dval <- dd(lmerMod2@theta[which(lmerMod2@theta!=0)])
  #
  # lmerMod2@pp$setTheta(lmerMod2@theta[which(lmerMod2@theta!=0)])
  #
  # lmerMod3 <- mkMerMod(environment(dd),
  #                           opt=list(par=theta,fval=dval,conv=0),
  #                           fr=model.frame(lmerMod2),
  #                           reTrms=getME(lmerMod2,
  #                                        c("Zt","Lambdat","Lind","theta",
  #                                          "lower","flist","cnms","Gp")))

  reTrmsHlp=list(Zt=lmerMod3@pp$Zt,theta= lmerMod3@pp$theta,Lambdat= lmerMod3@pp$Lambdat, Lind=lmerMod3@pp$Lind,Gp=lmerMod3@Gp,lmerMod3@lower,flist=lmerMod3@flist, cnms=lmerMod3@cnms)
  return(lme4::mkLmerDevfun(lmerMod3@frame, lmerMod3@pp$X, reTrmsHlp, REML = TRUE, start = NULL))
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

  #Old:
  #return(chol2inv(sechol(as.matrix(vcovInv)))*sigma^2)

  #Estimated variance-covariance matrix vcov_sigma:
  vcov_sigma <- tryCatch(as.matrix(Matrix::solve(vcovInv))*sigma^2, error=function(e){
    return(vcovInv*NA)
  })

  rownames(vcov_sigma) <- colnames(vcovInv)
  colnames(vcov_sigma) <- rownames(vcovInv)

  return(vcov_sigma)
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


########################################################################
##########End Satterthwaite degrees of freedom calculation##############
########################################################################

########################################################################
##########Between-within degrees of freedom calculation#################
########################################################################

#Should be done via custom_dfs because we also need the proteins as input!

dfBetweenWithin <- function(protein, model, L2, anova){

  #De zaken die je van protein moet hebben als attribute meegeven!!!!

  if(isTRUE(anova)){stop("Between-within degrees of freedom not yet implemented for ANOVA.")}

  #vertrekken van proteinS, modelS

  comp_preds <- getComp_preds(model, L2)
  options <- names(attr(model,"MSqRob_cnms"))

  ns <- getN(protein, options, comp_preds, L2)

  ps <- getPs(protein, ns, options, comp_preds, L2)

  elements <- rep(FALSE,length(names(getBetaVcovDf(model)$df_pars)))

  #Escape special characters such as the "(" in "(Intercept)"
  library(stringr)
  quotemeta <- function(string) {
    str_replace_all(string, "(\\W)", "\\\\\\1")
  }

  for(i in 1:length(ps[[1]])){
    elements <- elements | grepl(paste0("^",quotemeta(ps[[1]][i])),names(getBetaVcovDf(model)$df_pars))
  }

  ####Probleem: groter dan 6... -> je kan ook niet zomaar de treatments gaan nemen die over je contrast gaan,
  #want is reference class codering...

  getBetaVcovDf(model)$df_pars[elements]

  sum(getBetaVcovDf(model)$df_pars[elements])

  ns

}


#!!!!!Make sure to execute this with L2!!! => a will also be different!!!!
getComp_preds <- function(model, L2){

  comp_preds <- vector("list", ncol(L2))
  for(i in 1:ncol(L2)){
    Lcol <- L2[,i]
    if(!any(is.na(Lcol))){
    #a: levels in L2 matrix that are not 0
    a <- names(Lcol[Lcol!=0])

    #cnms: names of all predictors
    cnms <- names(attr(model,"MSqRob_cnms"))
    #df contains only TRUE/FALSE and has as many rows as a and columns as cnms
    df <- data.frame(lapply(cnms, function(x){return(grepl(paste0("^",x),a))}))
    #If any rowSums(df)!=1 => one or more elements of a start with the name of more than one predictor e.g. if there would be predictors like "sample" and "sample2"
    if(any(rowSums(df)!=1)){stop("Unable to determine which levels come from which predictor. Please use non-overlapping factor names!")} #Is also the case for lme4!
    #All columns with at least one TRUE -> these are predictors that are present in the column names of L2
    comp_preds[[i]] <- cnms[which(colSums(df)!=0)]
    }
  }

  return(comp_preds)
}



getN <- function(protein, options, comp_preds, L2){
  #Start with all independent
  n <- rep(Inf, ncol(L2))
  n_option <- rep(NA, ncol(L2))

  for(i in 1:ncol(L2)){

    #If the contrast is NA => don't try to calculate the partset and all...
    if(!any(is.na(L2[,i]))){

    partset <- getPartSet(protein, Lcol=L2[,i], comp_pred=comp_preds[[i]])
    options2 <- options[options!=comp_preds[[i]]]

    for(k in 1:length(options2)){

      #Er mag geen overlap zijn!
      #Is enkel voldaan als ze ze na paste0 dezelfde lengte hebben!
      if(length(unique(partset[,options2[k]]))==length(unique(apply(partset[,c(comp_preds[[i]],options2[k]), drop=FALSE], 1, paste0, collapse="")))){
        if(length(unique(protein[,options2[k]]))<n[i]){
          n[i] <- length(unique(protein[,options2[k]]))
          n_option[i] <- options2[k]
        }
      }
    }

    #als je contrast op muis is, is je n ook muis
    #Wat als er binnen muis nog een onderverdeling is? -> ook goed, je wilt dan weten hoe de ene muis van de andere verschilt, dus je wilt onafhankelijke herhaalde metingen op elke muis!

    if(is.infinite((n[i]))){n[i] <- length(protein$treat)
    }

  }
  }

  #This means the contrast was NA
  if(is.infinite((n[i]))){n[i] <- NA}

  names(n) <- n_option


  return(n)
}

getPartSet <- function(protein, Lcol, comp_pred){
  a <- names(Lcol[Lcol!=0])

  protein <- .adjustNames(list(protein), comp_pred)[[1]] #Function to be found in file fit_model.R: adjust names of factor variables: paste them with their name

  #All levels of a for all predictors comp_pred should be in the partset
  boollist <- vector("list", length(comp_pred))
  for(j in 1:length(comp_pred)){
    boollist[[j]] <- protein[[comp_pred[j]]] %in% a
  }

  partset <- protein[rowSums(data.frame(boollist))>0,]
  return(partset)
}



getPs <- function(protein, ns, options, comp_preds, L2){

  p_options <- vector("list", ncol(L2))

  for(i in 1:ncol(L2)){

    p_options[[i]] <- c("(Intercept)", comp_preds[[i]])
    partset <- getPartSet(protein, Lcol=L2[,i], comp_pred=comp_preds[[i]])
    n_option <- names(ns[i])
    options2 <- options[options!=n_option]

    for(k in 1:length(options2)){

      #Moeten n omvatten (i.e. mogen geen verschillende levels hebben binnen 1 level van n, zijn anders pseudoreplicaten)
      #Mogen ook niet meer levels hebben dan n
      if(length(unique(protein[,n_option]))==length(unique(apply(protein[,c(n_option,options2[k]), drop=FALSE], 1, paste0, collapse="")))){
        #OUD: mogen geen 2 of meer levels van het contrast omvatten (i.e. moeten volledig tussen contrasten zitten)
        #Betekent opnieuw dat als je ze paste met het contrast, ze nog altijd evenveel levels moeten hebben
        #if(length(unique(partset[,options2[k]]))==length(unique(apply(partset[,c(comp_preds[[i]],options2[k]), drop=FALSE], 1, paste0, collapse="")))){
        #-> klopt niet! -> effecten boven n zoals intercept moeten er ook bij! -> enkel orthogonale mogen niet!

        #Beter: mogen geen 2 of meer levels hebben in EEN level van de n_option (anders orthogonaal)
        #Als er maar 1 Sequence is wordt ze meegenomen! (feitelijk deel van het intercept dan!)
        if(length(unique(protein[,n_option]))==length(unique(paste0(protein[,options2[k]],protein[,n_option])))){
          p_options[[i]] <- c(p_options[[i]],options2[k])
        }
      }
    }

  }

  return(p_options)

}


########################################################################
##########End Between-within degrees of freedom calculation#############
########################################################################
