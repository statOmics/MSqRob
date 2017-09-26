
setGeneric (
  name= "getBetaVcovDf",
  def=function(model, exp_unit=NULL, pars_between=NULL, ...){standardGeneric("getBetaVcovDf")}
)

.getBetaVcovDfLm <- function(model, exp_unit=NULL, pars_between=NULL){
  beta <- as.matrix(model$coefficients,ncol=1)
  vcov <- summary(model)$cov.unscaled
  df <- getDf(model) #model$df.residual
  sigma <- getSigma(model) #We use the MSqRob-defined sigma function instead of #summary(model)$sigma

  #If you declare the experimental units, you can use a more conservative way to estimate the df, which we will give as df_exp
  df_exp <- NULL
  if(!is.null(exp_unit)){

    #If the effects that consume dfs are not given, all effects are used
    #This is not really advisable, mostly within-treatment effects (such as a Sequence effect) are present in the model:
    if(is.null(pars_between)){
      #Number of experimental units minus total number of parameters
      df_exp <- pmax(length(unique(eval(model$call$data)[,exp_unit]))-length(model$assign),0)

      #Else, use the given pars_between
    } else{
      #Number of experimental units minus total number of parameters between treatments "pars_between" (to be declared by the user!)
      #Zal afhangen van L => overweeg verplaatsen!!!
      df_exp <- pmax(length(unique(eval(model$call$data)[,exp_unit]))-sum(model$assign %in% which(names(model$xlevels) %in% pars_between) | model$assign==0),0)
    }
  }

  returnlist <- list(beta=beta,vcov=vcov,df=df,sigma=sigma,df_exp=df_exp)

  return(returnlist)
}

#' Get beta, vcov, df and sigma from an ordinary linear model
#'
#' @description This function returns a list containing the parameter estimates \code{beta}, the variance-covariance matrix \code{vcov},
#' the residual degrees of freedom \code{df} and the residual standard deviation \code{sigma} based on a general linear model fitted by the \code{lm} function from the stats package.
#' This function will only rarely be called by the end-user. When calculating these values for \code{\link[=protLM-class]{protLM}} objects, we recommend using the function \code{getbetaVcovDfList}.
#' @param model A general linear model object of class \code{\link[=lm-class]{lm}}.
#' @param exp_unit The effect in the model that corresponds to the experimental unit. Only needed when one would like to calculate a more conservative way of estimating the degrees of freedom.
#' The default way (\code{exp_unit=NULL}) estimates the degrees of freedom by substracting the total number of observations by the number of parameters. However, often, observations are not completely independent. A more conservative way (\code{df_exp}) is defining on which level the treatments were executed and substracting all degrees of freedom lost due to between-treatement effects (\code{pars_between}) from the number of treatments.
#' @param pars_between Only used if \code{exp_unit} is not \code{NULL}. Character vector indicating all parameters in the model that are between-treatment effects in order to calculate a more conservative degrees of freedom (\code{df_exp}). If left to default (\code{NULL}), all parameters in the model will be asumed to be between-treatment effects (this is not adviced as the result will mostly be too conservative).
#' @return A list containing (1) a named column matrix beta containing the parameter estimates, (2) a named square variance-covariance matrix, (3) a numeric value equal to the residual degrees of freedom, (4) a numeric value equal to the residual standard deviation of the model and (5) \code{NULL} if \code{exp_unit} is left to its default value \code{NULL}, else: a conservative estimate of the degrees of freedom based on the number of experimental units and the degrees of freedom lost due to between-treatment effects.
#' @examples
#' data(proteinsCPTAC, package="MSqRob")
#' lmmodel <- lm(formula="value ~ 1 + conc + instrlab + Sequence",data=getData(proteinsCPTAC[2]))
#' getBetaVcovDf(lmmodel)
#' @include protdata.R
#' @include protLM.R
#' @export
setMethod("getBetaVcovDf", "lm", .getBetaVcovDfLm)

.getBetaVcovDfMermod <- function(model, exp_unit=NULL, pars_between=NULL, Ginvoffset = 1e-18){
  # condVar <- function(model) {
  #   s2 <- sigma(model)^2
  #   Lamt <- getME(model, "Lambdat")
  #   L <- getME(model, "L")
  #
  #   ## never do it this way! fortune("SOOOO")
  #   #V <- solve(L, system = "A")
  #   #V <- chol2inv(L)
  #   #s2*crossprod(Lamt, V) %*% Lamt
  #
  #   LL <- solve(L, Lamt, system = "A")
  #   s2 * crossprod(Lamt, LL)
  # }
  #
  # Kunnen misschien nuttig zijn:
  # attributes(ranef(model, condVar = TRUE)$treat)$postVar
  # attributes(VarCorr(model))$sc^2
  # VarCorr(model)
  # vcov(model)
  #
  # #Deze zijn gelijk!!!!:(als er geen intercept is...)
  # vcov-crossprod(Lamt, LL)

  # #Deze ook:
  # -(sum((resid(model)*gew)^2)/sigma^2-ncol(Ct))
  # sum(Matrix::t(Ct)%*%vcov*Matrix::t(Ct_w)*gew)

  #Design matrix for the fixed effects:
  X <- lme4::getME(model,"X")
  
  if(!is.null(attr(X,"scaled:center")) & !is.null(attr(X,"scaled:scale"))){
  X <- X * attr(X, 'scaled:scale') + attr(X, 'scaled:center')
  }

  #If "Intercept_MSqRob" is present (sometimes if other fixed effects are shrunken in fit.model), change to "(Intercept)"
  colnames(X)[colnames(X)=="Intercept_MSqRob"] <- "(Intercept)"

  #Design matrix for the random effects:
  Z <- lme4::getME(model,"Z")

  R_fix <- attr(model, "MSqRob_R_fix")
  Zt_indices <- attr(model, "MSqRob_Zt_indices")
  unshr_pos <- attr(model, "MSqRob_unshr_pos")

  p <- model@devcomp$dims["p"] #equals ncol(X)
  q <- model@devcomp$dims["q"] #equals ncol(Z)

  #Transpose of Lambda, the relative covariance factor of the random effects:
  lambdat <- lme4::getME(model,"Lambdat")

  #Ginv, the inverse of the estimated variance-covariance matrix of the random effects. We add a small offset Ginvoffset to the diagonal of Ginv to prevent near-singularity:
  Ginv=Matrix::solve(Matrix::crossprod(lambdat)+Matrix::Diagonal(q,Ginvoffset))

  #Transpose of the design matrix C:
  Ct=rbind2(t(X),Matrix::t(Z)) ## t(C), make use of sparse matrices

  #Create gew: the sqrt of the weigths
  if(is.null(model@frame$"(weights)")){gew <- rep(1,ncol(Ct))} else{
    gew <- sqrt(model@frame$"(weights)")}

  #Sqrt-weighted design matrix Ct_w (sparse matrix):
  Ct_w <- Matrix::tcrossprod(Ct,Matrix::Diagonal(length(gew))*gew)

  vcovInv <- Matrix::tcrossprod(Ct_w)
  vcovInv[((p+1):(q+p)),((p+1):(q+p))]=vcovInv[((p+1):(q+p)),((p+1):(q+p))]+Ginv

  if(!is.null(R_fix)){

    R_fix_large <- Matrix::Matrix(diag(nrow(vcovInv)), sparse=TRUE)
    R_fix_large[c(1:p,Zt_indices+p),c(1:p,Zt_indices+p)] <- R_fix
    # R_fix_large[1:p,1:p] <- R_fix[unshr_pos,unshr_pos]
    # R_fix_large[Zt_indices+p,Zt_indices+p] <- R_fix[-unshr_pos,-unshr_pos]
    # R_fix_large[1:p,Zt_indices+p] <- R_fix[unshr_pos,-unshr_pos]#R_fix[-unshr_pos,unshr_pos] #
    # R_fix_large[Zt_indices+p,1:p] <- R_fix[-unshr_pos,unshr_pos]#R_fix[unshr_pos,-unshr_pos] #

    colnames(R_fix_large) <- colnames(vcovInv)
    rownames(R_fix_large) <- rownames(vcovInv)

    #Put the shrunken fixed effects back to their original scale:
    vcovInv <- Matrix::t(R_fix_large)%*%vcovInv%*%R_fix_large
    # Ct[Zt_indices+p,] <- R_fix[-unshr_pos,-unshr_pos]%*%Ct[Zt_indices+p,]
    # Ct_w[Zt_indices+p,] <- R_fix[-unshr_pos,-unshr_pos]%*%Ct_w[Zt_indices+p,]
    # Ct[1:p,] <- R_fix[unshr_pos,unshr_pos]%*%Ct[1:p,]
    # Ct_w[1:p,] <- R_fix[unshr_pos,unshr_pos]%*%Ct_w[1:p,]
    Ct <- Matrix::t(R_fix_large)%*%Ct
    Ct_w <- Matrix::t(R_fix_large)%*%Ct_w
  }

  #if scaled==TRUE:
  #Voor toevoegen ridge penalty:
  #Matrix::tcrossprod(Ct_w*c(1,colSums(Z)))==vcovInv*tcrossprod(c(1,colSums(Z)))
  # if(scaled){vcovInv <- vcovInv*Matrix::tcrossprod(c(1,Matrix::colSums(Z)))
  # Ct_w <- Ct_w*c(1,Matrix::colSums(Z))
  # Ct <- Ct*c(1,Matrix::colSums(Z))}

  #When you only need a part of the output (e.g. for performance reasons): doesn't work!
  # if(!is.null(subset)){
  #   if(is.character(subset)){subset <- which(rownames(Ct) %in% subset)}
  #   vcovInv <- vcovInv[subset,subset,drop=FALSE]
  #   Ct <- Ct[subset,,drop=FALSE]
  #   Ct_w <- Ct_w[subset,,drop=FALSE]
  # }

  #Sometimes, when you work with interaction terms, some factor combinations that are not present in the data might creep into to vcovInv matrix generating rows with only zeros, making it uninvertible
  defined <- rowSums(as.matrix(vcovInv==0))!=ncol(vcovInv)
  defined[is.na(defined)] <- TRUE
  vcovInv <- vcovInv[defined,defined,drop=FALSE]
  Ct <- Ct[defined,,drop=FALSE]
  Ct_w <- Ct_w[defined,,drop=FALSE]

  #Estimated variance-covariance matrix vcov:
  vcov <- tryCatch(as.matrix(Matrix::solve(vcovInv)), error=function(e){
    return(vcovInv*NA)
  })

  #!!!!!Perhaps better alternative where rownames are unnecessary: vcov <- as.matrix(chol2inv(chol(vcovInv)))
  rownames(vcov) <- colnames(vcovInv)
  colnames(vcov) <- rownames(vcovInv)

  #Parameter estimates beta:
  beta <- as.matrix(vcov%*%Ct_w%*%(lme4::getME(model, "y")*gew)) ##getME(model, "y")  equals  model@frame$value  equals  dataset$value

  #Residual standard deviation sigma:
  sigma <- getSigma(model) # We use the MSqRob-defined sigma, could be squeezed sigma in attribute as well!
  ##sigma(model) ## equals: yhat <- t(Ct)%*%beta  and  sigma <- sqrt(sum(((getME(model, "y")-yhat)*gew)^2)/df)

  #Residual degrees of freedom df:
  df <- getDf(model)
  #df <- sum((resid(model)*gew)^2)/sigma^2
  #Equal to: df <- ncol(Ct)-sum(Matrix::t(Ct)%*%vcov*Matrix::t(Ct_w)*gew) ## cheaper version of: Hat <- Matrix::t(Ct)%*%vcov%*%Ct_w*gew  and  df <- ncol(Ct)-sum(Matrix::diag(Hat))

  #Degrees of freedom that you lose due to each parameter
  df_pars <- Matrix::colSums(Matrix::t(Ct)%*%vcov*Matrix::t(Ct_w)*gew)

  #If you declare the experimental units, you can use a more conservative way to estimate the df
  df_exp <- NULL
  if(!is.null(exp_unit)){

    #If the effects that consume dfs are not given, all fixed effects (shrunken and unshrunken) are used
    if(is.null(pars_between)){
      ridgeFix <- which(grepl("ridgeGroup.",names(model@flist)))
      df_indices <- c(1:p,unlist(sapply(ridgeFix, function(x){return((model@Gp[x]+1):model@Gp[x+1])}))+p)

      #Else, use the given pars_between
    } else{
      pars_between2 <- which(attr(model,"MSqRob_fullcnms") %in% pars_between)
      df_indices <- unlist(sapply(pars_between2, function(x){return((attr(model,"MSqRob_XGp")[x]+1):attr(model,"MSqRob_XGp")[x+1])}))
    }
    #!!!!GEEFT PROBLEMEN ALS EXP. UNIT IN DE RIDGEGROUP ZIT!!!! => modelframe aanpassen zodat alles van protein erin zit!!!
    df_exp <- length(unique(model@frame[,exp_unit]))-sum(df_pars[df_indices])  #sum(Matrix::diag(Matrix::t(Ct[df_indices, ,drop=FALSE])%*%vcov[df_indices, ,drop=FALSE]%*%Ct_w*gew))
  }

  #Return beta, vcov, df and sigma as a list:
  returnlist <- list(beta=beta, vcov=vcov, df=df, sigma=sigma, df_exp=df_exp, df_pars=df_pars)

  return(returnlist)
}

#' Get beta, vcov, df and sigma from a linear mixed model
#'
#' @description This function returns a list containing the parameter estimates \code{beta}, the variance-covariance matrix \code{vcov},
#' the residual degrees of freedom \code{df} and the residual standard deviation \code{sigma} based on a mixed effects model fitted by the \code{lmer} function from the lme4 package.
#' This function will only rarely be called by the end-user. When calculating these values for \code{\link[=protLM-class]{protLM}} objects, we recommend using the function \code{\link{getBetaVcovDfList}}.
#' The variance covariance matrix is bias-adjusted, the degrees of freedom are calculated as the trace of the Hat matrix (Ruppert et al., 2003).
#' @param model A linear mixed effects model object of class \code{\link[=lmerMod-class]{lmerMod}}.
#' @param exp_unit The effect in the model that corresponds to the experimental unit. Only needed when one would like to calculate a more conservative way of estimating the degrees of freedom.
#' The default way of estimating the degrees of freedom (\code{exp_unit=NULL}) subtracts the total number of observations by the trace of the Hat matrix. However, often, observations are not completely independent. A more conservative way (\code{df_exp}) is defining on which level the treatments were executed and substracting all degrees of freedom lost due to between-treatement effects (\code{pars_between}) from the number of treatments.
#' @param pars_between Only used if exp_unit is not \code{NULL}. Character vector indicating all parameters in the model that are between-treatment effects in order to calculate a more conservative degrees of freedom (\code{df_exp}). If left to default (\code{NULL}), all parameters in the model will be asumed to be between-treatment effects (this is not adviced as the result will mostly be too conservative).
#' @param Ginvoffset A numeric value indicating the offset added the the diagonal of the G matrix to prevent near-singularity. Defaults to \code{1e-18}.
#' @return A list containing (1) a named column matrix beta containing the parameter estimates, (2) a named square variance-covariance matrix, (3) a numeric value equal to the residual degrees of freedom and (4) a numeric value equal to the residual standard deviation of the model.
#' @examples
#' data(proteinsCPTAC, package="MSqRob")
#' mixedmodel <- lmer(formula="value ~ 1 + (1|conc) + (1|instrlab) + (1|Sequence)",data=getData(proteinsCPTAC[2]))
#' getBetaVcovDf(mixedmodel)
#' @references David Ruppert, M.P. Want and R.J. Carroll.
#' Semiparametric Regression.
#'  Cambridge University Press, 2003.
#' @include protdata.R
#' @include protLM.R
#' @export
setMethod("getBetaVcovDf", "lmerMod", .getBetaVcovDfMermod)
