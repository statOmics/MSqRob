#' Squeeze variances and other parameters
#'
#' @description Squeeze the standard deviations of the models towards a common value using empirical Bayes moderation similar to the \code{\link[=limma]{limma}} package.
#' Information is borrowed over all accessions in the \code{\link[=protLM-class]{protLM}} object in order to stabilize variances. These stabilized variances are called posterior variances. Options are provided to also perform the same squeezing over the random effects and for the penalty of the shrunken fixed effects.
#' @param protLM A \code{\link[=protLM-class]{protLM}} object of which residual standard deviations and/or model parameters need to be squeezed.
#' @param par_squeeze Character vector indicating which model parameters need to be squeezed. When squeezing random effects, provide their names. Fixed effects are present in shrinkage groups, e.g. ridgeGroup.1. If you want them to be squeezed as well, provide the names of the shrinkage groups that need to be squeezed. The default \code{NULL} indicates that no parameters will be squeezed.
#' @param squeezeVar A logical indicating whether the residual standard deviation of all models should be squeezed towards a common value. Defaults to \code{TRUE}. If set to \code{FALSE}, residual standard deviations will not be squeezed.
#' @param min_df A numeric value indicating the minimal degrees of freedom that will be taken into account for calculating the prior degrees of freedom and prior variance.
#' @param robust_var A logical indicating wheter the estimation of the prior degrees of freedom and the prior variance (needed for shrinkage) should be robustified against outlier variances. Defaults to \code{TRUE}.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message_thetas Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the extraction of the variances, or \code{NULL} to hide the current message (if any).
#' @param message_squeeze Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the squeezing of the variances, or \code{NULL} to hide the current message (if any).
#' @param message_update Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the updating of the models, or \code{NULL} to hide the current message (if any).
#' @param ... Other arguments to be passed to the \code{update_protLM} function internally.
#' @return An updated protLM object with squeezed variances and/or squeezed parameters (or the input \code{\link[=protLM-class]{protLM}} object if \code{par_squeeze=NULL} and \code{squeezeVar=FALSE}).
#' @examples ....
#' @references ....
#' @include protLM.R
#' @include squeezeVarRob.R
#' @export
squeezePars <- function(protLM, par_squeeze=NULL, squeezeVar=TRUE, min_df=1, robust_var=TRUE, printProgress=FALSE, shiny=FALSE, message_thetas=NULL, message_squeeze=NULL, message_update=NULL, ...){

  if(!is.null(par_squeeze) || isTRUE(squeezeVar)){
  #Get variances on the parameters you want to squeeze and their df
  thetaVars <- getThetaVars(protLM, par_names=par_squeeze, printProgress=printProgress, shiny=shiny, message=message_thetas)

  thetas <- thetaVars$thetas
  df_thetas <- thetaVars$df_thetas
  vars <- thetaVars$vars
  df_vars <- thetaVars$df_vars

  #Squeeze variances and extract new thetas
  thetasVarsDf_post <- squeezeThetas(thetas, df_thetas, vars, df_vars, squeezeVar=squeezeVar, min_df=min_df, robust_var=robust_var, printProgress=printProgress, shiny=shiny, message=message_squeeze)
  thetas_new <- thetasVarsDf_post$thetas_post

  #update protLM with new thetas, new sigmas and new df_sigmas

  #!!Check volgorde van squeezen!!!
  protLM <- update_protLM(protLM,thetas_new,sigmas=sqrt(thetasVarsDf_post$vars_post),df_sigmas=thetasVarsDf_post$df_vars_post, robust_var=robust_var, printProgress=printProgress, shiny=shiny, message=message_update, ...)

  }

  return(protLM)
}


#' Get variances and degrees of freedom of model parameters
#'
#' @description This function extracts the estimated variances of specified random effects and shrunken fixed effects as well as their associated degrees of freedom (based on the trace of the Hat matrix). It also always returns the residual variance and the residual degrees of freedom.
#' @param protLM A \code{\link[=protLM-class]{protLM}} object of which residual variances and/or model parameters need to be returned.
#' @param par_names Character vector indicating of which model parameters the variance needs to be returned. When squeezing random effects, provide their names. Fixed effects are present in shrinkage groups, e.g. ridgeGroup.1. If you want to return their variance as well, provide the names of the shrinkage groups that need to be squeezed. If \code{par_names=NULL}, \code{NA} will be returned in the \code{thetas} and \code{df_thetas} slots of the output.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the extracting of the variances and degrees of freedom, or \code{NULL} to hide the current message (if any).
#' @return A named list with 4 slots. The first slot \code{thetas} contains a matrix with in each column the estimated variances for an effect specified in the \code{par_names} argument, each row corresponds to a different accession in the \code{\link[=protLM-class]{protLM}} object.
#' The second slot \code{df_thetas} contains a matrix of similar structure to \code{thetas} but containing the degrees of freedom corresponding to the estimated variances. The third slot \code{vars} contains a vector of residual variances for each accession and the fourth slot \code{df_vars} contains a vector of residual degrees of freedom.
#' @examples ....
#' @references ....
#' @export
getThetaVars <- function(protLM, par_names, printProgress=FALSE, shiny=FALSE, message=NULL){

  #Progress bar for extracting thetas
  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  df_thetas <- matrix(NA, nrow=length(protLM), ncol=length(par_names), dimnames=list(NULL,par_names))
  thetas <- matrix(NA, nrow=length(protLM), ncol=length(par_names), dimnames=list(NULL,par_names))
  df_vars <- rep(NA,length(protLM))
  vars <- rep(NA,length(protLM))

for(i in 1:length(protLM)){

  updateProgress(progress=progress, detail=paste0("Extracting variance for model ",i," of ",length(protLM),"."), n=length(protLM), shiny=shiny, print=isTRUE(printProgress))

lmerMod <- getModels(protLM[i])

if(is.null(lmerMod@frame$"(weights)")){gew <- rep(1,nobs(lmerMod))} else{
  gew <- sqrt(lmerMod@frame$"(weights)")}

theta_names <- names(lmerMod@flist)

#thetas
thetas[i,] <- vapply(par_names,function(x) {
  y <- lmerMod@pp$theta[which(theta_names==x)]
  if(length(y)==0) y <- NA
  return(y)
},0)

#df thetas
df_thetas[i,] <- vapply(par_names,function(x) {
  y <- sum(lme4::ranef(lmerMod)[[x]]^2)/unlist(lme4::VarCorr(lmerMod)[x])

  #If not specified => NA
  if(length(y)==0){ y <- NA
  #If NaN => 0 df
  } else if(is.na(y)){ y <- 0}

  return(y)
  },0)

#residual variances
vars[i] <- sigma(lmerMod)^2

#df residual variances
df_vars[i] <- sum((resid(lmerMod)*gew)^2)/vars[i]
}
  #df_thetas[is.na(df_thetas)] <- NA
  thetavars <- list(thetas=thetas, df_thetas=df_thetas, vars=vars, df_vars=df_vars)
  return(thetavars)
}


#' Squeeze variances and other parameters
#'
#' @description Squeeze a set of sample variances towards a common value using empirical Bayes moderation similar to the \code{\link[=limma]{limma}} package.
#' This can be done for residual variances, random effects or shrunken fixed effects. The resulting variances are called posterior variances.
#' @param thetas A numeric matrix wherein each column contains variances for a different parameter.
#' @param df_thetas A numeric matrix wherein each column contains degrees of freedom corresponding to the \code{thetas}.
#' @param vars A numeric vector containing residual variances.
#' @param df_vars A numeric vector containing residual degrees of freedom.
#' @param squeezeVar A logical indicating whether the residual variances as given by the argument \code{vars} should be squeezed towards a common value. Defaults to \code{TRUE}. If set to \code{FALSE}, \code{vars} will not be squeezed and the \code{vars_post} slot will contain the same values as the \code{vars} input argument.
#' @param min_df A numeric value indicating the minimal degrees of freedom that will be taken into account for calculating the prior degrees of freedom and prior variance.
#' @param robust_var A logical indicating wheter the estimation of the prior degrees of freedom and the prior variance (needed for shrinkage) should be robustified against outlier variances. Defaults to \code{TRUE}.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the squeezing of the variances and degrees of freedom, or \code{NULL} to hide the current message (if any).
#' @param ... Other arguments to be passed to the \code{squeezeVarRob} function internally.
#' @return A named list with 4 slots. The first slot \code{thetas_post} contains a matrix with in each column the squeezed variances corresponding to the \code{thetas} input argument.
#' The second slot \code{df_thetas_post} contains a matrix of similar structure to \code{thetas_post} but containing the posterior degrees of freedom corresponding to the \code{df_thetas} input argument. If \code{squeezeVar=TRUE}, the third slot \code{vars_post} contains a vector of posterior residual variances for each accession. If \code{squeezeVar=FALSE}, it contains the residual variances given in the \code{vars} input argument.
#' If \code{squeezeVar=TRUE}, the fourth slot \code{df_vars_post} contains a vector of posterior residual degrees of freedom. If \code{squeezeVar=FALSE}, it contains the residual degrees of freedom given in the \code{df_vars} input argument.
#' @examples ....
#' @references ....
#' @export
squeezeThetas <- function(thetas, df_thetas, vars, df_vars, squeezeVar=TRUE, min_df=1, robust_var=TRUE, printProgress=FALSE, shiny=FALSE, message=NULL, ...){

  #Progress bar for extracting thetas
  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  df_thetas_post <- matrix(NA, nrow=nrow(df_thetas), ncol=ncol(df_thetas), dimnames=dimnames(df_thetas))
  thetas_post <- matrix(NA, nrow=nrow(thetas), ncol=ncol(thetas), dimnames=dimnames(thetas))

  #Remove NaN
  df_thetas[is.na(df_thetas)] <- NA

  if(isTRUE(squeezeVar)){
    sqVarObj <- squeezeVarRob(vars, df=df_vars, min_df=min_df, robust=robust_var, ...)
    df_vars_post <- df_vars+sqVarObj$df.prior
    vars_post <- sqVarObj$var.post
  } else{
    vars_post <- vars
    df_vars_post <- df_vars
  }

  if(ncol(thetas)!=0){
    for(j in 1:ncol(thetas)){

      updateProgress(progress=progress, detail=paste0("Squeezing variances for model ",j," of ",ncol(thetas),"."), n=ncol(thetas), shiny=shiny, print=isTRUE(printProgress))

      sqVarObj <- tryCatch(
      squeezeVarRob(thetas[,j]*vars, df=df_thetas[,j], min_df=min_df, robust=robust_var, ...)
      , error=function(e){
        sqVarObj <- data.frame(var.post=thetas[,j]*vars, df.prior=0)

        if(!is.null(colnames(thetas)[j])){
          failedTheta <- colnames(thetas)[j]
        } else{paste0(failedTheta <- "theta[",j,"]")}

        warning(paste0("Could not squeeze thetas for \"",failedTheta,"\".")) #Maybe make error more informative.

        return(sqVarObj)
        })
      df_thetas_post[,j] <- df_thetas[,j]+sqVarObj$df.prior
      thetas_post[,j] <- sqVarObj$var.post/vars
      }
  } else{df_thetas_post <- df_thetas
    df_thetas_post <- thetas_post}

  thetasVarsDf_post <- list(thetas_post=thetas_post, df_thetas_post=df_thetas_post, vars_post=vars_post, df_vars_post=df_vars_post)
  return(thetasVarsDf_post)
}

#lol2 <- update.lmerMod(lol,theta)


#' Update lmerMod object with new thetas without changing the pointers of other objects
#'
#' @description This function changes the \code{theta} parameter of a \code{\link[=lmerMod-class]{lmerMod}} object to a given vector and updates the model accordingly.
#' @param object A \code{\link[=lmerMod-class]{lmerMod}} object of which residual standard deviations and/or model parameters need to be squeezed.
#' @param theta A numeric vector of theta values (parameter variances) that will replace the theta values in the \code{\link[=lmerMod-class]{lmerMod}} object.
#' @param ... Other arguments to be passed to the \code{\link{update.default}} function if the \code{theta} input argument would be missing.
#' @return A \code{\link[=lmerMod-class]{lmerMod}} object of which the parameter variances are replaced by the ones given in the \code{theta} input argument.
#' @examples ....
#' @references Bolker (2016). Wald errors of variances, https://rpubs.com/bbolker/waldvar
#' @export
update.lmerMod <- function(object, theta, ...) {
  if (missing(theta)) return(update.default(object, ...))

  #If all thetas are NA => don't waste any time -> leads only to errors (e.g. in dd)
  if(all(is.na(theta))){
    output_object <- object

  #If normal thetas => update!
  } else{
  object2 <- object
  ## deep-copy the (only) reference-class slots ...
  object2@pp <- object@pp$copy()
  object2@resp <- object@resp$copy()
  object2@pp$setTheta(theta)
  dd <- as.function(object2)
  dval <- dd(theta)  ## update internal structures
  ## object2@resp$updateMu(object2@resp$mu)  ## ?? not helping/not necessary

  output_object <- mkMerMod(environment(dd),
           opt=list(par=theta,fval=dval,conv=0),
           fr=model.frame(object2),
           reTrms=getME(object2,
                        c("Zt","Lambdat","Lind","theta",
                          "lower","flist","cnms","Gp"))
  )
  }

  attr(output_object,"MSqRob_R_fix") <- attr(object,"MSqRob_R_fix")
  attr(output_object,"MSqRob_Q_fix") <- attr(object,"MSqRob_Q_fix")
  attr(output_object,"MSqRob_Zt_indices") <- attr(object,"MSqRob_Zt_indices")
  attr(output_object,"MSqRob_unshr_pos") <- attr(object,"MSqRob_unshr_pos")
  attr(output_object,"MSqRob_levels") <- attr(object,"MSqRob_levels")
  attr(output_object,"MSqRob_cnms") <- attr(object,"MSqRob_cnms")
  attr(output_object,"MSqRob_Gp") <- attr(object,"MSqRob_Gp")
  attr(output_object,"MSqRob_XGp") <- attr(object,"MSqRob_XGp")
  attr(output_object,"MSqRob_Xp") <- attr(object,"MSqRob_Xp")
  attr(output_object,"MSqRob_fullcnms") <- attr(object,"MSqRob_fullcnms")
  #Only used afterwards in a typical workflow, still added just in case
  attr(output_object,"MSqRob_sigma") <- attr(object,"MSqRob_sigma")
  attr(output_object,"MSqRob_df_sigma") <- attr(object,"MSqRob_df_sigma")

  return(output_object)
}

#' Update protLM object with new parameter variances and residual variances
#'
#' @description This function changes the parameter variances \code{theta} of all \code{\link[=lmerMod-class]{lmerMod}} objects present in the \code{model} slot of a \code{\link[=protLM-class]{protLM}} object to a given vector and updates these models accordingly.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the  all of its \code{\link[=lmerMod-class]{lmerMod}} object in the \code{model} slot should have one or more parameter parameter variances \code{theta} replaced by given theta values. Optionally residual variances can also be updated.
#' @param theta A numeric matrix wherein each column contains parameter variances (theta values) that will replace the theta values in all \code{\link[=lmerMod-class]{lmerMod}} objects present in the \code{model} slot of the \code{\link[=protLM-class]{protLM}} object. Column names should be corresponding to random effects and/or ridge groups (e.g. "ridgeGroup.1") present in the model.
#' @param sigmas A vector of length equal to the \code{\link[=protLM-class]{protLM}} object containing residual variances that should replace the existing residual variances in each object present in the \code{model} slot of the \code{\link[=protLM-class]{protLM}} object.
#' @param df_sigmas A vector of length equal to the \code{\link[=protLM-class]{protLM}} object containing residual degrees of freedom that should replace the existing residual degrees of freedom in each object present in the \code{model} slot of the \code{\link[=protLM-class]{protLM}} object.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the updating of the models, or \code{NULL} to hide the current message (if any).
#' @param ... Other arguments to be passed to the \code{\link{update.lmerMod}} function.
#' @return A \code{\link[=lmerMod-class]{lmerMod}} object of which the parameter variances are replaced by the ones given in the \code{theta} input argument.
#' @examples ....
#' @references Bolker (2016). Wald errors of variances, https://rpubs.com/bbolker/waldvar
#' @export
update_protLM <- function(protLM, theta, sigmas=NULL, df_sigmas=NULL, printProgress=FALSE, shiny=FALSE, message=NULL, ...) {

  #Progress bar for extracting thetas
  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  modellist <- getModels(protLM, simplify=FALSE)
  theta_names <- colnames(theta)

  for(i in 1:length(protLM)){

    updateProgress(progress=progress, detail=paste0("Updating model ",i," of ",length(protLM),"."), n=length(protLM), shiny=shiny, print=isTRUE(printProgress))

    x <- modellist[[i]]
    if(class(x)=="lmerMod"){
    theta_updated <- x@theta
    index <- match(theta_names,names(x@flist))
    if(any(is.na(index))){
      theta_updated <- theta_updated[-which(is.na(index))]
      index <- na.omit(index)
    }
    theta_updated[index] <- theta[i,]

    #If Downdated VtV is not positive definite error => then insert zeros
    #Check for other datasets how this works!!!
    modellist[[i]] <-
      tryCatch(update.lmerMod(x,theta_updated,...), error=function(e){
        return(update.lmerMod(x,(theta_updated*0),...))
        })
    }
    #Update sigma and df_sigma, pass them as attributes
    attr(modellist[[i]],"MSqRob_sigma") <- sigmas[i]
    attr(modellist[[i]],"MSqRob_df_sigma") <- df_sigmas[i]

  }
return(new("protLM", accession=protLM@accession, model=modellist, annotation=protLM@annotation))
}

