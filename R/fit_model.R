#' Fit peptide-based models
#'
#' @description Fits a model to each protein of a \code{\link[=protdata-class]{protdata}} object and returns a corresponding \code{\link[=protLM-class]{protLM}} object.
#' In its standard settings, the function returns a \code{\link[=protLM-class]{protLM}} object containing robust ridge models as described in Goeminne et al. (2015).
#' However, the user can also specify to turn off the ridge regression and fit the models by ordinary least squares (OLS) and/or to turn off the down-weighing of outliers by M estimation with Huber weights and/or to turn off the Empirical Bayes squeezing of variances.
#' @param protdata A \code{\link[=protdata-class]{protdata}} object to which peptide-based models must be fitted. Note that factors should be coded as factors and numeric variables as numeric in each data slot.
#' @param response The name of the column in the data slot of the \code{\link[=protdata-class]{protdata}} object that contains the response variable for the model, mostly the column containing the log2 transformed intensity values.
#' @param fixed Either a vector of names corresponding to the columns in the data slot of the \code{\link[=protdata-class]{protdata}} object containing the predictor variables, or a right-hand sided fixed effects formula without intercept, which should be indicated in the argument \code{add.intercept}. \code{NULL} (the default) indicates that no fixed effects other than a possible fixed intercept are present in the model.
#' @param random Either a vector of names corresponding to the columns in the data slot of the \code{\link[=protdata-class]{protdata}} object containing the predictor variables or a right-hand sided random effects formula. Adding the peptide sequences as one of the random effects predictors is highly recommended as individual peptide effects are often quite strong. \code{NULL} (the default) indicates that no random effects are present in the model.
#' @param add.intercept A logical value indicating whether the fitted models should contain a fixed intercept. If missing, the value is set to \code{TRUE}, indicating the intercept should be present in the model.
#' @param shrinkage.fixed A numeric vector containing only 0 and/or 1 of length equal to the number of fixed effects, potential intercept included. The nth element in the shrinkage.fixed vector indicates whether the nth fixed effect should be shrunken (value 1) or not (value 0). If \code{add.intercept=TRUE}, the first element of the vector indicates the intercept. \code{shrinkage.intercept = NULL} (default) indicates all fixed effects except the intercept should be shrunken.
#' @param weights The type of weighing that should be performed. Supported weighing methods incluce \code{"Huber"} (the default) for M estimation with Huber weights and \code{NULL} when no weighing should be applied. One can also supply a list of weights with length equal to the number of proteins in the \code{\link[=protdata-class]{protdata}} object. Each element of the list should either contain \code{"Huber"} or \code{NULL} or a numeric vector containing weights with length equal to the number of observations for that protein.
#' @param k The tuning constant for the Huber mean when weighing down outliers. The default (\code{k = 1.345}) produces 95 \% efficiency relative to the sample mean when the population is normal but provides substantial resistance to outliers when it is not.
#' @param par_squeeze A character vector indicating which model parameters need to be squeezed. When squeezing random effects, provide their names. Fixed effects are present in shrinkage groups, e.g. ridgeGroup.1. If you want them to be squeezed as well, provide the names of the shrinkage groups that need to be squeezed. The default \code{NULL} indicates that no parameters will be squeezed.
#' @param squeezeVar A logical indicating whether the residual standard deviation of all models should be squeezed towards a common value. Defaults to \code{TRUE}. If set to \code{FALSE}, residual standard deviations will not be squeezed.
#' @param min_df A numeric value indicating the minimal degrees of freedom that will be taken into account for calculating the prior degrees of freedom and prior variance. Only used when \code{par_squeeze=TRUE} or \code{squeezeVar} is not \code{NULL}.
#' @param robust_var A logical indicating wheter the estimation of the prior degrees of freedom and the prior variance (needed for shrinkage) should be robustified against outlier variances. Only used when \code{par_squeeze=TRUE} or \code{squeezeVar} is not \code{NULL}. Defaults to \code{TRUE}.
#' @param robustM A logical indicating wheter the weighing of the observations in the M estimation step should be based on the Huber weights of the residuals devided by the median absolute deviation (mad) of the residuals instead of the residuals devided by the residual standard deviation. Older MSqRob implementations used the approach with the residual standard deviation as a default. With the new approach, downweighing outliers will be more efficient, as the mad is not strongly influenced by outliers. Defaults to \code{TRUE}.
#' @param modfiedGS A logical indicating wheter the Modified Gram-Schmidt algorithm should be used for orthogonalizing the fixed effects design matrix instead of the regular QR decomposition algorithm. Defaults to \code{FALSE}.
#' @param tolPwrss A numeric value indicating the maximally tolerated deviation on the penalized weighted residual sum of squares when iteratively estimating the weights by M estimation.
#' @param verbose A logical value indicating whether the models should be printed out. Defaults to \code{FALSE}.
#' @param printProgress A logical indicating whether the R should print a message before fitting each model. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message_fitting Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the fitting of the models, or \code{NULL} to hide the current message (if any).
#' @param message_thetas Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the extraction of the variances, or \code{NULL} to hide the current message (if any).
#' @param message_squeeze Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the squeezing of the variances, or \code{NULL} to hide the current message (if any).
#' @param message_update Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user during the updating of the models, or \code{NULL} to hide the current message (if any).
#' @param ... Other parameters to be passed to the model fitting function internally.
#' @return A \code{\link[=protLM-class]{protLM}} object containing the names of all proteins in the input \code{\link[=protdata-class]{protdata}} object and their corresponding fitted models.
#' @examples data(proteinsCPTAC, package="MSqRob")
#' #Fitting models for the first 10 proteins in the protdata object proteinsCPTAC using the robust ridge approach of Goeminne et al. (2015):
#' modelRRCPTAC <- fit.model(protdata=proteinsCPTAC[1:10], response="value", fixed="conc", random=c("Sequence","instrlab"), add.intercept=TRUE, shrinkage.fixed=NULL, weights="Huber", k = 1.345, tolPwrss = 1e-10, verbose=FALSE)
#' #Fitting models for the first 10 proteins in the protdata object proteinsCPTAC using an ordinary least squares approach (i.e. no shrinkage and no M estimation):
#' modelLmCPTAC <- fit.model(protdata=proteinsCPTAC[1:10], response="value", fixed=c("conc","Sequence","instrlab"), random=NULL, add.intercept=TRUE, shrinkage.fixed=c(0,0,0), weights=NULL, k = 1.345, tolPwrss = 1e-10, verbose=FALSE)
#' @references Ludger Goeminne, Kris Gevaert and Lieven Clement.
#' Peptide-level robust ridge regression improves estimation, sensitivity and specificity in data-dependent quantitative label-free shotgun proteomics.
#'  Molecular & Cellular Proteomics, 2016.
#' @include protdata.R
#' @include protLM.R
#' @include dummyVars_MSqRob.R
#' @include updateProgress.R
#' @export
fit.model=function(protdata, response=NULL, fixed=NULL, random=NULL, add.intercept=TRUE, shrinkage.fixed=NULL, weights="Huber", k = 1.345, par_squeeze=NULL, squeezeVar=TRUE, min_df=1, robust_var=TRUE, robustM = TRUE, scaleUnshrFix = FALSE, modfiedGS = FALSE, tolPwrss = 1e-10, verbose=FALSE, printProgress=FALSE, shiny=FALSE, message_fitting = NULL, message_thetas = NULL, message_squeeze = NULL, message_update = NULL, ...)
{
  datalist <- getData(protdata, simplify=FALSE)

  fixed_input <- .makeFormulaPredictors(fixed, add.intercept, effect="fixed")
  random_input <- .makeFormulaPredictors(random, add.intercept, effect="random")

  fixed <- fixed_input[[1]]
  random <- random_input[[1]]

  formula_fix <- fixed_input[[2]]
  formula_ran <- random_input[[2]]

  #If weights is a single element, turn it into a list
  if(!is.list(weights)){weights <- rep(list(weights), length(datalist))}

  #Do checks for appropriate input values
  .errorCheckFit(possible_vars=colnames(datalist[[1]]), response, fixed, random, weights, n_prot=length(datalist))

  ### Initialize variables ###

  #Variables for when the model cannot be fit:
  intercept <- TRUE

  #######

  modellist <- vector("list", length(datalist))

  #1. We need a mixed model (lme4 - lmer)
  if(!is.null(random) || !(all(shrinkage.fixed==0) && is.numeric(shrinkage.fixed))){

    #Adjust names of predictors to make them similar to those of "lm"
    #In case of factors, the name of the factor is concatenated with the level of the factor
    #Also takes a bit of time with big datasets:
    datalist <- adjustNames(datalist, random)

    #Fit the list of ridge models
    modellist <- .createRidgeList(datalist = datalist, weights = weights, response = response, fixed = fixed, shrinkage.fixed = shrinkage.fixed, formula_fix = formula_fix, random = random, formula_ran = formula_ran, add.intercept = add.intercept, intercept = intercept, intercept_name = "(Intercept)", k = k, robustM = robustM, scaleUnshrFix = scaleUnshrFix, modfiedGS = modfiedGS, tolPwrss = tolPwrss, verbose = verbose, printProgress=printProgress, shiny=shiny, message_fitting=message_fitting, ...)

    #2. We need a simple linear regression model (lm)
  } else if(is.null(random) && (all(shrinkage.fixed==0) && is.numeric(shrinkage.fixed))){

    #Initialize variables specifically for OLS models
    formula <- formula(paste0(response,paste0(as.character(formula_fix), collapse="")))

    #Fit the list of lm models
    modellist <- .createLmList(datalist = datalist, weights = weights, formula = formula, predictors = fixed, response = response, add.intercept = add.intercept, intercept_name = "(Intercept)", k = k, robustM = robustM, tolPwrss = tolPwrss, verbose = verbose, printProgress=printProgress, shiny=shiny, message_fitting=message_fitting, ...)

  }

  protLM <- new("protLM", accession=protdata@accession, model=modellist, annotation=protdata@annotation)
  protLM <- squeezePars(protLM, par_squeeze=par_squeeze, squeezeVar=squeezeVar, min_df=min_df, robust_var=robust_var, printProgress=printProgress, shiny=shiny, message_thetas=message_thetas, message_squeeze=message_squeeze, message_update=message_update)

  return(protLM)

}

.errorCheckFit <- function(possible_vars, response, fixed, random, weights, n_prot){

  #Error control: fixed and random effects should be present as colnames in the data slot ("possible_vars")

  errorMsg <- NULL

  #Control: fixed and random connot be NULL at the same time!
  if(is.null(response)){errorMsg <- paste0(errorMsg, "\n\n", "Please specify a response variable.")}
  if(is.null(fixed) && is.null(random)){errorMsg <- paste0(errorMsg, "\n\n", "Please specify appropriate fixed and/or random effects.")}

  #Control: fixed and random effects must be completely different: no overlaps, otherwise problems with "adjustNames" function!
  #If you really want to try something this crazy, just duplicate the effect.
  if(any(fixed %in% random)){errorMsg <- paste0(errorMsg, "\n\n", "Fixed and random effects must be different from each other.")}

  #Correct for possible interactions (only used in error control that follows)
  errorMsg <- .checkPredictors(errorMsg, fixed, possible_vars, "fixed")
  errorMsg <- .checkPredictors(errorMsg, random, possible_vars, "random")

  if(length(response)!=1){errorMsg <- paste0(errorMsg, "\n\n", "Please provide exactly one response variable.")}
  if(!(response %in% possible_vars)){
    errorMsg <- paste0(errorMsg, "\n\n", paste0("The response variable should be one of: \"",paste0(possible_vars, collapse="\", \""),"\"."))
  }

  #Check if length of weights equals length of protdata object
  if(length(weights)!=n_prot){errorMsg <- paste0(errorMsg, "\n\n", "The length of list \"weights\" should be equal to the number of proteins in the dataset.")}

  if(!is.null(errorMsg)){stop(errorMsg)}
}

.checkPredictors <- function(errorMsg, predictors, possible_vars, type_predictor){

  if(!is.null(predictors)){
    predictor_all <- unique(unlist(strsplit(predictors,":")))
    if(!all(predictor_all %in% possible_vars)){
      not_present <- which(!(predictor_all %in% possible_vars))
      errorMsg <- paste0(errorMsg, "\n\n",paste0("The following ",type_predictor," effects are no possible predictors: \"",paste0(predictor_all[not_present], collapse="\", \""),"\". Please choose from the following predictors: \"",paste0(possible_vars, collapse="\", \""),"\"."))
    }
  }
  return(errorMsg)
}

.removeSingleLevels <- function(fixed, shrinkage.fixed, formula_fix, x, add.intercept){

  morelvls <- .getMoreLevels(fixed, x)
  fixedOld <- fixed
  fixed <- fixedOld[morelvls]

  if(isTRUE(add.intercept)){
    #If add.intercept==TRUE, add an intercept
    terms <- c(1,fixed)
    #Keep the shrinkage term for the intercept, add the other terms
    shrinkage.fixed <- c(shrinkage.fixed[1],(shrinkage.fixed[-1])[morelvls])
    #If add.intercept==FALSE, no intercept
  }    else{
    terms <- c(-1,fixed)
    shrinkage.fixed <- shrinkage.fixed[morelvls]
  }

  #formula_fix <- paste0(c("~",rep("",length(terms)-1)),terms, collapse="+")
  #Remove the ones with only one level from the formula. Advantage of this approach: also works for other formulae!
  if(length(fixed)!=length(fixedOld)){
  formula_fix <- do.call("update", list(as.formula(formula_fix),paste0("~ . - ",paste0(fixedOld[!morelvls], collapse = "-"))))
  }

  return(list(fixed, shrinkage.fixed, formula_fix)) #as.formula(formula_fix)
}

.getMoreLevels <- function(predictors, x){
  #return(unlist(lapply(predictors, function(y) length(levels(x[,y]))!=1)))
  if(is.null(predictors)){return(FALSE)
  }else{
    return(unlist(lapply(strsplit(predictors,":"), function(y) prod(unlist(lapply(x[,y,drop=FALSE], function(z) prod(length(unique(z))))))!=1)))
  }
}


#A function to turn predictor variables into a formula
.makeFormulaPredictors <- function(input, add.intercept, effect){

  #Check the input type: TRUE if formula, FALSE if vector, NA if NULL
  type <- grepl("~",input)[1]

  #If it's a formula (type==TRUE), keep the formula
  if(isTRUE(type)){

    formula <- input

    #Predictors: deconstruct the formula
    #For fixed effects, we also want to catch interaction effects with ":" and account for possible "*"
    #For random effects, we want the terms as if they were given manually
    if(effect=="fixed"){
      predictors <- attr(terms(as.formula(formula)),"term.labels")
    } else{
      predictors <- all.vars(as.formula(formula))
    }

    #If it's NULL (type==NA), only an intercept if there is one, else keep NULL
  }  else if(is.na(type)){

    #If add.intercept==TRUE, add an intercept to the fixed effects, don't do it for random effects
    if(isTRUE(add.intercept) && effect=="fixed"){
      formula <- "~1"
      #We return our predictors WITHOUT a possible intercept term
      predictors <- NULL
      #If add.intercept==FALSE, no intercept
    }    else{
      formula <- NULL
      predictors <- NULL}

    #If it's a vector (type==FALSE), turn it into a formula
  }  else {

    #We return our predictors WITHOUT a possible intercept term
    predictors <- input

    #Make the formula

    #1. If it's fixed
    if(effect=="fixed"){

      #If add.intercept==TRUE, add an intercept
      if(isTRUE(add.intercept)){input <- c(1,input)
      #If add.intercept==FALSE, no intercept
      }    else{input <- c(-1,input)}

      formula <- paste0(c("~",rep("",length(input)-1)),input, collapse="+")

      #2. If it's random
    }  else if(effect=="random"){formula <- paste0(c("~(1|",rep("(1|",length(input)-1)),input, ")", collapse="+")}

  }
  return(list(predictors,as.formula(formula)))
}

###Replace fixed terms by (1|ridgeGroups) in a formula
.replaceFixed <- function(formula, fixed){

  form_new <- gsub("\ ","",formula)
  sub_fix <- paste0("\\+\\(1\\|",fixed,"\\)")

  #Replace the first "fixed" term by (1|ridgeGroups)
  form_new <- gsub(sub_fix[1],"+(1|ridgeGroups)",form_new)

  #If there are other fixed terms, remove them from the formula, they will be captured in ridgeGroups
  if(length(sub_fix)>1){

    for(i in 1:length(sub_fix)){
      form_new <- gsub(sub_fix[i],"",form_new)
    }

  }

  return(form_new)
}

#Create a list with fitted lm regression models
.createLmList=function(datalist, weights, formula, predictors, response, add.intercept, intercept_name = "(Intercept)", k, robustM, tolPwrss, verbose, printProgress=NULL, shiny=FALSE, message_fitting=NULL, ...){

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message_fitting, value = 0)
  }

  count <- 0

  #Return a fitted lmerMod model with M estimation or fitted with given weights or NULL weights
  mapply(function(x,y){

    count <<- count+1
    updateProgress(progress=progress, detail=paste0("Fitting model ",count," of ",length(datalist),"."), n=length(datalist), shiny=shiny, print=isTRUE(printProgress))

    n <- nrow(x)
    #If the weighs for this particular protein are of length 1, duplicate them to the correct length
    if(length(y)==1){y <- rep(y,n)}

    #Remove fixed effects with only one level so that these models can be fit (interpretation is done when fitting contrasts!)
    fixedTmp <- .makeFormulaPredictors(formula, add.intercept, effect="fixed")[[1]]
    templist <- .removeSingleLevels(fixed=fixedTmp, shrinkage.fixed=rep(0, length(fixedTmp)), formula, x, add.intercept)
    formulaAdj <- templist[[3]]

    return(
      tryCatch(
        .lmEstM3(formula=formulaAdj,data=x,k=k,weights=y,robustM=robustM,tolWrss=tolPwrss,verbose=verbose,...)
        , error=function(e){
          return(.emptyLM(formulaAdj,x,y,predictors,response,add.intercept,intercept_name,n))})
    )

  }, datalist, weights, SIMPLIFY = FALSE)

}

#Function for OLS regression
.lmEstM3=function(formula,data,k,weights,robustM=TRUE,tolWrss = 1e-10,verbose=FALSE,...)
{

  #If OLS regression with Huber weights:
  if(!is.null(weights) && weights[1]=="Huber"){lmFit <- .huberLmFit(formula,data,k,robustM,tolWrss,verbose,...)}

  #Additional weighing options could be added here

  #Else: OLS regression with given weights (or NULL):
  else{lmFit <- .weightsLmFit(formula,data,weights,verbose,...)}

  return(lmFit)
}


.weightsLmFit=function(formula,data,weights,verbose,...)
{
  #environment(formula) <- environment(weights)
  .weight. <- weights
  e <- new.env(parent=environment(formula))
  assign(".weight.", weights, envir=e)
  environment(formula) <- e

  lmFit <- lm(formula=formula,data=data,weights=.weight.,...)
  if (verbose) print(lmFit)

  return(lmFit)
}


#Create a list with fitted ridge regression models
.createRidgeList=function(datalist, weights, response, fixed, shrinkage.fixed, formula_fix, random, formula_ran, add.intercept, intercept, intercept_name = "(Intercept)", k, robustM, scaleUnshrFix, modfiedGS,  tolPwrss, verbose, printProgress=NULL, shiny=FALSE, message_fitting=NULL, ...){

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message_fitting, value = 0)
  }

  count <- 0

  modellist <- mapply(function(x,y){

    count <<- count+1
    updateProgress(progress=progress, detail=paste0("Fitting model ",count," of ",length(datalist),"."), n=length(datalist), shiny=shiny, print=isTRUE(printProgress))

    n <- nrow(x)
    #If the weighs for this particular protein are of length 1, duplicate them to the correct length
    if(length(y)==1){y <- rep(y,n)}

    parsedFormula <- .createParsedFormula(x,y,response,fixed,shrinkage.fixed,random,formula_ran,add.intercept,formula_fix,scaleUnshrFix,modfiedGS)

    #For models that cannot be fitted, factors with only one level should not be removed (as is being done in .createParsedFormula), as the model cannot be fitted anyways!
    predictors2 <- unique(unlist(strsplit(c(fixed,random),":")))
    ridgeModel <-

      tryCatch(

        #Return a fitted lmerMod model with M estimation or fitted with given weights or NULL weights
        .ridgeEstM3(parsedFormula,
                    k = k,
                    robustM = robustM,
                    tolPwrss = tolPwrss,
                    verbose = verbose, ...)
        ,
        #If the fitting fails, return an empty lmerMod object
        error=function(e){
          .emptylmerMod(parsedFormula,x,y,predictors2,response,intercept,intercept_name,n)
        })

    #Pass on attributes from parsedFormula
    attr(ridgeModel,"MSqRob_R_fix") <- attr(parsedFormula,"MSqRob_R_fix")
    attr(ridgeModel,"MSqRob_Q_fix") <- attr(parsedFormula,"MSqRob_Q_fix")
    attr(ridgeModel,"MSqRob_Zt_indices") <- attr(parsedFormula,"MSqRob_Zt_indices")
    attr(ridgeModel,"MSqRob_unshr_pos") <- attr(parsedFormula,"MSqRob_unshr_pos")
    attr(ridgeModel,"MSqRob_levels") <- attr(parsedFormula,"MSqRob_levels")
    attr(ridgeModel,"MSqRob_cnms") <- attr(parsedFormula,"MSqRob_cnms")
    attr(ridgeModel,"MSqRob_Xp") <- attr(parsedFormula,"MSqRob_Xp")
    attr(ridgeModel,"MSqRob_Gp") <- attr(parsedFormula,"MSqRob_Gp")
    attr(ridgeModel,"MSqRob_XGp") <- attr(parsedFormula,"MSqRob_XGp")
    attr(ridgeModel,"MSqRob_fullcnms") <- attr(parsedFormula,"MSqRob_fullcnms")

    return(ridgeModel)
  }, datalist, weights, SIMPLIFY = FALSE)
  return(modellist)
}




.createParsedFormula=function(x,y,response,fixed,shrinkage.fixed,random,formula_ran,add.intercept,formula_fix,scaleUnshrFix,modfiedGS,nObsPep=NULL,nNonObsPep=NULL,v=NULL,familyfun=NULL){

  #,nObsPep=NULL,nNonObsPep=NULL,v=NULL need to be given if a count formula needs to be made, even though these arguments are never used,
  #they need to be set in order for the formula not to fail!

  n <- nrow(x)

  x <- droplevels(x)

  error <- FALSE

  if(n==0){error <- TRUE}

  # if(!is.null(fixed) && any(sapply(strsplit(fixed, ":"), function(z){return(length(levels(factor(do.call(paste, x[,z, drop=FALSE])))))})==1)){
  #   warning("Fixed effects should have at least 2 levels before a model can be fit.")
  #   error <- TRUE
  #   #AANPASSEN: die factor eruit gooien i.p.v. error te geven: is hetzelfde als van 3 naar 2 levels gaan.
  # }

  #Remove fixed effects with only one level so that these models can be fit (interpretation is done when fitting contrasts!)
  templist <- .removeSingleLevels(fixed, shrinkage.fixed, formula_fix, x, add.intercept)
  fixed <- templist[[1]]
  shrinkage.fixed <- templist[[2]]
  formula_fix <- templist[[3]]

  #Default: shrink all fixed effects together except the intercept
  if(!is.null(fixed) && is.null(shrinkage.fixed)){
    shrinkage.fixed <- rep(1, length(fixed))
    if(isTRUE(add.intercept)){ #If shrinkage.fixed is NULL AND there is an intercept!!!
      shrinkage.fixed <- c(0,shrinkage.fixed)
    }
  }

  fixed2 <- fixed

  if(isTRUE(add.intercept)){
    fixed2 <- c(1,fixed)
    #Always add 0 to fixed2, needed for the formula
  } else{fixed2 <- c(0,fixed)}

  #If there are no shrunken fixed effects
  if(all(shrinkage.fixed==0) || isTRUE(error)){
    ridgeGroupsForm <- NULL
    fixedUnshrForm <- paste0(fixed2, collapse="+")
    fixedUnshrGroups <- NULL

    #If there is at least one shrunken fixed effect
  } else if(any(shrinkage.fixed!=0)){

    #We don't want any singularities in our fixed effects!!!
    #Remove those singularities from the data!
    oldContr <- options("contrasts")$contrasts
    newContr <- oldContr
    newContr["unordered"] <- "contr.ltfr_MSqRob"
    options(contrasts = newContr)

    MMFull <- model.matrix(formula_fix, data=x)

    options(contrasts = oldContr)

    if(any(duplicated(cov(MMFull), MARGIN = 2))){
    lol <- lapply(adjustNames(list(x), fixed), function(z) {
      colnames(MMFull)[duplicated(cov(MMFull), MARGIN = 2)]==z}
      )

    x <- x[which(rowSums(lol[[1]])==0),]
    x <- droplevels(x)
    n <- nrow(x)
    y <- y[which(rowSums(lol[[1]])==0)]

    ###COPY PASTE VAN HIERBOVEN!!!####
    ###Recall functie???###

    #Remove fixed effects with only one level so that these models can be fit (interpretation is done when fitting contrasts!)
    templist <- .removeSingleLevels(fixed, shrinkage.fixed, formula_fix, x, add.intercept)
    fixed <- templist[[1]]
    shrinkage.fixed <- templist[[2]]
    formula_fix <- templist[[3]]

    #Default: shrink all fixed effects together except the intercept
    if(!is.null(fixed) && is.null(shrinkage.fixed)){
      shrinkage.fixed <- rep(1, length(fixed))
      if(isTRUE(add.intercept)){ #If shrinkage.fixed is NULL AND there is an intercept!!!
        shrinkage.fixed <- c(0,shrinkage.fixed)
      }
    }

    fixed2 <- fixed

    if(isTRUE(add.intercept)){
      fixed2 <- c(1,fixed)
      #Always add 0 to fixed2, needed for the formula
    } else{fixed2 <- c(0,fixed)}

    ############################
    ############################
    ############################

    }

    #Make fixed design model matrix MM, this step is only because we need its attributes
    #If you find a bug here, check options("contrasts")$contrasts: unordered should be "contr.treatment" and definitely not "contr.ltfr_MSqRob"!

    #Catch error in case x would have 0 rows after removing singularities
    if(n==0){
      MM <- x
    #Default model matrix for the fixed effects
    } else{MM <- model.matrix(formula_fix, data=x)}

    #Make fixed design matrix
    #design_fix <- Matrix::Matrix(MM)

    #Wat we nodig hebben:
    pred_assigned <- attr(MM, "assign")
    pred_names <- attr(MM, "dimnames")[[2]]

    groups <- unique(shrinkage.fixed[shrinkage.fixed!=0]) #select the groups that should get shrinkage
    indices <- which(shrinkage.fixed %in% groups)
    shrink <- pred_assigned %in% unique(pred_assigned)[indices]

    #Decompose fixed design matrix that needs shrinkage
    #qr_matrix <- Matrix::qr(design_fix[,shrink])
    # Q_fix=tryCatch(Matrix::Matrix(Matrix::qr.Q(qr_matrix), sparse=TRUE), error=function(e){
    #   return(Matrix(nrow=n, ncol=0, sparse=TRUE))
    # })
    # R_fix=suppressWarnings(tryCatch(Matrix::Matrix(Matrix::qr.R(qr_matrix), sparse=TRUE), error=function(e){ #qrR lukt niet altijd!
    #   return(Matrix(nrow=0, ncol=0, sparse=TRUE))
    # }))

    if(!isTRUE(modfiedGS)){
    Q_fix=tryCatch(Matrix::qr.Q(Matrix::qr(MM)), error=function(e){
      return(matrix(nrow=n, ncol=0))
    })
    R_fix=tryCatch(Matrix::qr.R(Matrix::qr(MM)), error=function(e){
      return(matrix(nrow=0, ncol=0))
    })
    } else{
      R_fix=gramschmidt(MM)$R
      Q_fix=gramschmidt(MM)$Q
    }

    #If R_fix has more columns than rows, make it square by adding zeros
    #This happens when some predictors are perfectly linearly dependent
    QR_fix <- addZerosQR(Q=Q_fix, R=R_fix)
    Q_fix <- QR_fix[["Q"]]
    R_fix <- QR_fix[["R"]]

    # R_fix2 <- R_fix
    # Q_fix2 <- Q_fix

    #Intercept will always be first in the design matrix!!!!
    #Try e.g.: model.matrix(~ -1 + lab + 1 + condition, data=x)
    #Check if there is an intercept!!!
    # if(isTRUE(add.intercept)){
    #
    # R_fix2[1,] <- R_fix2[1,]/R_fix2[1,1]
    # Q_fix2[,1] <- 1 #Q_fix2[,1]/Q_fix2[1,1]
    #
    # R_fix <- R_fix2
    # Q_fix <- Q_fix2
    # }

    # R_fix <- NULL
    # Q_fix <- MM

    ridgeGroupsForm <- NULL
    if(length(groups)!=0){
      ridgeGroups <- vector()
      for(i in 1:length(groups)){
        index <- which(shrinkage.fixed==groups[i])
        fixed.names <- pred_names[pred_assigned %in% unique(pred_assigned)[index]]
        if(n>=length(fixed.names)){
          x[,paste0("ridgeGroup.", i)] <- factor(rep(fixed.names,length=n),levels=fixed.names) #moet gewoon als levels het aantal groepen hebben
        }else{ #The case where you have more levels in your ridge group than observations; this sometimes happens when you include mutliple numeric variables, set ridgeGroup to NA, which will make the fit invalid!
          x[,paste0("ridgeGroup.", i)] <- factor(rep(NA,length=n),levels=fixed.names) #moet gewoon als levels het aantal groepen hebben
        }
        ridgeGroups[i] <- paste0("ridgeGroup.", i)
      }
      ridgeGroupsForm <- paste0("+",paste0("(1 | ",ridgeGroups,")", collapse="+"))
    }

    #If there are fixed effects that should not be shrunken on top of fixed effects that should be shrunken:

    #If there are no unshrunken fixed effects:
    fixedUnshrForm <- "-1"

    unshr_pos <- which(!shrink)

    #If there are also unshrunken fixed effects next to shrunken fixed effects
    #position of unshrunken effects in Q_fix:
    if(any(shrinkage.fixed==0) && nrow(x)>0){
      fixedUnshrGroups <- vector()

      for(i in 1:length(unshr_pos)){
        unshrunken_fixed_name <- fixed2[unshr_pos[i]]
        if(unshrunken_fixed_name=="1"){unshrunken_fixed_name <- "Intercept_MSqRob"}
        x[[unshrunken_fixed_name]] <- Q_fix[,unshr_pos[i],drop=FALSE]
        fixedUnshrGroups[i] <- unshrunken_fixed_name
      }

      fixedUnshrForm <- paste0("-1+",paste0(fixedUnshrGroups,collapse="+"))
    }

  }

  if(length(formula_ran)==0){part_ran <- NULL
  } else{part_ran <- paste0(c("+",formula_ran[[2]]), collapse="")}

  formula <- formula(paste0(response,"~",paste0(c(fixedUnshrForm,ridgeGroupsForm,part_ran), collapse="")))

  #environment(formula) <- environment(weights)
  .weight. <- y
  e <- new.env(parent=environment(formula))
  assign(".weight.", y, envir=e)
  environment(formula) <- e

  parsedFormula <- suppressWarnings(tryCatch(lme4::lFormula(formula, data=x, weights=.weight., control=lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                                                                                                         check.nobs.vs.rankZ = "ignore",
                                                                                                                         check.nobs.vs.nRE="ignore", check.nlev.gtr.1 = "ignore")),
                                             error=function(e){
                                               #warning(e)  #Outputting errors as warnings is not compatible with Shiny!!!
                                               parsedFormula <- list(formula=formula, fr=x)
                                               attr(parsedFormula,"MSqRob_fail") <- TRUE
                                               return(parsedFormula)
                                             }))

  #Proposal for scaling
  if(isTRUE(scaleUnshrFix) && !is.null(fixedUnshrGroups)){

  # parsedFormula$fr[,fixedUnshrGroups] <- scale(parsedFormula$fr[,fixedUnshrGroups,drop=FALSE])
  # parsedFormula$X <- scale(parsedFormula$X)

    origX <- parsedFormula$X[1,"Intercept_MSqRob"]

    X1 <- as.matrix(parsedFormula$fr[,fixedUnshrGroups,drop=FALSE])
    X2 <- parsedFormula$X

    X1 <- X1/origX
    X1[,"Intercept_MSqRob"] <- 1
    attr(X1,"scaled:center") <- 0
    attr(X1,"scaled:scale") <- origX

    X2 <- X2/origX
    X2[,"Intercept_MSqRob"] <- 1
    attr(X2,"scaled:center") <- 0
    attr(X2,"scaled:scale") <- origX

    parsedFormula$fr[,fixedUnshrGroups] <- X1
    parsedFormula$X <- X2

  }

  if(isTRUE(attr(parsedFormula,"MSqRob_fail"))){error <- TRUE}

  #Execute only once: add missing columns in data frame (e.g. the ones that are now in ridgeGroup) to "fr" slot of parsedFormula
  #-> bad idea: fitting doesn't work anymore then!
  #parsedFormula$fr <- cbind(parsedFormula$fr,x[,!(colnames(x) %in% colnames(parsedFormula$fr)), drop=FALSE])

  if(any(shrinkage.fixed!=0)  && !isTRUE(error)){
    num_indices <- as.list(which(names(parsedFormula$reTrms$cnms) %in% ridgeGroups))
    Zt_indices <- unlist(lapply(as.list(num_indices), function(z){return((parsedFormula$reTrms$Gp+1)[z]:parsedFormula$reTrms$Gp[z+1])}))

    #Fixed effects together, peptide sequences together
    parsedFormula$reTrms <- try(within(parsedFormula$reTrms, {
      #cnms$ridgeGroups <- "fixed"
      if(length(unshr_pos)!=0){
        Zt[Zt_indices,] <- Matrix::t(Q_fix[,-unshr_pos,drop=FALSE]) #t(design_fix_shrink)
      } else{Zt[Zt_indices,] <- Matrix::t(Q_fix[,,drop=FALSE])}
    }), silent=TRUE)

    if(!inherits(parsedFormula$reTrms,'try-error')){
      #Pass Q_fix, R_fix and Zt_indices as attributes
      attr(parsedFormula,"MSqRob_R_fix") <- R_fix
      attr(parsedFormula,"MSqRob_Q_fix") <- Q_fix
      attr(parsedFormula,"MSqRob_Zt_indices") <- Zt_indices
      attr(parsedFormula,"MSqRob_unshr_pos") <- unshr_pos

      # oldContr <- options("contrasts")$contrasts
      # contr.ltfr_old <- tryCatch(contr.ltfr, error=function(e){
      #   return(contr.treatment)})
      #
      # contr.ltfr <- caret::contr.ltfr
      # environment(contr.ltfr) <- environment(caret::contr.ltfr)
      # newContr <- oldContr
      # newContr["unordered"] <- "contr.ltfr"
      # options(contrasts = newContr)

      dummies <- dummyVars_MSqRob(formula_fix, data=x)
      dummies$sep <- ""
      attr(parsedFormula,"MSqRob_levels") <-  c(colnames(predict(dummies, newdata=x)),rownames(parsedFormula$reTrms$Zt[-Zt_indices,]))

      #
      #colnames(model.matrix(dummies, data=x))

      # options(contrasts = oldContr)
      # contr.ltfr <- contr.ltfr_old

      if(isTRUE(add.intercept)){attr(parsedFormula,"MSqRob_levels") <- c("(Intercept)",attr(parsedFormula,"MSqRob_levels"))}


      ####Piece of code only needed for attributes####
      #Determine real cnms: replace ridgeGroups
      cnms <- parsedFormula$reTrms$cnms
      for(i in 1: length(groups)){

        index <- which(shrinkage.fixed==groups[i])

        target <- which(names(cnms)==paste0("ridgeGroup.", i))

        test <- rep("(Intercept)", length(index))
        names(test) <- fixed2[index]

        cnms <- c(cnms[(1:length(cnms))<target,drop=FALSE], test, cnms[(1:length(cnms))>target,drop=FALSE])
      }
      attr(parsedFormula,"MSqRob_cnms") <-  cnms

      #Determine real Gp (no ridgeGroups)
      #pred_assigned is about everything in Z

      #Give names to pred_assigned
      pred_counter <- pred_assigned[1]
      assign_counter <- 1
      for(k in 1: length(pred_assigned)){

        if(pred_counter==pred_assigned[k]){
        } else{
          assign_counter <- assign_counter+1
        }
        names(pred_assigned)[k] <- fixed[assign_counter]
        pred_counter <- pred_assigned[k]
      }

      test <- pred_assigned[shrink][order(match(names(pred_assigned[shrink]),names(cnms)))]+1
      Gp <- numeric()
      rest_indices <- Zt_indices

      for(k in 1:length(indices)){
        Gp[k] <- rest_indices[1]-1
        rest_indices <- rest_indices[-(1:sum(test==test[1]))]
        test <- test[-(1:sum(test==test[1]))]
      }

      Gp <- unique(sort(c(parsedFormula$reTrms$Gp, Gp)))
      attr(parsedFormula,"MSqRob_Gp") <-  Gp
      ##########

    }

    #If all fixed effects are zero, there still needs to be an attribute called "MSqRob_levels"
  } else if(!isTRUE(error)){
    dummies <- dummyVars_MSqRob(formula_fix, data=x)
    dummies$sep <- ""
    attr(parsedFormula,"MSqRob_levels") <-  c(colnames(predict(dummies, newdata=x)),rownames(parsedFormula$reTrms$Zt))

    if(isTRUE(add.intercept)){attr(parsedFormula,"MSqRob_levels") <- c("(Intercept)",attr(parsedFormula,"MSqRob_levels"))}
  }

  if(is.null(attr(parsedFormula,"MSqRob_Gp")) && !inherits(parsedFormula$reTrms,'try-error')){
    attr(parsedFormula,"MSqRob_cnms") <-  parsedFormula$reTrms$cnms
    attr(parsedFormula,"MSqRob_Gp") <-  parsedFormula$reTrms$Gp
  }

  Xp <- NULL
  fixed_X <- NULL
  if(!is.null(fixed) & !isTRUE(error)){
    #Unshrunken fixed effects

    fixed_X <- fixed2[shrinkage.fixed==0]
    fixed_X[fixed_X=="1"] <- "(Intercept)"

    #X_assigned is about everything in X
    X_assigned <- attr(parsedFormula$X,"assign")

    #Give names to X_assigned
    pred_counter <- X_assigned[1]
    assign_counter <- 1
    if(length(X_assigned)>0){
      for(k in 1:length(X_assigned)){

        if(pred_counter==X_assigned[k]){
        } else{
          assign_counter <- assign_counter+1
        }
        names(X_assigned)[k] <- fixed_X[assign_counter]
        pred_counter <- X_assigned[k]
      }
    }

    Xp <- which(!duplicated(names(X_assigned)))-1

    #If you want the "real" Xp:
    attr(parsedFormula,"MSqRob_Xp") <- c(Xp,ncol(parsedFormula$X))
  }

  attr(parsedFormula,"MSqRob_XGp") <- c(Xp,attr(parsedFormula,"MSqRob_Gp")+ncol(parsedFormula$X))
  attr(parsedFormula,"MSqRob_fullcnms") <- c(fixed_X,names(attr(parsedFormula,"MSqRob_cnms")))

  return(parsedFormula)
}


#This function adds zero columns to the Q matrix and zero rows to the R matrix in a QR decomposition so that both matrices become square.
#Can also be used if only an R matrix is provided
#Is only used when some predictors are perfectly linearly dependent
addZerosQR <- function(Q=NULL, R){
  #If R_fix has more columns than rows, make it square by adding zeros
  #Check before: Q_fix%*%R_fix
  missingcols <- diff(dim(R))
  if(missingcols!=0){
    if(!is.null(Q)){Q=cbind(Q, matrix(0, nrow=nrow(Q), ncol=missingcols))}
    R=rbind(R, matrix(0, nrow=missingcols, ncol=ncol(R)))
  }
  #Check after: Q_fix%*%R_fix (idem!)
  return(list(Q=Q, R=R))
}


#This function creates an empty LM object
.emptyLM=function(formula,x,weights,predictors,response,add.intercept,intercept_name,n){

  terms <- terms(formula)

  flist <- lapply(as.list(x[,predictors]),"factor")
  times <- vapply(flist, function(x){return(length(unique(x)))}, 0)

  coefs <- rep(NA,sum(times))

  if(isTRUE(add.intercept)){
    names(coefs) <- c(intercept_name,names(unlist(sapply(flist, function(x){return(unique(x)[-1])}))))
  }else {names(coefs) <- names(unlist(sapply(flist, function(x){return(unique(x))})))}

  emptyLM <- list(coefficients=coefs,residuals=rep(NA, n),fitted.values=rep(NA, n),effects=rep(NA, n),rank=sum(times)-1+as.numeric(add.intercept),weights=rep(1, n), qr=list(qr=matrix(NA, nrow(x), length(coefs)), qraux=rep(NA, length(coefs)), pivot=1:length(coefs), tol=NA, rank=length(coefs)), df.residual=n-(sum(times)-1+as.numeric(add.intercept)),call=call("lm(formula = formula, data = x, weights = weights)"),terms=terms,contrasts=NA,xlevels=NA,offset=NA,y=NA,x=NA,model=data.frame(x[c(response,predictors)], "(weights)"=rep(1, n), check.names=FALSE),na.action=NA)
  class(emptyLM) <- "lm"

  return(emptyLM)
}


#This function creates an empty lmerMod object
.emptylmerMod=function(formula,x,y,predictors,response,intercept,intercept_name,n){

  #Initialize variables specifically for ridge models when the model cannot be fit:
  beta <- as.numeric(NA)

  cnms <- rep(list("(Intercept)"),length(predictors))
  names(cnms) <- predictors
  flist <- lapply(as.list(x[,predictors,drop=FALSE]),"factor")
  attr(flist, "assign") <- seq(1,length(predictors))
  nrandom <- sum(sapply(predictors, function(y){return(length(unique(do.call('$',list(x,y)))))}))
  y <- do.call('$', list(x,response))
  X <- matrix(1,nrow=length(y), dimnames=list(NULL,intercept_name))
  Lind <- rep(1:length(predictors),times=vapply(flist, function(y){return(length(unique(y)))},0))
  emptylmerMod <- new("lmerMod", resp=new("lmerResp",Ptr="<externalptr>", mu=rep(NA,n),offset=rep(0,n),sqrtXwt=rep(1,n),sqrtrwt=rep(1,n),weights=rep(1,n),y=y),Gp=as.integer(0),call=call(paste0("lme4::lmer(formula = ","formula",", data = x, weights = weights)")),frame=x,flist=flist,cnms=cnms,lower=as.numeric(NA),theta=as.numeric(rep(NA, length(cnms))),beta=beta,u=rep(as.numeric(NA),nrandom),devcomp=list(dims=c(N=n,n=length(y),p=as.numeric(intercept),nmp=0,nth=1,q=nrandom,nAGQ=NA,compDev=TRUE,useSc=TRUE,reTrms=1,REML=1,GLMM=FALSE,NLMM=FALSE),cmp=c(ldL2=NA,ldRX2=NA,wrss=NA,ussq=NA,pwrss=NA,drsum=NA,REML=1,dev=NA,sigmaML=NA,sigmaREML=NA,tolPwrss=NA)),pp=new("merPredD",X=X,Zt=do.call("rbind", lapply(x[,predictors,drop=FALSE], Matrix::fac2sparse)),Lambdat=Matrix::Matrix(diag(nrandom),sparse=TRUE),Lind=Lind,theta=as.numeric(rep(NA, length(cnms))),n=1),optinfo=list(optimizer="bobyqa",control=list(iprint=0),derivs=list(gradient=NA,Hessian=NA),conv=list(opt=0,lme4=list()),feval=NA,warnings=list(),val=NA))
  return(emptylmerMod)
}

#' Adjust the names of the elements in a list of dataframes
#'
#' @description Given a list of dataframes, this function pastes the colnames of all factors and characters to the elements in these columns.
#' This is necessary to produces consisten results when fitting with lm or lmer. This function is borderline internal.
#' @param datalist A list of dataframes.
#' @param predictors A vector of predictors corresponding to some of the columnames in the list of dataframes for which the given adjustment will be performed if the columns corresponding to these predictors are characters or factors.
#' @return The list of dataframes with adusted elements.
#' @export
adjustNames=function(datalist, predictors){

  if(!is.null(predictors)){
  datalist_adj <- lapply(datalist, function(x){

    for(i in 1:length(x[,predictors, drop=FALSE]))
    {
      if(is.factor(x[,predictors, drop=FALSE][[i]]) || is.character(x[,predictors, drop=FALSE][[i]]))
      {
        x[,predictors[i]] <- factor(paste0(colnames(x[,predictors, drop=FALSE])[i],x[,predictors, drop=FALSE][[i]]), levels=unique(paste0(colnames(x[,predictors, drop=FALSE])[i],x[,predictors, drop=FALSE][[i]])))
      }
    }

    return(x)
  })
  } else{datalist_adj <- datalist}
  return(datalist_adj)
}


#Function for ridge regression
.ridgeEstM3=function(parsedFormula, k, robustM = robustM, tolPwrss = 1e-10, verbose=FALSE, ...)
{

  #If ridge regression with Huber weights:
  if(!is.null(parsedFormula$fr$`(weights)`[1]) && parsedFormula$fr$`(weights)`[1]=="Huber"){ridgeFit <- .huberRidgeFit(parsedFormula, k, robustM = robustM, tolPwrss = tolPwrss, verbose=verbose, ...)

  #Additional weighing options could be added here using if else

  #Else: ridge regression with given weights (or NULL):
  } else{ridgeFit <- .weightsRidgeFit(parsedFormula,verbose=verbose,...)}

  return(ridgeFit)
}

#Function for ridge regression with given weights or NULL weights
.weightsRidgeFit=function(parsedFormula, verbose=FALSE, ...){

  #!!!TryCatch aanpasssen enzo!!!
  devianceFunction <- tryCatch(do.call(lme4::mkLmerDevfun, parsedFormula), error=function(e){
    return(NULL)})

  optimizerOutput <- suppressWarnings(tryCatch(lme4::optimizeLmer(devianceFunction), error=function(e){return(NULL)}))

  #ridgeFit=lme4::lmer(formula, data=data,weights=weights,...)
  ridgeFit=lme4::mkMerMod(
    rho = environment(devianceFunction),
    opt = optimizerOutput,
    reTrms = parsedFormula$reTrms,
    fr = parsedFormula$fr)

  attr(ridgeFit,"MSqRob_lmer_Deviance") <- environment(devianceFunction)$lmer_Deviance

  if (isTRUE(verbose)) print(ridgeFit)

  return(ridgeFit)
}

#Function for Ridge regression with M estimation
.huberRidgeFit=function(parsedFormula, k, robustM=TRUE, tolPwrss = 1e-10, verbose=FALSE, ...){

  #Remove "Huber" from the `(weights)` column
  parsedFormula$fr$`(weights)` <- NULL

  pwrssnew <- Inf
  weights <- NULL
  start <- NULL ###NEW!!!

  devianceFunction <- tryCatch(do.call(lme4::mkLmerDevfun, parsedFormula), error=function(e){
    return(NULL)})

  n <- nrow(parsedFormula$X)
  p <- ncol(parsedFormula$X)

  repeat{
    pwrssold <- pwrssnew

    #!!!TryCatch aanpasssen enzo!!!
    parsedFormula$start <- start ###NEW!!!

    environment(devianceFunction)$resp <- lme4:::mkRespMod(parsedFormula$fr, REML = p)

    #devianceFunction <- lme4:::mkdevfun(environment(devianceFunction), 0L)

    #environment(devianceFunction)$pp$updateDecomp()
    #environment(devianceFunction)$pp$updateL()
    #environment(devianceFunction)$pp$updateLamtUt()
    #environment(devianceFunction)$pp$updateRes(wtres=resid)
    #environment(devianceFunction)$pp$updateXwts(wts=MASS::psi.huber(resid/sigma,k=k))

    #devianceFunction(environment(devianceFunction)$pp$theta) # one evaluation to ensure all values are set

    optimizerOutput <- suppressWarnings(tryCatch(lme4::optimizeLmer(devianceFunction), error=function(e){return(NULL)})) #, optimizer = "nloptwrap"

    #ridgeFit=lme4::lmer(formula, data=data,weights=weights,...)
    # ridgeFit=lme4::mkMerMod(
    #   rho = environment(devianceFunction),
    #   opt = optimizerOutput,
    #   reTrms = parsedFormula$reTrms,
    #   fr = parsedFormula$fr)

    ###NEW

    pp <- environment(devianceFunction)$pp
    resp <- environment(devianceFunction)$resp

    sqrLenU <- pp$sqrL(1)
    wrss    <- resp$wrss()

    #pwrss
    pwrssnew <- wrss + sqrLenU
    #wrss + sqrLenU equals:
    #sd(test@resp$wtres)^2*(n-1)+mean(test@resp$wtres)^2*n+test@pp$sqrL(1)

    if(pwrssold-pwrssnew <= tolPwrss){
      break
    }

    # n <- nrow(pp$V)
    # p <- ncol(pp$V)

    sigmaML <- pwrssnew/n #Or a robust derivative of sigmaML...!!!
    #sigma
    sigma <- sqrt(unname(sigmaML*(n/(n-p)))) #Or a robust derivative of sigma...!!!
    #theta
    theta <- pp$theta

    #residuals
    resid <- residuals(resp, "response")

    ###

    #pwrssnew <- lme4::getME(ridgeFit,"devcomp")$cmp["pwrss"]

    #parsedFormula$fr$`(weights)` <- MASS::psi.huber(resid(ridgeFit)/sigma(ridgeFit),k=k)
    #parsedFormula$fr$`(weights)` <- MASS::psi.huber(resid/sigma,k=k)

    #Niet correct: environment(devianceFunction)$pp$updateXwts(wts=MASS::psi.huber(resid/sigma,k=k))
    #environment(devianceFunction)$resp$setWeights(ww=MASS::psi.huber(resid/sigma,k=k))

    if(isTRUE(robustM)){
      parsedFormula$fr$`(weights)` <- MASS::psi.huber(resid/mad(resid),k=k)
    }else{
      parsedFormula$fr$`(weights)` <- MASS::psi.huber(resid/sigma,k=k)
    }

    #start <- lme4::getME(ridgeFit, "theta")
    start <- theta
  }

  ridgeFit=lme4::mkMerMod(
    rho = environment(devianceFunction),
    opt = optimizerOutput,
    reTrms = parsedFormula$reTrms,
    fr = parsedFormula$fr)

  attr(ridgeFit,"MSqRob_lmer_Deviance") <- environment(devianceFunction)$lmer_Deviance

  if (isTRUE(verbose)) print(ridgeFit)

  return(ridgeFit)
}

#Function for linear model with M estimation
.huberLmFit=function(formula, data, k=1.345, robustM=TRUE, tolWrss = 1e-10, verbose=FALSE,...)
{
  formula <- as.formula(formula)
  wrssnew <- Inf
  .weight. <- NULL
  repeat{
    e2 <- new.env(parent=environment(formula))
    assign(".weight.",.weight., envir=e2)
    environment(formula) <- e2
    wrssold <- wrssnew
    lmFit=lm(formula, data=data,weights=.weight.)
    summary <- summary(lmFit)
    if(isTRUE(robustM)){
      .weight.=MASS::psi.huber(summary$residuals/mad(summary$residuals),k=k) #The new, robust implementation
      #Note:
      #summary(test)$sigma==(sd(residuals)*sqrt(length(residuals)-1))/sqrt(test$df.residual)
    } else {.weight.=MASS::psi.huber(summary$residuals/summary$sigma,k=k)} #The old implementation
    #Weighted residual sum of squares
    wrssnew <- sum(summary$residuals^2*.weight.) #136.4315

    #The "is.na" is for the cases where some, but not all factors are NA => model fitted but no residual dfs => no sigma
    if(is.na(wrssold-wrssnew) || (wrssold-wrssnew <= tolWrss)){
      break
    }
  }
  if (verbose) print(lmFit)

  return(lmFit)
}

gramschmidt <- function(x) {
  x <- as.matrix(x)
  # Get the number of rows and columns of the matrix
  n <- ncol(x)
  m <- nrow(x)

  # Initialize the Q and R matrices
  q <- matrix(0, m, n)
  r <- matrix(0, n, n)

  for (j in 1:n) {
    v = x[,j] # Step 1 of the Gram-Schmidt process v1 = a1
    # Skip the first column
    if (j > 1) {
      for (i in 1:(j-1)) {
        r[i,j] <- t(q[,i]) %*% x[,j] # Find the inner product (noted to be q^T a earlier)
        # Subtract the projection from v which causes v to become perpendicular to all columns of Q
        v <- v - r[i,j] * q[,i]
      }
    }
    # Find the L2 norm of the jth diagonal of R
    r[j,j] <- sqrt(sum(v^2))
    # The orthogonalized result is found and stored in the ith column of Q.
    q[,j] <- v / r[j,j]
  }

  # Collect the Q and R matrices into a list and return
  qrcomp <- list('Q'=q, 'R'=r)
  return(qrcomp)
}



