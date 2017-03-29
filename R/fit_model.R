#' Fit peptide-based models
#'
#' @description Fits a model to each protein of a \code{\link[=protdata-class]{protdata}} object and returns a corresponding \code{\link[=protLM-class]{protLM}} object.
#' In its standard settings, the function returns a \code{\link[=protLM-class]{protLM}} object containing robust ridge models as described in Goeminne et al. (2015).
#' However, the the user can also specify to turn off the ridge regression and fit the models by ordinary least squares (OLS) and/or to turn off the down-weighing of outliers by M estimation with Huber weights and/or to turn off the Empirical Bayes squeezing of variances.
#' @param protdata A \code{\link[=protdata-class]{protdata}} object to which peptide-based models must be fitted. Note that factors should be coded as factors and numeric variables as numeric in each data slot.
#' @param response The name of the column in the data slot of the \code{\link[=protdata-class]{protdata}} object that contains the response variable for the model, mostly the column containing the log2 transformed intensity values.
#' @param fixed Either a vector of names corresponding to the columns in the data slot of the \code{\link[=protdata-class]{protdata}} object containing the predictor variables, or a right-hand sided fixed effects formula without intercept, which should be indicated in the argument \code{add.intercept}. \code{NULL} (the default) indicates that no fixed effects other than a possible fixed intercept are present in the model.
#' @param random Either a vector of names corresponding to the columns in the data slot of the \code{\link[=protdata-class]{protdata}} object containing the predictor variables or a right-hand sided random effects formula. Adding the peptide sequences as one of the random effects predictors is highly recommended as individual peptide effects are often quite strong. \code{NULL} (the default) indicates that no random effects are present in the model.
#' @param add.intercept A logical value indicating whether the fitted models should contain a fixed intercept. If missing, the value is set to \code{TRUE}, indicating the intercept should be present in the model.
#' @param shrinkage.fixed A numeric vector containing only 0 and/or 1 of length equal to the number of fixed effects, potential intercept included. The nth element in the shrinkage.fixed vector indicates whether the nth fixed effect should be shrunken (value 1) or not (value 0). If \code{add.intercept=TRUE}, the first element of the vector indicates the intercept. \code{shrinkage.intercept = NULL} (default) indicates all fixed effects except the intercept should be shrunken.
#' @param weights The type of weighing that should be performed. Supported weighing methods incluce \code{"Huber"} (the default) for M estimation with Huber weights and \code{NULL} when no weighing should be applied. One can also supply a list of weights with length equal to the number of proteins in the \code{\link[=protdata-class]{protdata}} object. Each element of the list should either contain \code{"Huber"} or \code{NULL} or a numeric vector containing weights with length equal to the number of observations for that protein.
#' @param k The tuning constant for the Huber mean when weighing down outliers. The default (\code{k = 1.345}) produces 95 \% efficiency relative to the sample mean when the population is normal but provides substantial resistance to outliers when it is not.
#' @param par_squeeze Character vector indicating which model parameters need to be squeezed. When squeezing random effects, provide their names. Fixed effects are present in shrinkage groups, e.g. ridgeGroup.1. If you want them to be squeezed as well, provide the names of the shrinkage groups that need to be squeezed. The default \code{NULL} indicates that no parameters will be squeezed.
#' @param squeezeVar A logical indicating whether the residual standard deviation of all models should be squeezed towards a common value. Defaults to \code{TRUE}. If set to \code{FALSE}, residual standard deviations will not be squeezed.
#' @param min_df A numeric value indicating the minimal degrees of freedom that will be taken into account for calculating the prior degrees of freedom and prior variance. Only used when \code{par_squeeze=TRUE} or \code{squeezeVar} is not \code{NULL}.
#' @param robust_var A logical indicating wheter the estimation of the prior degrees of freedom and the prior variance (needed for shrinkage) should be robustified against outlier variances. Only used when \code{par_squeeze=TRUE} or \code{squeezeVar} is not \code{NULL}. Defaults to \code{TRUE}.
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
fit.model=function(protdata, response=NULL, fixed=NULL, random=NULL, add.intercept=TRUE, shrinkage.fixed=NULL, weights="Huber", k = 1.345, par_squeeze=NULL, squeezeVar=TRUE, min_df=1, robust_var=TRUE, tolPwrss = 1e-10, verbose=FALSE, printProgress=FALSE, shiny=FALSE, message_fitting=NULL, message_thetas=NULL, message_squeeze=NULL, message_update=NULL, ...)
{

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message_fitting, value = 0)
  }

  #Control: fixed en random connot be NULL at the same time!
  if(is.null(response)){stop("Please specify a response variable.")}
  if(is.null(fixed) && is.null(random)){stop("Please specify appropriate fixed and/or random effects.")}

  fixed_input <- makeFormulaPredictors(fixed, add.intercept, effect="fixed")
  random_input <- makeFormulaPredictors(random, add.intercept, effect="random")

  fixed <- fixed_input[[1]]
  random <- random_input[[1]]

  formula_fix <- fixed_input[[2]]
  formula_ran <- random_input[[2]]

  #Error control: fixed and random effects should be present as colnames in the data slot
  #Random effects only checked if random is not NULL (see below)
  possible_vars <- colnames(getData(protdata, simplify=FALSE)[[1]])

  #Correct for possible interactions (only used in error control that follows)
  fixed_all <- unique(unlist(strsplit(fixed,":")))
  if(!all(fixed_all %in% possible_vars)){
    not_present <- which(!(fixed_all %in% possible_vars))
    stop(paste0("The following fixed effects are no possible predictors: \"",paste0(fixed_all[not_present], collapse="\", \""),"\". Please choose from the following predictors: \"",paste0(possible_vars, collapse="\", \""),"\"."))
  }

  if(length(response)!=1){stop("Please provide exactly one response variable.")}
  if(!(response %in% possible_vars)){
    stop(paste0("The response variable should be one of: \"",paste0(possible_vars, collapse="\", \""),"\"."))
  }

  #If an intercept is specified in either the fixed or the random part, an intercept should be added!
  #add.intercept <- fixed_input[[3]] || random_input[[3]]

  ### Initialize variables ###

  #Variables for when the model cannot be fit:
  intercept_name <- "(Intercept)"
  intercept <- TRUE

  modellist <- vector("list", length(getAccessions(protdata)))
  data <- getData(protdata,simplify=FALSE)

  #######

  #If weights is a single element, turn it into a list
  if(!is.list(weights)){weights <- rep(list(weights), length(protdata))}

  #Check if length of weights equals length of protdata object
  if(length(weights)!=length(protdata)){stop("The length of list \"weights\" should be equal to the number of proteins in the dataset.")}

  if(!is.null(random)){

    random_all <- unique(unlist(strsplit(random,":")))
    if(!all(random_all %in% possible_vars)){
      not_present <- which(!(random_all %in% possible_vars))
      stop(paste0("The following random effects are no possible predictors: \"",paste0(random_all[not_present], collapse="\", \""),"\". Please choose from the following predictors: \"",paste0(possible_vars, collapse="\", \""),"\"."))
    }

    #Initialize variables specifically for ridge models
    beta <- as.numeric(NA)
    intercept_val <- 1

    #Adjust names of predictors to make them similar to those of "lm"
    #In case of factors, the name of the factor is concatenated with the level of the factor
    data <- .adjustNames(data, random)

    #Fit the list of ridge models
    modellist <- .createRidgeList(data,weights,response,fixed,shrinkage.fixed,formula_fix,random,formula_ran,add.intercept,intercept,intercept_val,intercept_name,k,tolPwrss,beta,verbose,progress=progress,printProgress=printProgress,shiny=shiny,...)

  } else if(is.null(random)){

    #Initialize variables specifically for OLS models
    formula <- formula(paste0(response,paste0(as.character(formula_fix), collapse="")))
    intercept_val <- NA
    terms <- formula
    attributes(terms) <- list(variables=c(response, fixed), factors=rbind(rep(0,length(fixed)),diag(1, length(fixed))), term.labels=fixed, order=rep(1, length(fixed)), intercept=1, response=1)

    #Fit the list of lm models
    modellist <- .createLmList(data,weights,formula,fixed,response,fixed,intercept,intercept_val,intercept_name,terms,k,tolPwrss,verbose,progress=progress,printProgress=printProgress,shiny=shiny,...)

  }

  protLM <- new("protLM", accession=protdata@accession, model=modellist, annotation=protdata@annotation)
  protLM <- squeezePars(protLM, par_squeeze=par_squeeze, squeezeVar=squeezeVar, min_df=min_df, robust_var=robust_var, printProgress=printProgress, shiny=shiny, message_thetas=message_thetas, message_squeeze=message_squeeze, message_update=message_update)

  return(protLM)

}

removeSingleLevels <- function(fixed, shrinkage.fixed, formula_fix, x, add.intercept){

  morelvls <- getMoreLevels(fixed, x)
  fixed <- fixed[morelvls]

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

  formula_fix <- paste0(c("~",rep("",length(terms)-1)),terms, collapse="+")

  return(list(fixed, shrinkage.fixed, as.formula(formula_fix)))
}

getMoreLevels <- function(predictors, x){
  #return(unlist(lapply(predictors, function(y) length(levels(x[,y]))!=1)))
  if(is.null(predictors)){return(FALSE)
    }else{
  return(unlist(lapply(strsplit(predictors,":"), function(y) prod(unlist(lapply(x[,y,drop=FALSE], function(z) prod(length(unique(z))))))!=1)))
    }
}


#A function to turn predictor variables into a formula
makeFormulaPredictors <- function(input, intercept, effect){

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

    #If intercept==TRUE, add an intercept to the fixed effects, don't do it for random effects
    if(isTRUE(intercept) && effect=="fixed"){
      formula <- "~1"
      #We return our predictors WITHOUT a possible intercept term
      predictors <- NULL
    #If intercept==FALSE, no intercept
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

      #If intercept==TRUE, add an intercept
      if(isTRUE(intercept)){input <- c(1,input)
      #If intercept==FALSE, no intercept
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
.createLmList=function(data,weights,formula,predictors,response,fixed,intercept,intercept_val,intercept_name,terms,k,tolPwrss,verbose,progress=NULL,printProgress=NULL,shiny=FALSE,...){

  count <- 0

  #Return a fitted lmerMod model with M estimation or fitted with given weights or NULL weights
  mapply(function(x,y){

    count <<- count+1
    updateProgress(progress=progress, detail=paste0("Fitting model ",count," of ",length(data),"."), n=length(data), shiny=shiny, print=isTRUE(printProgress))

    n <- nrow(x)
    #If the weighs for this particular protein are of length 1, duplicate them to the correct length
    if(length(y)==1){y <- rep(y,n)}

    return(
      tryCatch(
        .lmEstM3(formula=formula,data=x,k=k,weights=y,tolWrss=tolPwrss,verbose=verbose,...)
        , error=function(e){
          return(.emptyLM(data,formula,x,y,predictors,response,intercept,intercept_val,intercept_name,n,terms))})
    )

  }, data, weights, SIMPLIFY = FALSE)

}

#Function for OLS regression
.lmEstM3=function(formula,data,k,weights,tolWrss = 1e-10,verbose=FALSE,...)
{

  #If OLS regression with Huber weights:
  if(!is.null(weights) && weights[1]=="Huber"){lmFit <- .huberLmFit(formula,data,k,tolWrss,verbose,...)}

  #Additional weighing options could be added here

  #Else: OLS regression with given weights (or NULL):
  else{lmFit <- .weightsLmFit(formula,data,k,weights,verbose,...)}

  return(lmFit)
}


.weightsLmFit=function(formula,data,k,weights,verbose,...)
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
.createRidgeList=function(data,weights,response,fixed,shrinkage.fixed,formula_fix,random,formula_ran,add.intercept,intercept,intercept_val,intercept_name,k,tolPwrss,beta,verbose,progress=NULL,printProgress=NULL,shiny=FALSE,...){

  count <- 0

  modellist <- mapply(function(x,y){

    count <<- count+1
    updateProgress(progress=progress, detail=paste0("Fitting model ",count," of ",length(data),"."), n=length(data), shiny=shiny, print=isTRUE(printProgress))

    n <- nrow(x)
    #If the weighs for this particular protein are of length 1, duplicate them to the correct length
    if(length(y)==1){y <- rep(y,n)}

      parsedFormula <- .createParsedFormula(x,y,response,fixed,shrinkage.fixed,random,formula_ran,add.intercept,formula_fix)

      #For models that cannot be fitted, factors with only one level should not be removed (as is being done in .createParsedFormula), as the model cannot be fitted anyways!
      predictors2 <- unique(unlist(strsplit(c(fixed,random),":")))
      ridgeModel <-

         tryCatch(

        #Return a fitted lmerMod model with M estimation or fitted with given weights or NULL weights
        .ridgeEstM3(parsedFormula,
                    k = k,
                    tolPwrss = tolPwrss,
                    verbose = verbose, ...)
       ,
        #If the fitting fails, return an empty lmerMod object
         error=function(e){
           .emptylmerMod(data,parsedFormula,x,y,predictors2,response,intercept,intercept_val,intercept_name,n,beta)
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
  }, data, weights, SIMPLIFY = FALSE)
  return(modellist)
}




.createParsedFormula=function(x,y,response,fixed,shrinkage.fixed,random,formula_ran,add.intercept,formula_fix){

  n <- nrow(x)

  x <- droplevels(x)

  error <- FALSE

  # if(!is.null(fixed) && any(sapply(strsplit(fixed, ":"), function(z){return(length(levels(factor(do.call(paste, x[,z, drop=FALSE])))))})==1)){
  #   warning("Fixed effects should have at least 2 levels before a model can be fit.")
  #   error <- TRUE
  #   #AANPASSEN: die factor eruit gooien i.p.v. error te geven: is hetzelfde als van 3 naar 2 levels gaan.
  # }

  #Remove fixed effects with only one level so that these models can be fit (interpretation is done when fitting contrasts!)
  templist <- removeSingleLevels(fixed, shrinkage.fixed, formula_fix, x, add.intercept)
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

  #If there is at least one shrunken fixed effect
  } else if(any(shrinkage.fixed!=0)){

  #Make fixed design model matrix MM, this step is only because we need its attributes
  MM <- model.matrix(formula_fix, data=x)

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

  Q_fix=tryCatch(Matrix::qr.Q(Matrix::qr(MM)), error=function(e){
    return(matrix(nrow=n, ncol=0))
  })
  R_fix=tryCatch(Matrix::qr.R(Matrix::qr(MM)), error=function(e){
    return(matrix(nrow=0, ncol=0))
  })

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

  #position of unshrunken effects in Q_fix:
  fixedUnshrForm <- NULL

  if(any(shrinkage.fixed==0)){
  unshr_pos <- which(!shrink)
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

  formula <- formula(paste0(response,"~",paste0(c(fixedUnshrForm,ridgeGroupsForm,"+",formula_ran[[2]]), collapse="")))

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
                              parsedFormula <- formula
                              attr(parsedFormula,"MSqRob_fail") <- TRUE
                              return(parsedFormula)
                            }))

  if(isTRUE(attr(parsedFormula,"MSqRob_fail"))){error <- TRUE}

  if(any(shrinkage.fixed!=0)  && !isTRUE(error)){
  num_indices <- as.list(which(names(parsedFormula$reTrms$cnms) %in% ridgeGroups))
  Zt_indices <- unlist(lapply(as.list(num_indices), function(z){return((parsedFormula$reTrms$Gp+1)[z]:parsedFormula$reTrms$Gp[z+1])}))

  #Fixed effects together, peptide sequences together
  parsedFormula$reTrms <- try(within(parsedFormula$reTrms, {
    #cnms$ridgeGroups <- "fixed"
    Zt[Zt_indices,] <- Matrix::t(Q_fix[,-unshr_pos,drop=FALSE]) #t(design_fix_shrink)
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
    names(pred_assigned)[k] <- fixed2[assign_counter]
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
    for(k in 1: length(X_assigned)){

      if(pred_counter==X_assigned[k]){
      } else{
        assign_counter <- assign_counter+1
      }
      names(X_assigned)[k] <- fixed_X[assign_counter]
      pred_counter <- X_assigned[k]
    }

    Xp <- which(!duplicated(names(X_assigned)))-1

    #If you want the "real" Xp:
    attr(parsedFormula,"MSqRob_Xp") <- c(Xp,ncol(parsedFormula$X))
  }

  attr(parsedFormula,"MSqRob_XGp") <- c(Xp,attr(parsedFormula,"MSqRob_Gp")+ncol(parsedFormula$X))
  attr(parsedFormula,"MSqRob_fullcnms") <- c(fixed_X,names(attr(parsedFormula,"MSqRob_cnms")))

  return(parsedFormula)
}



#This function creates an empty LM object
.emptyLM=function(data,formula,x,weights,predictors,response,intercept,intercept_val,intercept_name,n,terms){
  flist <- lapply(as.list(x[,predictors]),"factor")
  times=sapply(flist, function(x){return(length(unique(x)))})

  coefs <- c(intercept_val,rep(NA,sum(times-1)))
  names(coefs) <- c(intercept_name,names(unlist(sapply(flist, function(x){return(unique(x)[-1])}))))

  emptyLM <- list(coefficients=coefs,residuals=rep(NA, n),fitted.values=rep(NA, n),effects=rep(NA, n),rank=sum(times-1)+as.numeric(intercept),weights=rep(1, n), qr=list(qr=matrix(NA,max(nrow(x),length(coefs)), length(coefs)), qraux=rep(NA, length(coefs)),pivot=rep(NA, length(coefs)),tol=NA,rank=length(coefs)),df.residual=n-(sum(times-1)+as.numeric(intercept)),call=call("lm(formula = formula, data = data, weights = weights)"),terms=terms,contrasts=NA,xlevels=NA,offset=NA,y=NA,x=NA,model=data.frame(x[c(response,predictors)], "(weights)"=rep(1, n), check.names=FALSE),na.action=NA)
  class(emptyLM) <- "lm"

  return(emptyLM)
}


#This function creates an empty lmerMod object
.emptylmerMod=function(data,formula,x,y,predictors,response,intercept,intercept_val,intercept_name,n,beta){
  cnms <- rep(list("(Intercept)"),length(predictors))
  names(cnms) <- predictors
  flist <- lapply(as.list(x[,predictors]),"factor")
  attr(flist, "assign") <- seq(1,length(predictors))
  nrandom <- sum(sapply(predictors, function(y){return(length(unique(do.call('$',list(x,y)))))}))
  y <- do.call('$', list(x,response))
  X <- matrix(intercept_val,nrow=length(y), dimnames=list(NULL,intercept_name))
  Lind <- rep(1:length(predictors),times=sapply(flist, function(y){return(length(unique(y)))}))
  emptylmerMod <- new("lmerMod", resp=new("lmerResp",Ptr="<externalptr>", mu=rep(NA,n),offset=rep(0,n),sqrtXwt=rep(1,n),sqrtrwt=rep(1,n),weights=rep(1,n),y=y),Gp=as.integer(0),call=call(paste0("lme4::lmer(formula = ","formula",", data = data, weights = weights)")),frame=x,flist=flist,cnms=cnms,lower=as.numeric(NA),theta=as.numeric(rep(NA, length(cnms))),beta=beta,u=rep(as.numeric(NA),nrandom),devcomp=list(dims=c(N=n,n=length(y),p=as.numeric(intercept),nmp=0,nth=1,q=nrandom,nAGQ=NA,compDev=TRUE,useSc=TRUE,reTrms=1,REML=1,GLMM=FALSE,NLMM=FALSE),cmp=c(ldL2=NA,ldRX2=NA,wrss=NA,ussq=NA,pwrss=NA,drsum=NA,REML=1,dev=NA,sigmaML=NA,sigmaREML=NA,tolPwrss=NA)),pp=new("merPredD",X=X,Zt=do.call("rbind", sapply(x[,predictors], Matrix::fac2sparse)),Lambdat=Matrix::Matrix(diag(nrandom),sparse=TRUE),Lind=Lind,theta=as.numeric(rep(NA, length(cnms))),n=1),optinfo=list(optimizer="bobyqa",control=list(iprint=0),derivs=list(gradient=NA,Hessian=NA),conv=list(opt=0,lme4=list()),feval=NA,warnings=list(),val=NA))
  return(emptylmerMod)
}


#This function adjusts the names for the factors in a data object
.adjustNames=function(data, predictors){
  data_adj <- lapply(data, function(x){

    for(i in 1:length(x[,predictors]))
    {
      if(is.factor(x[,predictors][[i]]) || is.character(x[,predictors][[i]]))
      {
        x[,predictors][[i]] <- factor(paste0(colnames(x[,predictors])[i],x[,predictors][[i]]), levels=unique(paste0(colnames(x[,predictors])[i],x[,predictors][[i]])))
      }
    }

    return(x)
  })
  return(data_adj)
}


#Function for ridge regression
.ridgeEstM3=function(parsedFormula, k, tolPwrss = 1e-10, verbose=FALSE, ...)
{

  #If ridge regression with Huber weights:
  if(!is.null(parsedFormula$fr$`(weights)`[1]) && parsedFormula$fr$`(weights)`[1]=="Huber"){ridgeFit <- .huberRidgeFit(parsedFormula, k, tolPwrss = tolPwrss, verbose=verbose, ...)

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

  optimizerOutput <- tryCatch(lme4::optimizeLmer(devianceFunction), error=function(e){return(NULL)})

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
.huberRidgeFit=function(parsedFormula, k, tolPwrss = 1e-10, verbose=FALSE, ...){

  #Remove "Huber" from the `(weights)` column
  parsedFormula$fr$`(weights)` <- NULL

  pwrssnew <- Inf
  weights <- NULL
  start <- NULL ###NEW!!!

  devianceFunction <- tryCatch(do.call(lme4::mkLmerDevfun, parsedFormula), error=function(e){
    return(NULL)})

  repeat{
    pwrssold <- pwrssnew

    #!!!TryCatch aanpasssen enzo!!!
    parsedFormula$start <- start ###NEW!!!

    p <- ncol(parsedFormula$X)
    environment(devianceFunction)$resp <- lme4:::mkRespMod(parsedFormula$fr, REML = p)

    #devianceFunction <- lme4:::mkdevfun(environment(devianceFunction), 0L)

    #environment(devianceFunction)$pp$updateDecomp()
    #environment(devianceFunction)$pp$updateL()
    #environment(devianceFunction)$pp$updateLamtUt()
    #environment(devianceFunction)$pp$updateRes(wtres=resid)
    #environment(devianceFunction)$pp$updateXwts(wts=MASS::psi.huber(resid/sigma,k=k))

    #devianceFunction(environment(devianceFunction)$pp$theta) # one evaluation to ensure all values are set

    optimizerOutput <- tryCatch(lme4::optimizeLmer(devianceFunction), error=function(e){return(NULL)})

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

    if(pwrssold-pwrssnew <= tolPwrss){
      break
    }

    n <- nrow(pp$V)
    p <- ncol(pp$V)

    sigmaML <- pwrssnew/n
    #sigma
    sigma <- sqrt(unname(sigmaML*(n/(n-p))))
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
    parsedFormula$fr$`(weights)` <- MASS::psi.huber(resid/sigma,k=k)

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
.huberLmFit=function(formula,data,k=1.345,tolWrss = 1e-10,verbose=FALSE,...)
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
    .weight.=MASS::psi.huber(summary$residuals/summary$sigma,k=k)
    #Weighted residual sum of squares
    wrssnew <- sum(summary$residuals^2*.weight.) #136.4315

    if(wrssold-wrssnew <= tolWrss){
      break
    }
  }
  if (verbose) print(lmFit)

  return(lmFit)
}


