#' Protein Linear Model Object - class
#'
#' An S4 class for storing the model fits of feature- or peptide-level data of a set of LC-MS proteomics experiments, grouped by protein.
#' Objects of this class are supplied as an argument to the \code{\link{test.protLMcontrast}} function.
#' @slot protein A factor vector with length \code{k} equal to the total number of proteins or protein groups in the data.
#' @slot model A list of length \code{k} containing the corresponding fitted models.
#' @include protdata.R

setClass("protLM",
         slots = list(accession="factor", model="list", annotation="matrix"),
         prototype = list(accession=factor(), model=list(), annotation=matrix(nrow=0,ncol=0)))

#Validity check

.protLM.valid <- function(object){
  len1 <- length(object@accession)
  len2 <- length(object@model)
  len3 <- nrow(object@annotation)
  if(length(unique(c(len1,len2,len3)))!=1){
    return("Length mismatch")
  } else if(!all(sapply(object@model,class) %in% c("lmerMod","lm"))) {
    return("Model slot can only contain fitted model objects of type \"lmerMod\" or \"lm\".")
  } else if(!is.matrix(object@annotation)) {
    return("Annotation slot should be a matrix.")
  } else {
    return(TRUE)
  }
}

setValidity("protLM", .protLM.valid)

#Uitprinten

.ProtLM.print <- function(x,n=6){

  if(length(x@accession)<n){n <- length(x@accession)}

  over <- length(x@accession)-n
  annotEW <- x@annotation
  if(ncol(annotEW)==0){annotEW <- matrix(rep("",nrow(annotEW)))}

  for(i in 1:n)
  {
    cat("@accession",
        "\n\n")
    print(as.character(x@accession[i]))

    cat("\n",
        "@annotation",
        "\n\n")
    print(annotEW[i,])

    cat("\n",
        "@model",
        "\n\n")
    print(x@model[[i]])
    cat("\n")
  }
  if(over>0){
    cat(paste0(over," more elements..."),
        "\n")
  }
}

print.protLM <- function(x, ...){
  jj <- .ProtLM.print(x, ...)
  #print(jj)
  #Verhindert dat de invisible NULL van cat geprint wordt
  jj
  return(invisible(jj))
}

setMethod("show", "protLM", function(object){print.protLM(object)})

#Subsetting protLM objects

# setGeneric (
#   name= "getAccessions",
#   def=function(object){standardGeneric("getAccessions")}
# )

#' Extract the accession slot of a protLM object
#'
#' \code{getAccessions} extracts the \code{accession} slot from a \code{\link[=protLM-class]{protLM}} object.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the \code{accession} slot needs to be extracted.
#' @return A factor vector containing the accessions present in the \code{\link[=protLM-class]{protLM}} object.
#' @export
setMethod("getAccessions", "protLM",
          function(object){
            return(object@accession)})

setGeneric (
  name= "getModels",
  def=function(object,...){standardGeneric("getModels")}
)

#' Extract the model slot of a protLM object
#'
#' \code{getModels} extracts the \code{model} slot from a \code{\link[=protLM-class]{protLM}} object.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the \code{model} slot needs to be extracted.
#' @param simplify Logical. Should the list of models be simplified to a model frame if the length of this list equals 1?
#' @return The list of models corresponding to the accessions present in the \code{\link[=protLM-class]{protLM}} object.
#' @export
setMethod("getModels", "protLM",
          function(object,simplify=TRUE){

            if(simplify & length(object@accession)==1){
              return(object@model[[1]])
            } else{
              return(object@model)
            }

          })

# setGeneric (
#   name= "getAnnotations",
#   def=function(object,...){standardGeneric("getAnnotations")}
# )

#' Extract the annotation slot of a protLM object
#'
#' \code{getAnnotations} extracts the \code{annotation} slot from a \code{\link[=protLM-class]{protLM}} object.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the \code{annotation} slot needs to be extracted.
#' @return A character matrix containing the annotations present in the \code{\link[=protLM-class]{protLM}} object. Rows correspond to accessions, while columns correspond to different annotations.
#' @export
setMethod("getAnnotations", "protLM",
          function(object){

            annotations <- object@annotation
            rownames(annotations) <- rownames(object@accession)
            return(annotations)

          })


#' Select elements from a protLM object
#'
#' Select certain elements from a \code{\link[=protLM-class]{protLM}} object based on one or more accession and/or numeric indices.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which certain elements will be selected.
#' @param i Either a character or a factor vector containing accessions present in the \code{accessions} slot of the \code{\link[=protLM-class]{protLM}} object or a numeric vector containing indices pointing to certain elements in the the \code{\link[=protLM-class]{protLM}} object or a combination of both (deprecated).
#' @return A reduced \code{\link[=protLM-class]{protLM}} object containing only the elements specified in \code{i}.
#' @export
setMethod("[", "protLM",
          function(x, i, j,  drop){
            if(!missing(j)){
              stop("No second argument to extractor function allowed!")
            }
            i <- index2Numeric(i,x@accession)
            #initialize(x, accession=x@accession[i], intensities=x@intensities[i,,drop=FALSE], properties=x@properties[i,])
            new("protLM", accession=x@accession[i], model=x@model[i,drop = FALSE], annotation=x@annotation[i,,drop = FALSE])
          })

#' Length of a protLM object
#'
#' The \code{length} method counts the number of elements in a \code{\link[=protLM-class]{protLM}} object.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the number of elements needs to be determined.
#' @return An integer value indicating the number of elements in a \code{\link[=protLM-class]{protLM}} object.
#' @export
setMethod("length", "protLM",
          function(x){
            return(length(x@accession)[1])})

#Set accessions, data, annotations

# setGeneric (
#   name= "setAccessions",
#   def=function(object,accessions){standardGeneric("setAccessions")}
# )

#' Change the accession slot of a protLM object
#'
#' \code{setAccessions} changes the \code{accession} slot from a \code{\link[=protLM-class]{protLM}} object to a prespecified factor vector.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the \code{accession} slot needs to be modified.
#' @param accessions A factor or character vector of accessions which will be inserted in the \code{accession} slot of the \code{\link[=protLM-class]{protLM}} object. Note that the length of this vector should be equal to the number of models in the \code{model} slot and the number of rows in the \code{annotation} slot.
#' @return A \code{\link[=protLM-class]{protLM}} object of which the \code{accession} slot is replaced by the \code{accessions} factor vector.
#' @export
setMethod("setAccessions", "protLM",
          function(object, accessions){
            return(new("protLM", accession=factor(accessions), model=object@model, annotation=object@annotation))          })

setGeneric (
  name= "setModels",
  def=function(object, model){standardGeneric("setModels")}
)

#' Change the model slot of a protLM object
#'
#' \code{setModels} changes the \code{model} slot from a \code{\link[=protLM-class]{protLM}} object to a prespecified list of models.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the \code{model} slot needs to be modified.
#' @param model A list of models which will be inserted in the \code{model} slot of the \code{\link[=protLM-class]{protLM}} object. Note that the length of this list should be equal to the number of elements in the \code{accession} slot and the number of rows in the \code{annotation} slot.
#' @return A \code{\link[=protLM-class]{protLM}} object of which the \code{model} slot is replaced by the list of models \code{models}.
#' @export
setMethod("setModels", "protLM",
          function(object, model){
            return(new("protLM", accession=object@accession, model=models, annotation=object@annotation))          })

# setGeneric (
#   name= "setAnnotations",
#   def=function(object, annotations){standardGeneric("setAnnotations")}
# )

#' Change the annotation slot of a protLM object
#'
#' \code{setAnnotations} changes the \code{annotation} slot from a \code{\link[=protLM-class]{protLM}} object to a prespecified matrix.
#' @param object A \code{\link[=protLM-class]{protLM}} object of which the \code{annotation} slot needs to be modified.
#' @param annotations A matrix which will be inserted in the \code{annotation} slot of the \code{\link[=protLM-class]{protLM}} object. Note that the number of rows of this matrix should be equal to the number of elements in both the \code{accession} and the \code{model} slot.
#' @return A \code{\link[=protLM-class]{protLM}} object of which the \code{annotation} slot is replaced by the matrix \code{annotations}.
#' @export
setMethod("setAnnotations", "protLM",
          function(object, annotations){
            return(new("protLM", accession=object@accession, model=object@model, annotation=annotations))          })



#Select specific accessions

# setGeneric (
#   name= "selectAccessions",
#   def=function(object,i,...){standardGeneric("selectAccessions")}
# )

#' @export
setMethod("selectAccessions", "protLM",
          function(object,i,keep=TRUE){
            if(keep==TRUE){
              #Makes use of previous method "["!
              return(object[i])
            } else{
              i <- index2Numeric(i,object@accession)
              return(new("protLM", accession=object@accession[-i], model=object@model[-i,], annotation=object@annotation[-i,]))
            }
          })


