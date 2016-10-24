#' Protein Data Object - class
#'
#' An S4 class for storing the feature- or peptide-level data of a set of LC-MS proteomics experiments, grouped by accession (typically protein accession).
#' Objects of this class are supplied as an argument to the \code{\link{fit.model}} function.
#' @slot accession A factor vector with length \code{k} equal to the total number of accessions (mostly proteins or protein groups) in the data.
#' @slot data A list of length \code{k} containing data frames with preprocessed intensity values in a single column and other covariates of interest in the other columns.

setClass("protdata",
         slots = list(accession="factor", data="list", annotation="matrix"),
         prototype = list(accession=factor(), data=list(), annotation=matrix(nrow=0,ncol=0)))

#Validity check

.protdata.valid <- function(object){
  len1 <- length(object@accession)
  len2 <- length(object@data)
  len3 <- nrow(object@annotation)
  if(length(unique(c(len1,len2,len3)))!=1){
    return("Length mismatch")
  } else if(!all(vapply(object@data,is.data.frame, FUN.VALUE = logical(1)))) {
    return("Data slot can only contain data frames")
  } else if(length(unique(lapply(object@data,colnames)))!=1) {
    return("Data frames in data slot must have the same number of columns and identical colnames")
  } else if(!is.matrix(object@annotation)) {
    return("Annotation slot should be a matrix.")
  } else {
    return(TRUE)
  }
}


setValidity("protdata", .protdata.valid)


#Printing protdata (3 stages)

.Protdata.print <- function(x,n=6,n2=6){

  n2init <- n2
  if(length(x@accession)<n){n <- length(x@accession)}

  over <- length(x@accession)-n
  annotEW <- x@annotation
  if(ncol(annotEW)==0){annotEW <- matrix(rep("",nrow(annotEW)))}

  for(i in 1:n)
  {
    dataEW <- x@data[[i]]
    if(nrow(dataEW)<n2){n2 <- nrow(dataEW)}
    over2 <- nrow(dataEW)-n2
    cat("@accession",
      "\n\n")
  print(as.character(x@accession[i]))

  cat("\n",
      "@annotation",
      "\n\n")
  print(annotEW[i,])

  cat("\n",
      "@data",
      "\n\n")
  print(utils::head(dataEW,n2))

  if(over2>0){
    cat("\n",
        paste0(over2," more rows..."),
        "\n\n")
  }
  n2 <- n2init
  }
  if(over>0){
    cat(paste0(over," more elements..."),
        "\n")
  }
}

print.protdata <- function(x, ...){
  jj <- .Protdata.print(x, ...)
  #print(jj)
  #Verhindert dat de invisible NULL van cat geprint wordt
  jj
  return(invisible(jj))
}

setMethod("show", "protdata", function(object){print.protdata(object)})


#Subsetting protdata objects

setGeneric (
  name= "getAccessions",
  def=function(object){standardGeneric("getAccessions")}
)

#' Extract the accession slot of a protdata object
#'
#' \code{getAccessions} extracts the \code{accession} slot from a \code{\link[=protdata-class]{protdata}} object.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which the \code{accession} slot needs to be extracted.
#' @return A factor vector containing the accessions present in the \code{\link[=protdata-class]{protdata}} object.
#' @export
setMethod("getAccessions", "protdata",
          function(object){
            return(object@accession)})

setGeneric (
  name= "getData",
  def=function(object,...){standardGeneric("getData")}
)

#' Extract the data slot of a protdata object
#'
#' \code{getData} extracts the \code{data} slot from a \code{\link[=protdata-class]{protdata}} object.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which the \code{data} slot needs to be extracted.
#' @param simplify Logical. Should the list of data frames be simplified to a data frame if the length of this list equals 1?
#' @return The list of data frames corresponding to the data present in the \code{\link[=protdata-class]{protdata}} object.
#' @export
setMethod("getData", "protdata",
          function(object,simplify=TRUE){

            if(simplify & length(object@accession)==1){
              return(object@data[[1]])
            } else{
              return(object@data)
            }

          })


setGeneric (
  name= "getAnnotations",
  def=function(object,...){standardGeneric("getAnnotations")}
)

#' Extract the annotation slot of a protdata object
#'
#' \code{getAnnotations} extracts the \code{annotation} slot from a \code{\link[=protdata-class]{protdata}} object.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which the \code{annotation} slot needs to be extracted.
#' @return A character matrix containing the annotations present in the \code{\link[=protdata-class]{protdata}} object. Rows correspond to accessions, while columns correspond to different annotations.
#' @export
setMethod("getAnnotations", "protdata",
          function(object){

            annotations <- object@annotation
            rownames(annotations) <- rownames(object@accession)
              return(annotations)

          })


#' Select elements from a protdata object
#'
#' Select certain elements from a \code{\link[=protdata-class]{protdata}} object based on one or more accession and/or numeric indices.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which certain elements will be selected.
#' @param i Either a character or a factor vector containing accessions present in the \code{accessions} slot of the \code{\link[=protdata-class]{protdata}} object or a numeric vector containing indices pointing to certain elements in the the \code{\link[=protdata-class]{protdata}} object or a combination of both (deprecated).
#' @return A reduced \code{\link[=protdata-class]{protdata}} object containing only the elements specified in \code{i}.
#' @export
setMethod("[", "protdata",
          function(x, i, j,  drop){
            if(!missing(j)){
              stop("No second argument to extractor function allowed!")
            }
            i <- index2Numeric(i,x@accession)
            #initialize(x, accessions=x@accessions[i], intensities=x@intensities[i,,drop=FALSE], properties=x@properties[i,])
            new("protdata", accession=x@accession[i], data=x@data[i,drop = FALSE], annotation=x@annotation[i,,drop = FALSE])
          })


#' Length of a protdata object
#'
#' The \code{length} method counts the number of elements in a \code{\link[=protdata-class]{protdata}} object.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which the number of elements needs to be determined.
#' @return An integer value indicating the number of elements in a \code{\link[=protdata-class]{protdata}} object.
#' @export
setMethod("length", "protdata",
          function(x){
            return(length(x@accession))})



#Set accessions, data, annotations

setGeneric (
  name= "setAccessions",
  def=function(object,accessions){standardGeneric("setAccessions")}
)

#' Change the accession slot of a protdata object
#'
#' \code{setAccessions} changes the \code{accession} slot from a \code{\link[=protdata-class]{protdata}} object to a prespecified factor vector.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which the \code{accession} slot needs to be modified.
#' @param accessions A factor or character vector of accessions which will be inserted in the \code{accession} slot of the \code{\link[=protdata-class]{protdata}} object. Note that the length of this vector should be equal to the number of elements in the \code{data} slot and the number of rows in the \code{annotation} slot.
#' @return A \code{\link[=protdata-class]{protdata}} object of which the \code{accession} slot is replaced by the \code{accessions} factor vector.
#' @export
setMethod("setAccessions", "protdata",
          function(object, accessions){
            return(new("protdata", accession=factor(accessions), data=object@data, annotation=object@annotation))          })

setGeneric (
  name= "setData",
  def=function(object, data){standardGeneric("setData")}
)

#' Change the data slot of a protdata object
#'
#' \code{setData} changes the \code{data} slot from a \code{\link[=protdata-class]{protdata}} object to a prespecified list of data frames.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which the \code{data} slot needs to be modified.
#' @param data A list of data frames which will be inserted in the \code{data} slot of the \code{\link[=protdata-class]{protdata}} object. Note that the length of this list should be equal to the number of elements in the \code{accession} slot and the number of rows in the \code{annotation} slot.
#' @return A \code{\link[=protdata-class]{protdata}} object of which the \code{data} slot is replaced by the list of data frames \code{data}.
#' @export
setMethod("setData", "protdata",
          function(object, data){
            return(new("protdata", accession=object@accession, data=data, annotation=object@annotation))          })

setGeneric (
  name= "setAnnotations",
  def=function(object, annotations){standardGeneric("setAnnotations")}
)

#' Change the annotation slot of a protdata object
#'
#' \code{setAnnotations} changes the \code{annotation} slot from a \code{\link[=protdata-class]{protdata}} object to a prespecified matrix.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which the \code{annotation} slot needs to be modified.
#' @param annotations A matrix which will be inserted in the \code{annotation} slot of the \code{\link[=protdata-class]{protdata}} object. Note that the number of rows of this matrix should be equal to the number of elements in both the \code{accession} and the \code{data} slot.
#' @return A \code{\link[=protdata-class]{protdata}} object of which the \code{annotation} slot is replaced by the matrix \code{annotations}.
#' @export
setMethod("setAnnotations", "protdata",
          function(object, annotations){
            return(new("protdata", accession=object@accession, data=object@data, annotation=annotations))          })


#Select specific accessions

setGeneric (
  name= "selectAccessions",
  def=function(object,i,...){standardGeneric("selectAccessions")}
)

#' Select accessions from a protdata object
#'
#' Select certain elements from a \code{\link[=protdata-class]{protdata}} object based on one or more accession and/or numeric indices. If \code{keep==FALSE}, these indices are omitted instead of selected. Note that for selecting, one could also use the "[" method.
#' @param object A \code{\link[=protdata-class]{protdata}} object of which certain elements will be selected or removed.
#' @param i Either a character or a factor vector containing accessions present in the \code{accessions} slot of the \code{\link[=protdata-class]{protdata}} object or a numeric vector containing indices pointing to certain elements in the the \code{\link[=protdata-class]{protdata}} object or a combination of both (deprecated).
#' @param keep A logical indicating whether the selected elements should be retained (\code{keep==TRUE}) or removed (\code{keep==FALSE}). Defaults to \code{TRUE}.
#' @return A reduced \code{\link[=protdata-class]{protdata}} object containing only the elements specified in \code{i} (if \code{keep==TRUE}) or all elements except those specified in \code{i} (if \code{keep==FALSE}).
#' @export
setMethod("selectAccessions", "protdata",
          function(object,i,keep=TRUE){
            if(isTRUE(keep)){
              #Makes use of previous method "["!
              return(object[i])
            } else{
                i <- index2Numeric(i,object@accession)
              return(new("protdata", accession=object@accession[-i], data=object@data[-i], annotation=object@annotation[-i,]))
            }
          })


#' Convert index to only numeric values
#'
#' This function converts a vector of indices, which are a potential mixture of numeric values and accessions (characters or even factors)
#' to a vector which retains the numeric indices, but converts accessions to their corresponding index in a provided \code{accessions} vector.
#' @param i A vector of indices (either numeric or a character or factor vector corresponding to elements in the \code{accessions} vector).
#' @param accessions A character or factor vector containing accessions.
#' @return A numeric vector containing either the original numeric values or the numeric values corresponding to the position of each element of \code{i} in the \code{accessions} vector.
index2Numeric <- function(i, accessions){
  i[i %in% accessions] <- which(accessions %in% i)
  i <- as.numeric(i)
  if(anyDuplicated(i)){warning("Some accessions are selected multiple times!")}
  if(anyNA(i)){stop("You tried to select some accessions that were not present in the \"Accessions\" slot or you tried to select out of bounds indices.")}
  return(i)
}


#'
#' setGeneric (
#'   name= "selectData",
#'   def=function(object,i,...){standardGeneric("selectData")}
#' )
#'
#' #' @export
#' setMethod("selectData", "protdata",
#'           function(object,i,keep=TRUE){
#'             if(is.character(i) | is.factor(i)){
#'               i <- which(object@accession %in% i)
#'             }
#'             if(keep==TRUE){
#'               return(new("protdata", accession=object@accession[i], data=object@data[i,drop=FALSE]))
#'             } else{
#'               return(new("protdata", accession=object@accession[i], data=object@data[-i,drop=FALSE]))
#'             }
#'           })
#'


#Add and remove data columns

