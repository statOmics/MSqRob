#' Add a combination of two or more factor variables to a protData object
#'
#' @description Adds a categorical variable to the data slot of a \code{\link[=protdata-class]{protdata}} that is a combination of two or more exisiting factor variables.
#' @param protdata A \code{\link[=protdata-class]{protdata}} object of which the data slot needs to be expanded.
#' @param basecols The name of the columns in the existing data on which you base your new column.
#' @param name The name of the new variable. If no name is provided (default: \code{name=NULL}), the names of the \code{basecols} will be concatenated separated by an "_".
#' @return A \code{\link[=protdata-class]{protdata}} object with an added column in its data slot containing the new variable.
#' @examples #.........
#' @include protdata.R
#' @export
addFactorInteractions <- function(protdata, basecols, name=NULL){

  protdatanew <- lapply(protdata@data, function(x){

    if(all(vapply(x[,basecols], "is.numeric",TRUE))){
      x <- data.frame(x,apply(x[,basecols], 1, prod))
    } else{
      x <- data.frame(x,apply(x[,basecols], 1, paste, collapse="_"))
    }

    if(is.null(name)){colnames(x)[ncol(x)] <- paste(colnames(x[,basecols]), collapse="_")
    } else{
    #Give the chosen name to the added column.
    colnames(x)[ncol(x)] <- name
    }

    return(x)

  })

  #Put the new list of data frames in the data slot of the new protdata object.
  return(new("protdata", accession=protdata@accession, data=protdatanew, annotation=protdata@annotation))
}
