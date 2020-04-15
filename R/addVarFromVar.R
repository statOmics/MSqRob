#' Add a variable to a protdata object based on an existing column
#'
#' @description Adds a categorical variable to the data slot of a \code{\link[=protdata-class]{protdata}} object based on an exisiting variable.
#' @param protdata A \code{\link[=protdata-class]{protdata}} object of which the data slot needs to be expanded.
#' @param basecol The name of the column in the existing data on which you base your new column.
#' @param name The name of the new variable.
#' @param vector A vector that contains the values of the new variable with names corresponding to all possible values of the \code{basecol}.
#' @return A \code{\link[=protdata-class]{protdata}} object with an added column in its data slot containing the new variable.
#' @examples #Here we show how the columns conc, rep and instrlab were added to the object proteinsCPTAC
#' data(proteinsCPTAC, package="MSqRob")
#' #Add a grouping factor for the spike-in condition (i.e. 6A, 6B, 6C, 6D and 6E)
#' conc <- factor(substr(levels(getData(proteinsCPTAC[1])$run), 11,12))
#' names(conc) <- levels(getData(proteinsCPTAC[1])$run)
#' proteinsCPTAC <- addVarFromVar(proteinsCPTAC,"run","conc",conc)
#'
#' #Add a grouping factor for each of the 9 technical repeats (i.e. instrument x run)
#' rep <- factor(substr(levels(getData(proteinsCPTAC[1])$run), 14,14))
#' proteinsCPTAC <- addVarFromVar(proteinsCPTAC,"run","rep",rep)
#'
#' #Add a grouping factor for the instrument effect (i.e. LTQ-Orbitrap at site 86, LTQ-Orbitrap O at site 65 and LTQ-Orbitrap W at site 56)
#' instrlab <- factor((as.numeric(rep)-1)%/%3+1)
#' proteinsCPTAC <- addVarFromVar(proteinsCPTAC,"run","instrlab",instrlab)
#' @include protdata.R
#' @export
addVarFromVar <- function(protdata, basecol, name, vector){

   protdatanew <- lapply(protdata@data, function(x){

     #Make a data frame equal to the existing one with an extra column
     #that places the values in "vector" next to their corresponding values in "basecol".
     x <- data.frame(x,vector[as.character(do.call('$',list(x,basecol)))], stringsAsFactors = TRUE)

     #Give the chosen name to the added column.
     colnames(x)[ncol(x)] <- name

     return(x)

     })

   #Put the new list of data frames in the data slot of the new protdata object.
   return(new("protdata", accession=protdata@accession, data=protdatanew, annotation=protdata@annotation))
}
