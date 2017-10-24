#' MSnSet to protein data
#'
#' @description Converts an \code{\link[=MSnSet-class]{MSnSet}} object to a a \code{\link[=protdata-class]{protdata}} object, which can be used for further analysis.
#' @param MSnSet An \code{\link[=MSnSet-class]{MSnSet}} object that needs to be converted to a \code{\link[=protdata-class]{protdata}} object.
#' @param accession A character string or numeric index indicating the column in the fData slot of object  \code{\link[=MSnSet-class]{MSnSet}} that contains the identifiers by which the data should be grouped in the resulting \code{\link[=protdata-class]{protdata}} object (typically the protein identifier).
#' @param annotations A vector of character strings or numeric indices indicating the columns in the fData slot of object  \code{\link[=MSnSet-class]{MSnSet}} that contain additional information on the accessions (typically protein names, gene names, gene ontologies,...) that should be retained. In case multiple values of the same annotation column would exist for a unique accession, these values are pasted together. Defaults to \code{NULL}, in which case no annotations will be added.
#' @param quant_name A character string indicating the name to be given to the column that will contain the quantitative values of interest (mostly peptide intensities or peptide areas under the curve). Defaults to \code{"quant_value"}.
#' @param printProgress A logical indicating whether the R should print a message before converting each accession. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @return A \code{\link[=protdata-class]{protdata}} object.
#' @examples #This example will convert MSnSet object peptides into a protdata object proteins.
#' library(MSnbase)
#' colInt <- grepEcols(system.file("extdata/CPTAC", "peptides.txt", package = "MSqRob"), pattern="Intensity.", split = "\t")
#' #Import the data as an MSnSet object
#' peptides <- readMSnSet2(system.file("extdata/CPTAC", "peptides.txt", package = "MSqRob"), ecol = colInt, sep = "\t")
#' fData(peptides) <- fData(peptides)[,c("Proteins","Sequence","PEP")]
#' #Do this only for the first 50 peptide to save execution time
#' proteins <- MSnSet2protdata(peptides[1:50], "Proteins")
#' @include updateProgress.R
#' @include protdata.R
#' @include preprocess_MaxQuant.R
#' @export
MSnSet2protdata <- function(MSnSet, accession=NULL, annotations=NULL, quant_name="quant_value", printProgress=FALSE, shiny=FALSE, message=NULL){

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  #If MSnSet is no MSnSet object, throw error
  if(class(MSnSet)!="MSnSet"){stop("MSnSet should be an object of class \"MSnSet\".")}
  #Extra check to see if MSnSet object is valid
  if(!all(rownames(Biobase::pData(MSnSet))==colnames(Biobase::exprs(MSnSet)))){stop("Possible problem with annotation! Rownames of pData slot are not equal to colnames of exprs slot!")}
  #An accession must be provided
  if(is.null(accession)){stop("Please provide an accession column.")}

  fData <- Biobase::fData(MSnSet)
  pData <- Biobase::pData(MSnSet)
  row.names(pData) <- NULL
  exprs <- Biobase::exprs(MSnSet)

  accessions <- unique(fData[,accession])

  sel <- lapply(accessions, function(x){return(which(fData[,accession] == x))})

  datalist <- lapply(sel, function(x){

    properties <- fData[x,-which(colnames(fData) %in% c(colnames(fData[,accession,drop=FALSE]),colnames(fData[,annotations,drop=FALSE]))), drop=FALSE]

    return(
      cbind(
        data.frame(c(exprs[x,,drop=FALSE])),
        lapply(properties, function(z){
          return(droplevels(rep(z, ncol(exprs))))
        }),
        pData[rep(row.names(pData), each=length(x)),,drop=FALSE]
      )
    )})

  datalist <- lapply(datalist, function(x){
    #Give name
    names(x)[1] <- quant_name

    #Removing NA
    x <- x[complete.cases(x[,1]),,drop=FALSE]
    return(x)
  })

  properties <- fData[,annotations,drop=FALSE]

  annotation_matrix <- matrix(nrow=length(accessions), ncol=length(annotations), dimnames=list(accessions,annotations))

  for(i in 1:nrow(annotation_matrix)){
    annotation_matrix[i,] <- vapply(properties[sel[[1]],, drop=FALSE], function(y){paste0(unique(y), collapse="")},"; ")
  }

  protdata <- new("protdata", accession=accessions, data=datalist, annotation=annotation_matrix, pData=pData)
  return(protdata)
}
