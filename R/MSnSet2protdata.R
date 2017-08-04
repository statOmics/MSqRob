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
MSnSet2protdata <- function(MSnSet, accession, annotations=NULL, quant_name="quant_value", printProgress=FALSE, shiny=FALSE, message=NULL){

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

  fData <- Biobase::fData(MSnSet)
  pData <- Biobase::pData(MSnSet)
  row.names(pData) <- NULL
  exprs <- Biobase::exprs(MSnSet)

  accessions <- unique(fData[,accession])
  datalist <- replicate(length(accessions), data.frame())
  annotation_matrix <- matrix(nrow=length(accessions), ncol=length(annotations), dimnames=list(accessions,annotations))

  for(i in 1:length(accessions)){

    updateProgress(progress=progress, detail=paste0("Converting protein ",i," of ",length(accessions),"."), n=length(accessions), shiny=shiny, print=isTRUE(printProgress))

    sel <- fData[,accession] == accessions[i]
    intensities <- exprs[sel,,drop=FALSE]
    properties <- fData[sel,,drop=FALSE]
    #If for the same accession in an annation column, there would be multiple different values, just paste them together.
    annotation_matrix[i,] <- matrix(as.character(unlist(lapply(properties[,annotations, drop=FALSE], function(x){paste0(unique(x), collapse="")}))),nrow=1)
    properties <- properties[,-which(colnames(properties) %in% c(colnames(properties[,accession,drop=FALSE]),colnames(properties[,annotations,drop=FALSE]))), drop=FALSE]

    quant_value <- as.vector(intensities)

    framelist <- replicate(ncol(intensities), droplevels(properties), simplify = FALSE)

    frame <- plyr::rbind.fill(framelist)

    frame2 <- cbind(data.frame(quant_value=quant_value,frame), pData[rep(row.names(pData), each=nrow(intensities)),,drop=FALSE])
    #Name the newly create intensity column.
    names(frame2)[1] <- quant_name

    #Removing NA
    frame3 <- frame2[complete.cases(frame2[,1]),,drop=FALSE]

    datalist[[i]] <- frame3
  }

  protdata <- new("protdata", accession=accessions, data=datalist, annotation=annotation_matrix)
  return(protdata)
}
