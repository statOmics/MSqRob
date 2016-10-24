#' MSnSet to protein data
#'
#' @description Converts an \code{\link[=MSnSet-class]{MSnSet}} object to a a \code{\link[=protdata-class]{protdata}} object, which can be used for further analysis.
#' @param MSnSet An \code{\link[=MSnSet-class]{MSnSet}} object that needs to be converted to a \code{\link[=protdata-class]{protdata}} object.
#' @param accession A character string or numeric index indicating the column in the fData slot of object  \code{\link[=MSnSet-class]{MSnSet}} that contains the identifiers by which the data should be grouped in the resulting \code{\link[=protdata-class]{protdata}} object (typically the protein identifier).
#' @param annotations A vector of character strings or numeric indices indicating the columns in the fData slot of object  \code{\link[=MSnSet-class]{MSnSet}} that contain additional information on the accessions (typically protein names, gene names, gene ontologies,...) that should be retained. In case multiple values of the same annotation column would exist for a unique accession, these values are pasted together. Defaults to \code{NULL}, in which case no annotations will be added.
#' @return A \code{\link[=protdata-class]{protdata}} object.
#' @examples #This example will convert MSnSet object peptides into a protdata object proteins.
#' library(MSnbase)
#' colInt <- grepEcols(system.file("extdata/CPTAC", "peptides.txt", package = "MSqRob"), pattern="Intensity.", split = "\t")
#' #Import the data as an MSnSet object
#' peptides <- readMSnSet2(system.file("extdata/CPTAC", "peptides.txt", package = "MSqRob"), ecol = colInt, sep = "\t")
#' fData(peptides) <- fData(peptides)[,c("Proteins","Sequence","PEP")]
#' #Do this only for the first 50 peptide to save execution time
#' proteins <- MSnSet2protdata(peptides[1:50], "Proteins")
#' @include protdata.R
#' @include preprocess_MaxQuant.R
#' @export
MSnSet2protdata <- function(MSnSet, accession, annotations=NULL){

  #If MSnSet is no MSnSet object, throw error
  if(class(MSnSet)!="MSnSet"){stop("MSnSet should be an object of class \"MSnSet\".")}

  fData <- Biobase::fData(MSnSet)
  pData <- Biobase::pData(MSnSet)
  row.names(pData) <- NULL
  exprs <- Biobase::exprs(MSnSet)

  accessions <- unique(fData[,accession])
  datalist <- replicate(length(accessions), data.frame())
  annotation_matrix <- matrix(nrow=length(accessions), ncol=length(annotations), dimnames=list(accessions,annotations))

  for(i in 1:length(accessions)){

    sel <- fData[,accession] == accessions[i]
    intensities <- exprs[sel,,drop=FALSE]
    properties <- fData[sel,]
    #If for the same accession in an annation column, there would be multiple different values, just paste them together.
    annotation_matrix[i,] <- apply(properties[,annotations], 2, function(x){paste0(unique(x))})
    properties <- properties[,-which(colnames(properties) %in% c(colnames(properties[,accession,drop=FALSE]),colnames(properties[,annotations,drop=FALSE]))), drop=FALSE]

    value <- as.vector(intensities)

    framelist <- replicate(ncol(intensities), droplevels(properties), simplify = FALSE)

    frame <- plyr::rbind.fill(framelist)

    frame2 <- cbind(data.frame(value=value,frame), pData[rep(row.names(pData), each=nrow(intensities)), ])

    #Removing NA
    frame3 <- frame2[complete.cases(frame2[,1]),]

    datalist[[i]] <- frame3
  }

  protdata <- new("protdata", accession=accessions, data=datalist, annotation=annotation_matrix)
  return(protdata)
}
