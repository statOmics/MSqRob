#' Import a MaxQuant peptides.txt file
#'
#' @description Imports a MaxQuant peptides.txt file and converts it into an \code{\link[=MSnSet-class]{MSnSet}} object.
#' @param file The name of a MaxQuant peptides.txt file. For more details about how this argument can be specified, see \code{\link[utils]{read.table}}.
#' @param pattern A character string containing a regular expression that will be matched to the file's header. The columns matching the expression should be the columns containing the peptide intensity values. Defaults to "Intensity.".
#' @param remove_pattern A logical indicating whether the expression in "pattern" should be removed from the column names in the resulting \code{\link[=MSnSet-class]{MSnSet}} object. Defaults to \code{TRUE}.
#' @return An object of class \code{\link[=MSnSet-class]{MSnSet}}.
#' @export
read_MaxQuant <- function(file, pattern="Intensity.", remove_pattern=TRUE)
{
  colInt <- MSnbase::grepEcols(file, pattern=pattern, split = "\t")
  peptides <- MSnbase::readMSnSet2(file, ecol = colInt, sep = "\t")

  if(isTRUE(remove_pattern)){
    #Remove pattern from colnames of exprs
    exprs <- Biobase::exprs(peptides)
    pData <- Biobase::pData(peptides)
    colnames(exprs) <- gsub(pattern,"",colnames(exprs))
    rownames(pData) <- gsub(pattern,"",rownames(pData))
    Biobase::sampleNames(Biobase::protocolData(peptides)) <- rownames(pData)
    Biobase::pData(peptides) <- pData
    Biobase::exprs(peptides) <- exprs
  }
  return(peptides)
}


#' Preprocess MaxQuant peptides.txt files
#'
#' @description Performs a standard preprocessing pipeline on \code{\link[=MSnSet-class]{MSnSet}} objects originating from MaxQuant peptides.txt files.
#' By default, intensity values are log2 transformed and then quantile normalized. Next, the \code{\link[=smallestUniqueGroups]{smallestUniqueGroups}} function is applied,
#' which removes proteins groups for which any of its member proteins is present in a smaller protein group. Then, contaminants and reverse sequences are removed.
#' Next, irrelevant columns are dropped. Then, peptide sequences that are identified only once in a single mass spec run are removed because with only 1 identification, the model will be perfectly confounded. Finally, potential experimental annotations are added to the data frame.
#' @param MSnSet An \code{\link[=MSnSet-class]{MSnSet}} object that contains data originating from MaxQuant's peptides.txt file.
#' @param accession A character indicating the column that contains the unit on which you want to do inference (typically the protein identifiers). Defaults to "Proteins".
#' @param exp_annotation Either the path to the file which contains the experiment annotation or a data frame containing the experiment annotation. Exactly one colum in the experiment annotation should contain the mass spec run names. Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}. As an error protection measurement, leading and trailing spaces in each column are trimmed off. The default, \code{NULL} indicates there is no annotation to be added.
#' @param type_annot If \code{exp_annotation} is a path to a file, the type of file. \code{type_annot} is mostly obsolete as supported files will be automatically recognized by their extension. Currently only \code{"tab-delim"} (tab-delimited file), \code{"xlsx"} (Office Open XML Spreadsheet file) and \code{NULL} (file type decided based on the extension) are supported. If the extension is not recognized, the file will be assumed to be a tab-delimited file. Defaults to \code{NULL}.
#' @param logtransform A logical value indicating whether the intensities should be log-transformed. Defaults to \code{TRUE}.
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Defaults to 2.
#' @param normalisation A character vector of length one that describes how to normalise the \code{\link[=MSnSet-class]{MSnSet}} object. See \code{\link[=normalise-methods]{normalise}} for details. Defaults to \code{"quantiles"}. If no normalisation is wanted, set \code{normalisation="none"}.
#' @param smallestUniqueGroups A logical indicating whether proteins groups for which any of its member proteins is present in a smaller protein group should be removed from the dataset. Defaults to \code{TRUE}.
#' @param useful_properties The columns of the \code{\link{featureData}} slot that are useful in the further analysis and/or inspection of the data and should be retained. Defaults to \code{c("Proteins","Sequence","PEP")}.
#' @param filter A vector of names corresponding to the columns in the \code{\link{featureData}} slot of the \code{\link[=MSnSet-class]{MSnSet}} object that contain a \code{filtersymbol} that indicates which rows should be removed from the data.
#' Typical examples are contaminants or reversed sequences. Defaults to \code{c("Contaminant","Reverse")}. Note that in earlier versions of MaxQuant the "Contaminant" column was called "Potential.contaminant", so one should use \code{c("Potential.contaminant","Reverse")} then instead.
#' @param filter_symbol A character indicating the symbol in the columns corresponding to the \code{filter} argument that is used to indicate rows that should be removed from the data. Defaults to "+".
#' @param minIdentified A numeric value indicating the minimal number of times a peptide sequence should be identified in the dataset in order not to be removed. Defaults to 2.
#' @param remove_only_site A logical indicating wheter proteins that are only identified by peptides carrying one or more modification sites should be removed from the data. This requires the extra input of a proteinGroups.txt file in the \code{file_proteinGroups} argument. Defaults to \code{FALSE}.
#' @param file_proteinGroups The name of the proteinGroups.txt file, which is used to remove proteins that are only identified by peptides carrying one or more modification sites. Only used when \code{remove_only_site} is set to \code{TRUE}.
#' @param colClasses character. Only used when the \code{exp_annotation} argument is a filepath. A vector of classes to be assumed for the columns of the experimental annotation data frame. Recycled if necessary. If named and shorter than required, names are matched to the column names with unspecified values are taken to be NA.
#' Possible values are \code{NA} (the default, when \code{type.convert} is used), \code{NULL} (when the column is skipped), one of the atomic vector classes (\code{logical}, \code{integer}, \code{numeric}, \code{complex}, \code{character}, \code{raw}), or \code{factor}, \code{Date} or \code{POSIXct}. Otherwise there needs to be an as method (from package \code{methods}) for conversion from \code{character} to the specified formal class.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @return A preprocessed \code{\link[=MSnSet-class]{MSnSet}} object that is ready to be converted into a \code{\link[=protdata-class]{protdata}} object.
#' @include preprocess_hlpFunctions.R
#' @export
preprocess_MaxQuant <- function(MSnSet, accession="Proteins", exp_annotation=NULL, type_annot=NULL, logtransform=TRUE, base=2, normalisation="quantiles", smallestUniqueGroups=TRUE, useful_properties=c("Proteins","Sequence","PEP"), filter=c("Contaminant","Reverse"), filter_symbol="+", minIdentified=2, remove_only_site=FALSE, file_proteinGroups=NULL, colClasses=NA, printProgress=FALSE, shiny=FALSE)
{

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Preprocessing...", value = 0)
  }

#Error control
if(!all(useful_properties %in% colnames(Biobase::fData(MSnSet)))){stop("Argument \"useful_properties\" must only contain column names of the featureData slot.")}
if(!all(filter %in% colnames(Biobase::fData(MSnSet)))){stop("One or more elements in the \"filter\" argument are no column names of the featureData slot of the MSnSet object.")}

#upDateProgress(progress=progress, detail="Log-transforming data", n=8, shiny=shiny, print=isTRUE(printProgress & logtransform))

if(isTRUE(logtransform)){
  #Log transform
  MSnSet <- log(MSnSet, base=base)
}

#Change -Inf values in the peptide intensities to NA
exprs <- Biobase::exprs(MSnSet)
exprs[is.infinite(exprs)] <- NA
Biobase::exprs(MSnSet) <- exprs

#upDateProgress(progress=progress, detail="Normalizing data", n=8, shiny=shiny, print=isTRUE(printProgress & (normalisation!="none")))

#Normalisation
if(normalisation!="none"){
  MSnSet <- tryCatch(
    MSnbase::normalise(MSnSet, normalisation), error=function(e){stop("\"normalisation\" argument should be one of \"none\", \"sum\", \"max\", \"center.mean\", \"center.median\", \"quantiles\", \"quantiles.robust\", \"vsn\"")})
}

#upDateProgress(progress=progress, detail="Removing overlapping protein groups", n=8, shiny=shiny, print=isTRUE(printProgress & smallestUniqueGroups))

#Our approach: a peptide can map to multiple proteins,
#as long as there is none of these proteins present in a smaller subgroup
if(isTRUE(smallestUniqueGroups)){
groups2 <- smallestUniqueGroups(Biobase::fData(MSnSet)[,accession], split = ";")
sel <- Biobase::fData(MSnSet)[,accession] %in% groups2
MSnSet <- MSnSet[sel]
}

#upDateProgress(progress=progress, detail="Removing contaminants and/or reverse sequences", n=8, shiny=shiny, print=isTRUE(printProgress & (length(filter)==0)))

#Remove contaminants and reverse sequences
if(!is.null(filter)){
filterdata <- Biobase::fData(MSnSet)[,filter, drop=FALSE]
#Sometimes, there are no contaminants or no reverse sequences, R then reads these empty columns as "NA"
filterdata[is.na(filterdata)] <- ""
sel <- rowSums(filterdata!= filter_symbol)==length(filter)
MSnSet <- MSnSet[sel]
}

#upDateProgress(progress=progress, detail="Removing proteins only identified by modified peptides", n=8, shiny=shiny, print=isTRUE(printProgress & remove_only_site))

#Remove only identified by site if proteinGroups.txt file is given
if(isTRUE(remove_only_site)){
  proteinGroups <- read.table(file_proteinGroups, sep="\t", header=TRUE, quote="", comment.char = "")
  only_site <- proteinGroups$Only.identified.by.site
  only_site[is.na(only_site)] <- ""
  removed_proteins <- proteinGroups$Protein.IDs[only_site==filter_symbol]

  sel <- !(as.character(Biobase::fData(MSnSet)[,accession]) %in% as.character(removed_proteins))
  MSnSet <- MSnSet[sel]
}

#upDateProgress(progress=progress, detail="Removing less usefull properties", n=8, shiny=shiny, print=isTRUE(printProgress))

#Retain only those properties in the fData slot that are useful (or might be useful) for our further analysis:
#This always includes the accession (protein) as well as the peptide identifier (almost always)
#If the accession is not present, add it
if(!(accession %in% useful_properties)){useful_properties <- c(accession,useful_properties)}
Biobase::fData(MSnSet) <- Biobase::fData(MSnSet)[,useful_properties]

#upDateProgress(progress=progress, detail=paste0("Removing peptides identified less than ", minIdentified," times"), n=8, shiny=shiny, print=isTRUE(printProgress))

#How many times shoud a peptide be identified?
#We require by default at least 2 identifications of a peptide sequence, as with 1 identification, the model will be perfectly confounded
keepers <- rowSums(!is.na(Biobase::exprs(MSnSet)))>=minIdentified
MSnSet <- MSnSet[keepers]

#upDateProgress(progress=progress, detail=paste0("Adding experimental annotation"), n=8, shiny=shiny, print=isTRUE(printProgress & !is.null(exp_annotation)))

#Add experiment annotation
if(!is.null(exp_annotation)){

  pData <- makeAnnotation(exp_annotation=exp_annotation, type_annot=type_annot, colClasses=colClasses)

  exprs <- Biobase::exprs(MSnSet)

  #Check which column of the given exp_annotation (pData) contains exactly the same elements as the mass spec run names in the data
  annotation_run <- getAnnotationRun(pData, colnames(exprs))

  #Error checks
  check_expAnn(pData, annotation_run)

  #Sort columns in exprs by annotation_run column
  exprs <- exprs[,match(as.character(pData[,annotation_run]), colnames(exprs))]
  rownames(pData) <- colnames(exprs)

  Biobase::sampleNames(Biobase::protocolData(MSnSet)) <- colnames(exprs)
  Biobase::exprs(MSnSet) <- exprs
  Biobase::pData(MSnSet) <- pData

}

return(MSnSet)
}




