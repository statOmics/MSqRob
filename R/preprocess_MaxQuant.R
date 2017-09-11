#' Import a file and convert it to an MSnSet object
#'
#' @description Imports a file and converts it into an \code{\link[=MSnSet-class]{MSnSet}} object (Gatto et al., 2012).
#' @param file The name of a file. For more details about how this argument can be specified, see \code{\link[utils]{read.table}}.
#' @param pattern A character string containing a regular expression that will be matched to the file's header. The columns matching the expression should be the columns containing the peptide intensity values.
#' @param remove_pattern A logical indicating whether the expression in "pattern" should be removed from the column names in the resulting \code{\link[=MSnSet-class]{MSnSet}} object. Defaults to \code{FALSE}.
#' @param sep A character indicating the field separator character. Values on each line of the file are separated by this character.The  Defaults to "\\t", indicating a tab-delimited input file.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @return An object of class \code{\link[=MSnSet-class]{MSnSet}}.
#' @references Gatto L, Lilley KS. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. \url{https://doi.org/10.1093/bioinformatics/btr645}.
#' @export
read2MSnSet <- function(file, pattern, remove_pattern=FALSE, sep="\t", shiny=FALSE, message=NULL)
{

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  colInt <- MSnbase::grepEcols(file, pattern=pattern, split = sep)
  peptides <- MSnbase::readMSnSet2(file, ecol = colInt, sep = sep)

  #Only now call make.names
  pattern <- make.names(pattern, unique = TRUE)

  if(isTRUE(remove_pattern)){
    #Remove pattern from colnames of exprs
    exprs <- Biobase::exprs(peptides)
    pData <- Biobase::pData(peptides)
    colnames(exprs) <- make.names(gsub(pattern,"",colnames(exprs)), unique = TRUE)
    rownames(pData) <- make.names(gsub(pattern,"",rownames(pData)), unique = TRUE)
    Biobase::sampleNames(Biobase::protocolData(peptides)) <- rownames(pData)
    Biobase::pData(peptides) <- pData
    Biobase::exprs(peptides) <- exprs
  }
  return(peptides)
}


#' Import a MaxQuant peptides.txt file
#'
#' @description Imports a MaxQuant (Cox and Mann, 2008) peptides.txt file and converts it into an \code{\link[=MSnSet-class]{MSnSet}} object (Gatto et al., 2012).
#' @param file The name of a MaxQuant peptides.txt file. For more details about how this argument can be specified, see \code{\link[utils]{read.table}}.
#' @param pattern A character string containing a regular expression that will be matched to the file's header. The columns matching the expression should be the columns containing the peptide intensity values. Defaults to "Intensity ".
#' @param remove_pattern A logical indicating whether the expression in "pattern" should be removed from the column names in the resulting \code{\link[=MSnSet-class]{MSnSet}} object. Defaults to \code{TRUE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @return An object of class \code{\link[=MSnSet-class]{MSnSet}}.
#' @references Gatto L, Lilley KS. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. \url{https://doi.org/10.1093/bioinformatics/btr645}.
#' @references Cox, J. and Mann, M. MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. Nat Biotechnol, 2008, 26, pp 1367-72. \url{http://www.nature.com/nbt/journal/v26/n12/full/nbt.1511.html}.
#' @export
read_MaxQuant <- function(file, pattern="Intensity ", remove_pattern=TRUE, shiny=FALSE, message=NULL)
{
  peptides <- read2MSnSet(file=file, pattern=pattern, remove_pattern=remove_pattern, sep="\t", shiny=shiny, message=message)
  return(peptides)
}


#' Import a peptides summary file produced by moFF
#'
#' @description Imports a moFF (Argentini et al., 2016) .tab file and converts it into an \code{\link[=MSnSet-class]{MSnSet}} object (Gatto et al., 2012).
#' @param file The name of a moFF .tab file. For more details about how this argument can be specified, see \code{\link[utils]{read.table}}.
#' @param pattern A character string containing a regular expression that will be matched to the file's header. The columns matching the expression should be the columns containing the peptide intensity values. Defaults to "?????".
#' @param remove_pattern A logical indicating whether the expression in "pattern" should be removed from the column names in the resulting \code{\link[=MSnSet-class]{MSnSet}} object. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @return An object of class \code{\link[=MSnSet-class]{MSnSet}}.
#' @references Gatto L, Lilley KS. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. \url{https://doi.org/10.1093/bioinformatics/btr645}.
#' @references Argentini A,	Goeminne LJE,	Verheggen K,	Hulstaert N,	Staes A, Clement L	& Martens L. moFF: a robust and automated approach to extract peptide ion intensities. Nature Methods. 2016 13:964–966.  \url{http://www.nature.com/nmeth/journal/v13/n12/full/nmeth.4075.html}.
#' @export
read_moFF <- function(file, pattern="sumIntensity_", remove_pattern=TRUE, shiny=FALSE, message=NULL)
{
  peptides <- read2MSnSet(file=file, pattern=pattern, remove_pattern=remove_pattern, sep="\t", shiny=shiny, message=message)
  return(peptides)
}


#' Preprocess MSnSet objects
#'
#' @description This function allows to perform a standard preprocessing pipeline on \code{\link[=MSnSet-class]{MSnSet}} objects (Gatto et al., 2012).
#' By default, intensity values are log2 transformed and then quantile normalized. Next, the \code{\link[=smallestUniqueGroups]{smallestUniqueGroups}} function is applied,
#' which removes proteins groups for which any of its member proteins is present in a smaller protein group. Then, peptides that need to be filtered out are removed.
#' Next, irrelevant columns are dropped. Then, peptide sequences that are identified only once in a single mass spec run are removed because with only 1 identification, the model will be perfectly confounded. Finally, potential experimental annotations are added to the data frame.
#' @param MSnSet An \code{\link[=MSnSet-class]{MSnSet}} object.
#' @param accession A character indicating the column that contains the unit on which you want to do inference (typically the protein identifiers).
#' @param exp_annotation Either the path to the file which contains the experiment annotation or a data frame containing the experiment annotation. Exactly one colum in the experiment annotation should contain the mass spec run names. Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}. As an error protection measurement, leading and trailing spaces in each column are trimmed off. The default, \code{NULL} indicates there is no annotation to be added.
#' @param type_annot If \code{exp_annotation} is a path to a file, the type of file. \code{type_annot} is mostly obsolete as supported files will be automatically recognized by their extension. Currently only \code{"tab-delim"} (tab-delimited file), \code{"xlsx"} (Office Open XML Spreadsheet file) and \code{NULL} (file type decided based on the extension) are supported. If the extension is not recognized, the file will be assumed to be a tab-delimited file. Defaults to \code{NULL}.
#' @param logtransform A logical value indicating whether the intensities should be log-transformed. Defaults to \code{TRUE}.
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Defaults to 2.
#' @param normalisation A character vector of length one that describes how to normalise the \code{\link[=MSnSet-class]{MSnSet}} object. See \code{\link[=normalise-methods]{normalise}} for details. Defaults to \code{"quantiles"}. If no normalisation is wanted, set \code{normalisation="none"}.
#' @param weights Only used when \code{normalisation} is set to or "rlr", "loess.fast", "loess.affy" or "loess.pairs". A numeric vector of weights for each row in the MSnSet object to be used for the fitting during the normalisation step. Defaults to \code{NULL}.
#' @param smallestUniqueGroups A logical indicating whether protein groups for which any of its member proteins is present in a smaller protein group should be removed from the dataset. Defaults to \code{TRUE}.
#' @param split A character string that indicates the separator between protein groups. Only used when \code{smallestUniqueGroups} is set to \code{TRUE}.
#' @param useful_properties The columns of the \code{\link{featureData}} slot that are useful in the further analysis and/or inspection of the data and should be retained. Defaults to \code{NULL}, in which case no additional columns will be retained.
#' @param filter A vector of names corresponding to the columns in the \code{\link{featureData}} slot of the \code{\link[=MSnSet-class]{MSnSet}} object that contain a \code{filtersymbol} that indicates which rows should be removed from the data.
#' Typical examples are contaminants or reversed sequences. Defaults to \code{NULL}, in which case no filtering will be performed.
#' @param filter_symbol Only used when \code{filter} is not \code{NULL}. A character indicating the symbol in the columns corresponding to the \code{filter} argument that is used to indicate rows that should be removed from the data. Defaults to \code{NULL}, which will throw an error if \code{filter} is not \code{NULL} to alert the user to specify a filter symbol.
#' @param minIdentified A numeric value indicating the minimal number of times a peptide sequence should be identified in the dataset in order not to be removed. Defaults to 2.
#' @param external_filter_file The name of an external protein filtering file. Sometimes, users want to filter out proteins based on a separate protein file. This file should contain at least a column with name equal to the value in \code{external_filter_accession} containing proteins, and one or more columns on which to filter, with names equal to the input in \code{external_filter_column}. Proteins that need to be filtered out should have the \code{filter_symbol} in their \code{external_filter_column}. Defaults to \code{NULL}, in which case no filtering based on an external protein file will be done.
#' @param external_filter_accession Only used when \code{external_filter_file} is not \code{NULL}. A character indicating the column that contains the protein identifiers in the \code{external_filter_file}. Defaults to \code{NULL}, which will throw an error if \code{external_filter_file} is not \code{NULL} to alert the user to specify a filter column.
#' @param external_filter_column Only used when \code{external_filter_file} is not \code{NULL}. A vector of names containing the column name(s) on which to filter in the \code{external_filter_file}. Defaults to \code{NULL}, which will throw an error if \code{external_filter_file} is not \code{NULL} to alert the user to specify a filter column.
#' @param colClasses character. Only used when the \code{exp_annotation} argument is a filepath. A vector of classes to be assumed for the columns of the experimental annotation data frame. Recycled if necessary. If named and shorter than required, names are matched to the column names with unspecified values are taken to be NA.
#' Possible values are \code{"keep"} (the default, when the colClasses are unchanged for data frames and \code{type.convert} is used for files),  \code{NA} (when \code{type.convert} is always used), \code{NULL} (when the column is skipped), one of the atomic vector classes (\code{"logical"}, \code{"integer"}, \code{"numeric"}, \code{"complex"}, \code{"character"}, \code{"raw"}), or \code{"factor"}, \code{"Date"} or \code{"POSIXct"}. Otherwise there needs to be an as method (from package \code{methods}) for conversion from \code{"character"} to the specified formal class.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @param details Only used when \code{shiny=TRUE} or \code{printProgress=TRUE}. A character vector containing the detail messages to be displayed to the user, or \code{NULL} to hide the current detail messages (if any). The detail messages will be shown with a de-emphasized appearance relative to the message.
#' @return A preprocessed \code{\link[=MSnSet-class]{MSnSet}} object that is ready to be converted into a \code{\link[=protdata-class]{protdata}} object.
#' @references Gatto L, Lilley KS. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. \url{https://doi.org/10.1093/bioinformatics/btr645}. PubMed PMID:22113085.
#' @include preprocess_hlpFunctions.R
#' @include updateProgress.R
#' @export
preprocess_MSnSet <- function(MSnSet, accession, exp_annotation=NULL, type_annot=NULL, logtransform=TRUE, base=2, normalisation="quantiles", weights=NULL, smallestUniqueGroups=TRUE, split=NULL, useful_properties=NULL, filter=NULL, filter_symbol=NULL, minIdentified=2, external_filter_file=NULL, external_filter_accession=NULL, external_filter_column=NULL, colClasses="keep", printProgress=FALSE, shiny=FALSE, message=NULL, details=NULL)
{
  #Error control

  accession <- make.names(accession, unique = TRUE)
  useful_properties <- make.names(useful_properties, unique = TRUE)

  if(isTRUE(smallestUniqueGroups) && is.null(split)){stop("Please provide the protein groups separator (split argument) or set the smallestUniqueGroups argument to FALSE.")}

  if(!(accession %in% useful_properties)){useful_properties <- c(accession,useful_properties)}

  if(!all(useful_properties %in% colnames(Biobase::fData(MSnSet)))){stop("Argument \"useful_properties\" must only contain column names of the featureData slot.")}
  if(!all(filter %in% colnames(Biobase::fData(MSnSet)))){stop("One or more elements in the \"filter\" argument are no column names of the featureData slot of the MSnSet object.")}

  n <- sum(isTRUE(logtransform),(normalisation!="none"),isTRUE(smallestUniqueGroups),
           !is.null(filter),!is.null(external_filter_file),
           !all(colnames(Biobase::fData(MSnSet)) %in% useful_properties),
           minIdentified>1,!is.null(exp_annotation)
           )

  progress <- NULL
  if(isTRUE(shiny) && n>0){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  if(isTRUE(logtransform)){
    updateProgress(progress=progress, detail=details[1], n=n, shiny=shiny, print=isTRUE(printProgress & logtransform))

    #Log transform
    MSnSet <- log(MSnSet, base=base)
  }

  #Change -Inf values in the peptide intensities to NA
  #-Inf is only to be found in data that is already transformed => can safely stay here, even if the data has not been transformed in advance!
  exprs <- Biobase::exprs(MSnSet)
  exprs[is.infinite(exprs)] <- NA
  Biobase::exprs(MSnSet) <- exprs

  #Normalisation
  if(normalisation!="none"){
    updateProgress(progress=progress, detail=details[2], n=n, shiny=shiny, print=isTRUE(printProgress & (normalisation!="none")))
    MSnSet <- .normaliseMSnSet(MSnSet, normalisation, weights)
  }

  #Our approach: a peptide can map to multiple proteins,
  #as long as there is none of these proteins present in a smaller subgroup
  if(isTRUE(smallestUniqueGroups)){
    updateProgress(progress=progress, detail=details[3], n=n, shiny=shiny, print=isTRUE(printProgress & smallestUniqueGroups))

    groups2 <- smallestUniqueGroups(Biobase::fData(MSnSet)[,accession], split = split)
    sel <- Biobase::fData(MSnSet)[,accession] %in% groups2
    MSnSet <- MSnSet[sel]
  }

  #Remove contaminants and reverse sequences
  if(!is.null(filter)){
    updateProgress(progress=progress, detail=details[4], n=n, shiny=shiny, print=isTRUE(printProgress & (length(filter)!=0)))

    filterdata <- Biobase::fData(MSnSet)[,filter, drop=FALSE]
    #Sometimes, there are no contaminants or no reverse sequences, R then reads these empty columns as "NA"
    filterdata[is.na(filterdata)] <- ""
    sel <- rowSums(filterdata!= filter_symbol)==length(filter)
    MSnSet <- MSnSet[sel]
  }

  #Remove only identified by site if proteinGroups.txt file is given
  if(!is.null(external_filter_file)){
    updateProgress(progress=progress, detail=details[5], n=n, shiny=shiny, print=isTRUE(printProgress))

    externalFilter <- read.table(external_filter_file, sep="\t", header=TRUE, quote="", comment.char = "")
    only_site <- externalFilter[[external_filter_column]]
    only_site[is.na(only_site)] <- ""
    removed_proteins <- externalFilter[[external_filter_accession]][only_site==filter_symbol]

    sel <- !(as.character(Biobase::fData(MSnSet)[,accession]) %in% as.character(removed_proteins))
    MSnSet <- MSnSet[sel]
  }

  if(!all(colnames(Biobase::fData(MSnSet)) %in% useful_properties)){
  updateProgress(progress=progress, detail=details[6], n=n, shiny=shiny, print=isTRUE(printProgress))

  #Retain only those properties in the fData slot that are useful (or might be useful) for our further analysis:
  #This always includes the accession (protein) as well as the peptide identifier (almost always)
  #If the accession is not present, add it
  Biobase::fData(MSnSet) <- Biobase::fData(MSnSet)[,useful_properties, drop=FALSE]
  }

  if(minIdentified>1){
  updateProgress(progress=progress, detail=details[7], n=n, shiny=shiny, print=isTRUE(printProgress))

  #How many times shoud a peptide be identified?
  #We require by default at least 2 identifications of a peptide sequence, as with 1 identification, the model will be perfectly confounded
  keepers <- rowSums(!is.na(Biobase::exprs(MSnSet)))>=minIdentified
  MSnSet <- MSnSet[keepers]
  }

  #Add experiment annotation
  if(!is.null(exp_annotation)){

    updateProgress(progress=progress, detail=details[8], n=n, shiny=shiny, print=isTRUE(printProgress & !is.null(exp_annotation)))

    exprs <- Biobase::exprs(MSnSet)

    pData <- makeAnnotation(exp_annotation=exp_annotation, run_names=colnames(exprs), type_annot=type_annot, colClasses=colClasses)

    #Sort columns in exprs by annotation_run column
    annotation_run <- getAnnotationRun(pData=pData, run_names=colnames(exprs))
    exprs <- exprs[,match(as.character(pData[,annotation_run]), colnames(exprs))]
    rownames(pData) <- colnames(exprs)

    # Biobase::sampleNames(Biobase::protocolData(MSnSet)) <- colnames(exprs)
    # environment(exprs) <- MSnSet@assayData
    # Biobase::exprs(MSnSet) <- exprs
    # #Important check:
    # if(!all(colnames(exprs)==colnames(Biobase::exprs(MSnSet)))){stop("Biobase error, cannot change exprs slot...")}
    # Biobase::pData(MSnSet) <- pData
    # rm(MSnSet)
    MSnSet <- MSnbase::MSnSet(exprs=exprs, fData=Biobase::fData(MSnSet), pData=pData)

  }

  return(MSnSet)
}


#' Preprocess MSnSet objects originating from MaxQuant peptides.txt files
#'
#' @description Performs a standard preprocessing pipeline on \code{\link[=MSnSet-class]{MSnSet}} objects (Gatto et al., 2012) originating from MaxQuant (Cox and Mann, 2008) peptides.txt files.
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
#' @param weights Only used when \code{normalisation} is set to or "rlr", "loess.fast", "loess.affy" or "loess.pairs". A numeric vector of weights for each row in the MSnSet object to be used for the fitting during the normalisation step. Defaults to \code{NULL}.
#' @param smallestUniqueGroups A logical indicating whether protein groups for which any of its member proteins is present in a smaller protein group should be removed from the dataset. Defaults to \code{TRUE}.
#' @param useful_properties The columns of the \code{\link{featureData}} slot that are useful in the further analysis and/or inspection of the data and should be retained. Defaults to \code{c("Proteins","Sequence","PEP")}.
#' @param filter A vector of names corresponding to the columns in the \code{\link{featureData}} slot of the \code{\link[=MSnSet-class]{MSnSet}} object that contain a \code{filtersymbol} that indicates which rows should be removed from the data.
#' Typical examples are contaminants or reversed sequences. Defaults to \code{c("Contaminant","Reverse")}. Note that in earlier versions of MaxQuant the "Contaminant" column was called "Potential.contaminant". If "Potential.contaminant" is mentioned in this argument but could not be found, this function automatically tries to filter on "Contaminant".
#' @param filter_symbol A character indicating the symbol in the columns corresponding to the \code{filter} argument that is used to indicate rows that should be removed from the data. Defaults to "+".
#' @param minIdentified A numeric value indicating the minimal number of times a peptide sequence should be identified in the dataset in order not to be removed. Defaults to 2.
#' @param remove_only_site A logical indicating wheter proteins that are only identified by peptides carrying one or more modification sites should be removed from the data. This requires the extra input of a proteinGroups.txt file in the \code{file_proteinGroups} argument. Defaults to \code{FALSE}.
#' @param file_proteinGroups The name of the proteinGroups.txt file, which is used to remove proteins that are only identified by peptides carrying one or more modification sites. Only used when \code{remove_only_site} is set to \code{TRUE}.
#' @param colClasses character. Only used when the \code{exp_annotation} argument is a filepath. A vector of classes to be assumed for the columns of the experimental annotation data frame. Recycled if necessary. If named and shorter than required, names are matched to the column names with unspecified values are taken to be NA.
#' Possible values are \code{"keep"} (the default, when the colClasses are unchanged for data frames and \code{type.convert} is used for files),  \code{NA} (when \code{type.convert} is always used), \code{NULL} (when the column is skipped), one of the atomic vector classes (\code{"logical"}, \code{"integer"}, \code{"numeric"}, \code{"complex"}, \code{"character"}, \code{"raw"}), or \code{"factor"}, \code{"Date"} or \code{"POSIXct"}. Otherwise there needs to be an as method (from package \code{methods}) for conversion from \code{"character"} to the specified formal class.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @return A preprocessed \code{\link[=MSnSet-class]{MSnSet}} object that is ready to be converted into a \code{\link[=protdata-class]{protdata}} object.
#' @references Gatto L, Lilley KS. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. \url{https://doi.org/10.1093/bioinformatics/btr645}. PubMed PMID:22113085.
#' @references Cox, J. and Mann, M. MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. Nat Biotechnol, 2008, 26, pp 1367-72. \url{http://www.nature.com/nbt/journal/v26/n12/full/nbt.1511.html}.
#' @include preprocess_hlpFunctions.R
#' @include updateProgress.R
#' @export
preprocess_MaxQuant <- function(MSnSet, accession="Proteins", exp_annotation=NULL, type_annot=NULL, logtransform=TRUE, base=2, normalisation="quantiles", weights=NULL, smallestUniqueGroups=TRUE, useful_properties=c("Proteins","Sequence","PEP"), filter=c("Potential.contaminant","Reverse"), filter_symbol="+", minIdentified=2, remove_only_site=FALSE, file_proteinGroups=NULL, colClasses="keep", printProgress=FALSE, shiny=FALSE, message=NULL)
{

  #Some older versions of MaxQuant use "Contaminant" instead of "Potential.contaminant"
  #Condition "Potential.contaminant" %in% filter added to prevent filter=NULL turn into character(0)!!!
  if("Potential.contaminant" %in% filter && !("Potential.contaminant" %in% colnames(Biobase::fData(MSnSet)))){filter[filter=="Potential.contaminant"] <- "Contaminant"}

  details <- c("Log-transforming data",
               "Normalizing data",
               "Removing overlapping protein groups",
               "Removing contaminants and/or reverse sequences",
               "Removing proteins only identified by modified peptides",
               "Removing less useful properties",
               paste0("Removing peptides identified less than ", minIdentified," times"),
               "Adding experimental annotation")

  external_filter_accession="Protein.IDs"
  external_filter_column="Only.identified.by.site"

  if(!isTRUE(remove_only_site)){file_proteinGroups <- NULL}

  MSnSet <- preprocess_MSnSet(MSnSet=MSnSet, accession=accession, exp_annotation=exp_annotation, type_annot=type_annot, logtransform=logtransform, base=base, normalisation=normalisation, weights=weights, smallestUniqueGroups=smallestUniqueGroups, split=";", useful_properties=useful_properties, filter=filter, filter_symbol=filter_symbol, minIdentified=minIdentified,
                              external_filter_file=file_proteinGroups, external_filter_accession=external_filter_accession, external_filter_column=external_filter_column, colClasses=colClasses, printProgress=printProgress, shiny=shiny, message=message, details=details)

  #If pData is completely empty due to a lack of annotation, at least add the runs
  #This feature is not included in preprocess_MSnSet, where you have more liberty to play with all kinds of preprocessing
  if(ncol(pData(MSnSet))==0){
  emptyPData <- data.frame(run=rownames(pData(MSnSet)))
  rownames(emptyPData) <- rownames(pData(MSnSet))
  pData(MSnSet) <- emptyPData
  }

  return(MSnSet)
}


#' Preprocess MSnSet objects originating from moFF .tab files
#'
#' @description This function allows to perform a standard preprocessing pipeline on \code{\link[=MSnSet-class]{MSnSet}} objects (Gatto et al., 2012) originating from moFF .tab files (Argentini et al., 2016).
#' By default, intensity values are log2 transformed and then quantile normalized. Next, the \code{\link[=smallestUniqueGroups]{smallestUniqueGroups}} function is applied,
#' which removes proteins groups for which any of its member proteins is present in a smaller protein group. Then, peptides that need to be filtered out are removed.
#' Next, irrelevant columns are dropped. Then, peptide sequences that are identified only once in a single mass spec run are removed because with only 1 identification, the model will be perfectly confounded. Finally, potential experimental annotations are added to the data frame.
#' @param MSnSet An \code{\link[=MSnSet-class]{MSnSet}} object that originates from a moFF .tab file.
#' @param accession A character indicating the column that contains the unit on which you want to do inference (typically the protein identifiers).
#' @param exp_annotation Either the path to the file which contains the experiment annotation or a data frame containing the experiment annotation. Exactly one colum in the experiment annotation should contain the mass spec run names. Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}. As an error protection measurement, leading and trailing spaces in each column are trimmed off. The default, \code{NULL} indicates there is no annotation to be added.
#' @param type_annot If \code{exp_annotation} is a path to a file, the type of file. \code{type_annot} is mostly obsolete as supported files will be automatically recognized by their extension. Currently only \code{"tab-delim"} (tab-delimited file), \code{"xlsx"} (Office Open XML Spreadsheet file) and \code{NULL} (file type decided based on the extension) are supported. If the extension is not recognized, the file will be assumed to be a tab-delimited file. Defaults to \code{NULL}.
#' @param logtransform A logical value indicating whether the intensities should be log-transformed. Defaults to \code{TRUE}.
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Defaults to 2.
#' @param normalisation A character vector of length one that describes how to normalise the \code{\link[=MSnSet-class]{MSnSet}} object. See \code{\link[=normalise-methods]{normalise}} for details. Defaults to \code{"quantiles"}. If no normalisation is wanted, set \code{normalisation="none"}.
#' @param weights Only used when \code{normalisation} is set to or "rlr", "loess.fast", "loess.affy" or "loess.pairs". A numeric vector of weights for each row in the MSnSet object to be used for the fitting during the normalisation step. Defaults to \code{NULL}.
#' @param smallestUniqueGroups A logical indicating whether protein groups for which any of its member proteins is present in a smaller protein group should be removed from the dataset. Defaults to \code{TRUE}.
#' @param useful_properties The columns of the \code{\link{featureData}} slot that are useful in the further analysis and/or inspection of the data and should be retained. Defaults to \code{NULL}, in which case no additional columns will be retained.
#' @param filter A vector of names corresponding to the columns in the \code{\link{featureData}} slot of the \code{\link[=MSnSet-class]{MSnSet}} object that contain a \code{filtersymbol} that indicates which rows should be removed from the data.
#' Typical examples are contaminants or reversed sequences. Defaults to \code{NULL}, in which case no filtering will be performed.
#' @param filter_symbol Only used when \code{filter} is not \code{NULL}. A character indicating the symbol in the columns corresponding to the \code{filter} argument that is used to indicate rows that should be removed from the data. Defaults to \code{NULL}, which will throw an error if \code{filter} is not \code{NULL} to alert the user to specify a filter symbol.
#' @param minIdentified A numeric value indicating the minimal number of times a peptide sequence should be identified in the dataset in order not to be removed. Defaults to 2.
#' @param external_filter_file The name of an external protein filtering file. Sometimes, users want to filter out proteins based on a separate protein file. This file should contain at least a column with name equal to the value in \code{external_filter_accession} containing proteins, and one or more columns on which to filter, with names equal to the input in \code{external_filter_column}. Proteins that need to be filtered out should have the \code{filter_symbol} in their \code{external_filter_column}. Defaults to \code{NULL}, in which case no filtering based on an external protein file will be done.
#' @param external_filter_accession Only used when \code{external_filter_file} is not \code{NULL}. A character indicating the column that contains the protein identifiers in the \code{external_filter_file}. Defaults to \code{NULL}, which will throw an error if \code{external_filter_file} is not \code{NULL} to alert the user to specify a filter column.
#' @param external_filter_column Only used when \code{external_filter_file} is not \code{NULL}. A vector of names containing the column name(s) on which to filter in the \code{external_filter_file}. Defaults to \code{NULL}, which will throw an error if \code{external_filter_file} is not \code{NULL} to alert the user to specify a filter column.
#' @param colClasses character. Only used when the \code{exp_annotation} argument is a filepath. A vector of classes to be assumed for the columns of the experimental annotation data frame. Recycled if necessary. If named and shorter than required, names are matched to the column names with unspecified values are taken to be NA.
#' Possible values are \code{"keep"} (the default, when the colClasses are unchanged for data frames and \code{type.convert} is used for files),  \code{NA} (when \code{type.convert} is always used), \code{NULL} (when the column is skipped), one of the atomic vector classes (\code{"logical"}, \code{"integer"}, \code{"numeric"}, \code{"complex"}, \code{"character"}, \code{"raw"}), or \code{"factor"}, \code{"Date"} or \code{"POSIXct"}. Otherwise there needs to be an as method (from package \code{methods}) for conversion from \code{"character"} to the specified formal class.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @param details Only used when \code{shiny=TRUE} or \code{printProgress=TRUE}. A character vector containing the detail messages to be displayed to the user, or \code{NULL} to hide the current detail messages (if any). The detail messages will be shown with a de-emphasized appearance relative to the message.
#' @return A preprocessed \code{\link[=MSnSet-class]{MSnSet}} object that is ready to be converted into a \code{\link[=protdata-class]{protdata}} object.
#' @references Gatto L, Lilley KS. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. \url{https://doi.org/10.1093/bioinformatics/btr645}. PubMed PMID:22113085.
#' @references Argentini A,	Goeminne LJE,	Verheggen K,	Hulstaert N,	Staes A, Clement L	& Martens L. moFF: a robust and automated approach to extract peptide ion intensities. Nature Methods. 2016 13:964–966.  \url{http://www.nature.com/nmeth/journal/v13/n12/full/nmeth.4075.html}.
#' @include preprocess_hlpFunctions.R
#' @include updateProgress.R
#' @export
preprocess_moFF <- function(MSnSet, accession="prot", exp_annotation=NULL, type_annot=NULL, logtransform=TRUE, base=2, normalisation="quantiles", weights=NULL, smallestUniqueGroups=TRUE, useful_properties="peptide", filter=NULL, filter_symbol=NULL, minIdentified=2, external_filter_file=NULL, external_filter_accession=NULL, external_filter_column=NULL, colClasses="keep", printProgress=FALSE, shiny=FALSE, message=NULL){

  details <- c("Log-transforming data",
               "Normalizing data",
               "Removing overlapping protein groups",
               "Filtering",
               "Filtering on external file",
               "Removing less useful properties",
               paste0("Removing peptides identified less than ", minIdentified," times"),
               "Adding experimental annotation")

  MSnSet <- preprocess_MSnSet(MSnSet=MSnSet, accession=accession, exp_annotation=exp_annotation, type_annot=type_annot, logtransform=logtransform, base=base, normalisation=normalisation, weights=weights, smallestUniqueGroups=smallestUniqueGroups, split=", ", useful_properties=useful_properties, filter=filter, filter_symbol=filter_symbol, minIdentified=minIdentified,
                              external_filter_file=external_filter_file, external_filter_accession=external_filter_accession, external_filter_column=external_filter_column, colClasses=colClasses, printProgress=printProgress, shiny=shiny, message=message, details=details)

  #If pData is completely empty due to a lack of annotation, at least add the runs
  #This feature is not included in preprocess_MSnSet, where you have more liberty to play with all kinds of preprocessing
  if(ncol(pData(MSnSet))==0){
    emptyPData <- data.frame(run=rownames(pData(MSnSet)))
    rownames(emptyPData) <- rownames(pData(MSnSet))
    pData(MSnSet) <- emptyPData
  }

  return(MSnSet)
}


.normaliseMSnSet <- function(MSnSet, normalisation, weights=NULL){

  if(!(normalisation) %in% c("none", "sum", "max", "center.mean", "center.median", "quantiles", "quantiles.robust", "vsn", "rlr", "loess.fast", "loess.affy", "loess.pairs")){
    stop("\"normalisation\" argument should be one of \"none\", \"sum\", \"max\", \"center.mean\", \"center.median\", \"quantiles\", \"quantiles.robust\", \"vsn\", \"rlr\", \"loess.fast\", \"loess.affy\", \"loess.pairs\"")
  } else if(normalisation %in% c("rlr", "loess.fast", "loess.affy", "loess.pairs")){ #WARNING: "loess.affy" and "loess.pairs" are deprecated because they remove all the NA values!!!

    exprs <- exprs(MSnSet)

    loess_choices <- c("fast", "affy", "pairs")
    names(loess_choices) <- c("loess.fast", "loess.affy", "loess.pairs")

    if(normalisation %in% names(loess_choices)){
      exprsNorm <- limma::normalizeCyclicLoess(exprs, method = loess_choices[normalisation], weights=weights)
    } else{ #if rlr
      exprsNorm <- .normalizeRLR(exprs)
    }

    dimnames(exprsNorm) <- dimnames(exprs)

    Biobase::exprs(MSnSet) <- exprsNorm

  } else{
    MSnSet <- MSnbase::normalise(MSnSet, normalisation)
  }
  return(MSnSet)
}

.normalizeRLR <- function(exprs, weights=NULL){
  mediandata <- apply(exprs, 1, "median", na.rm = TRUE)
  flag1 = 1
  for (j in 1:ncol(exprs)) {
    LRfit <- MASS::rlm(as.matrix(exprs[, j]) ~ mediandata, weights=weights, na.action = na.exclude)
    Coeffs <- LRfit$coefficients
    a <- Coeffs[2]
    b <- Coeffs[1]
    if (flag1 == 1) {
      globalfittedRLR <- (exprs[, j] - b)/a
      flag1 = 2
    }
    else {
      globalfittedRLR <- cbind(globalfittedRLR, (exprs[, j] - b)/a)
    }
  }
  return(globalfittedRLR)
}



