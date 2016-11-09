#' Preprocess data in "wide" format
#'
#' @description Performs a standard preprocessing pipeline on data frames in "wide" format (i.e. the data frame has one observation row per subject with each measurement present as a different variable).
#' By default, intensity values are log2 transformed and then quantile normalized. Next, the \code{\link[=smallestUniqueGroups]{smallestUniqueGroups}} function is applied,
#' which removes proteins groups for which any of its member proteins is present in a smaller protein group. Then, unwanted sequences (such as reverse sequences or unwanted sequences) are filtered out.
#' Next, irrelevant columns are dropped. Then, peptide sequences that are identified only once in a single mass spec run are removed because with only 1 identification, the model will be perfectly confounded. Finally, potential experimental annotations are added to the data frame.
#' @param df A data frame that contains data in "wide" format.
#' @param accession A character indicating the column that contains the unit on which you want to do inference (typically the protein identifiers).
#' @param split A character indicating which string is used to separate accession groups.
#' @param exp_annotation Either the path to the file which contains the experiment annotation or a data frame containing the experiment annotation. Exactly one colum in the experiment annotation should contain the mass spec run names. Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}. As an error protection measurement, leading and trailing spaces in each column are trimmed off. The default, \code{NULL} indicates there is no annotation to be added.
#' @param type_annot If \code{exp_annotation} is a path to a file, the type of file. \code{type_annot} is mostly obsolete as supported files will be automatically recognized by their extension. Currently only \code{"tab-delim"} (tab-delimited file), \code{"xlsx"} (Office Open XML Spreadsheet file) and \code{NULL} (file type decided based on the extension) are supported. If the extension is not recognized, the file will be assumed to be a tab-delimited file. Defaults to \code{NULL}.
#' @param quant_cols Either a character or numeric vector indicating the columns that contain the quantitative values of interest (mostly peptide intensities or peptide areas under the curve) or a character string of length one indicating a pattern that is unique for the column names of the columns that contain the quantitative values of interest.
#' @param aggr_by A character indicating the column by which the data should be aggregated. We advise to aggregate the data by peptide sequence (thus aggregate over different charge states and modification statuses of the same peptide). If you only want to aggregate over charge states, set \code{aggr_by} to the column corresponding to the modified sequences. If no aggregation at all is desired, create a new column in which you paste together the modified sequences and the charge states.
#' @param aggr_function The function used to aggregate intensity data. Defaults to \code{"sum"}.
#' @param logtransform A logical value indicating whether the intensities should be log-transformed. Defaults to \code{TRUE}.
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Defaults to 2.
#' @param normalisation A character vector of length one that describes how to normalise the data frame \code{df}. See \code{\link[=normalise-methods]{normalise}} for details. Defaults to \code{"quantiles"}. If no normalisation is wanted, set \code{normalisation="none"}.
#' @param smallestUniqueGroups A logical indicating whether proteins groups for which any of its member proteins is present in a smaller protein group should be removed from the dataset. Defaults to \code{TRUE}.
#' @param useful_properties The columns of the data frame \code{df} that are useful in the further analysis and/or inspection of the data and should be retained. Defaults to \code{NULL}, indicating no extra columns should be retained.
#' @param filter A vector of names corresponding to the columns in the data frame \code{df} that contain a \code{filtersymbol} that indicates which rows should be removed from the data.
#' Typical examples are contaminants or reversed sequences. Defaults to \code{NULL}, indicating no filtering should be applied.
#' @param filter_symbol A character indicating the symbol in the columns corresponding to the \code{filter} argument that is used to indicate rows that should be removed from the data. Defaults to \code{NULL}.
#' @param minIdentified A numeric value indicating the minimal number of times a peptide sequence should be identified in the dataset in order not to be removed. Defaults to 2.
#' @param colClasses character. Only used when the \code{exp_annotation} argument is a filepath. A vector of classes to be assumed for the columns of the experimental annotation data frame. Recycled if necessary. If named and shorter than required, names are matched to the column names with unspecified values are taken to be NA.
#' Possible values are \code{NA} (the default, when \code{type.convert} is used), \code{NULL} (when the column is skipped), one of the atomic vector classes (\code{logical}, \code{integer}, \code{numeric}, \code{complex}, \code{character}, \code{raw}), or \code{factor}, \code{Date} or \code{POSIXct}. Otherwise there needs to be an as method (from package \code{methods}) for conversion from \code{character} to the specified formal class.
#' @param printProgress A logical indicating whether the R should print a message before performing each preprocessing step. Defaults to \code{FALSE}.
#' @param shiny A logical indicating whether this function is being used by a Shiny app. Setting this to \code{TRUE} only works when using this function in a Shiny app and allows for dynamic progress bars. Defaults to \code{FALSE}.
#' @param message Only used when \code{printProgress=TRUE} and \code{shiny=TRUE}. A single-element character vector: the message to be displayed to the user, or \code{NULL} to hide the current message (if any).
#' @param ... Extra parameters to be passed to the normalisation functions.
#' @return A preprocessed \code{\link[=MSnSet-class]{MSnSet}} object that is ready to be converted into a \code{\link[=protdata-class]{protdata}} object.
#' @include preprocess_hlpFunctions.R
#' @export
preprocess_wide <- function(df, accession, split, exp_annotation=NULL, type_annot=NULL, quant_cols, aggr_by, aggr_function="sum", logtransform=TRUE, base=2, normalisation="quantiles", smallestUniqueGroups=TRUE, useful_properties=NULL, filter=NULL, filter_symbol=NULL, minIdentified=2, colClasses=NA, printProgress=FALSE, shiny=FALSE, message=NULL, ...)
{

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  #If quant_cols has length one => it is a pattern => select all columns corresponding to this pattern
  if(is.character(quant_cols) & length(quant_cols==1)){quant_cols <- colnames(df)[grepl("Intensity.",colnames(df))]}

  #Error control
  if(!all(useful_properties %in% colnames(df))){stop("Argument \"useful_properties\" must only contain column names of the data frame df.")}
  if(!all(filter %in% colnames(df))){stop("One or more elements in the \"filter\" argument are no column names of the data frame df.")}

  #1. Aggregate peptides with the same sequence (parameter "aggr_by") (but maybe different charges and/or modifications)

  updateProgress(progress=progress, detail="Aggregating peptides", n=8, shiny=shiny, print=isTRUE(printProgress))

  #New progress bar for aggregation of peptides!

  progress2 <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress2 <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress2$close())
    progress2$set(message = "Aggregating peptides", value = 0)
  }


  df$MSqRob_ID <- apply(df[,c(aggr_by,filter), drop=FALSE], 1, paste , collapse = "_")

  doubleIDs <- df$MSqRob_ID[duplicated(df$MSqRob_ID)]

  if(length(doubleIDs)>0){

    classes <- lapply(df, class)
    uniqueIDs <- unique(doubleIDs)
    n <- length(uniqueIDs)

    #Everything to character
    df2 <- data.frame(lapply(df[0,], as.character), stringsAsFactors=FALSE)
    df2[n,] <- NA

    #22 minutes on our system...
    for(i in 1:n){

      updateProgress(progress=progress2, detail=paste0("Aggregating peptide ",i," of ",n,"."), n=n, shiny=shiny, print=isTRUE(printProgress))
      tmp <- df[df$MSqRob_ID==uniqueIDs[i], , drop=FALSE]
      df2[i,] <- apply(tmp, 2, function(x){paste0(unique(x), collapse = split)})
      df2[i,quant_cols] <- apply(tmp[,quant_cols, drop=FALSE],2,aggr_function)

    }

    #Set classes back as they were (if possible)
    for(j in 1:length(classes)){

      df2[,j] <- tryCatch(as(df2[,j], unlist(classes[j])), error=function(e){
        return(as.factor(df2[,j]))},
        warning=function(w){
          return(as.factor(df2[,j]))})
    }

    #rbind: "Factors have their levels expanded as necessary"
    df <- rbind(df[!(df$MSqRob_ID %in% doubleIDs), , drop=FALSE],df2)
    #Controle:
    #length(df$MSqRob_ID)==length(unique(df$MSqRob_ID)) #TRUE

  }


  intensities <- as.matrix(df[,quant_cols])
  df <- df[ , -which(names(df) %in% quant_cols)]

  #2. Log-transform

  updateProgress(progress=progress, detail="Log-transforming data", n=8, shiny=shiny, print=isTRUE(printProgress & logtransform))

  if(isTRUE(logtransform)){
    #Log transform
    intensities <- log(intensities, base=base)
  }

  #Change -Inf values in the peptide intensities to NA
  intensities[is.infinite(intensities)] <- NA

  #3,4. Normalisation

  updateProgress(progress=progress, detail="Normalizing data", n=8, shiny=shiny, print=isTRUE(printProgress & (normalisation!="none")))

  intcolnames <- colnames(intensities)

  if(normalisation=="none"){
    intensities <- intensities
  } else if (normalisation == "vsn") {

    #Remove rows that only contain NA prior to normalisation
    allNA <- apply(intensities, 1, function(x){all(is.na(x))})
    df <- df[!allNA,]
    intensities <- intensities[!allNA,]

    intensities <- exprs(vsn::vsn2(intensities, ...))
  } else if (normalisation == "quantiles") {
    intensities <- preprocessCore::normalize.quantiles(intensities, ...)
    #Reset colnames:
    colnames(intensities) <- intcolnames
  } else if (normalisation == "quantiles.robust") {
    intensities <- preprocessCore::normalize.quantiles.robust(intensities, ...)
    #Reset colnames:
    colnames(intensities) <- intcolnames
  } else if (normalisation == "center.mean") {
    center <- colMeans(intensities, na.rm = TRUE)
    intensities <- sweep(intensities, 2L, center, check.margin = FALSE, ...)
  } else if (normalisation == "center.median") {
    center <- apply(intensities, 2L, median, na.rm = TRUE)
    intensities <- sweep(intensities, 2L, center, check.margin = FALSE, ...)
  } else {
    switch(normalisation,
           max = suppressWarnings(div <- apply(intensities, 1L, max, na.rm = TRUE)), #Gives warning if no non-missing argument => -Inf => doesn't matter: NA divided by -Inf will be NA!
           sum = div <- rowSums(intensities, na.rm = TRUE))
    intensities <- intensities/div
  }

  #5. Our approach: a peptide can map to multiple proteins,
  #as long as there is none of these proteins present in a smaller subgroup

  updateProgress(progress=progress, detail="Removing overlapping protein groups", n=8, shiny=shiny, print=isTRUE(printProgress & smallestUniqueGroups))

  if(isTRUE(smallestUniqueGroups)){
    groups2 <- smallestUniqueGroups(df[,accession], split = split)
    sel <- df[,accession] %in% groups2
    df <- df[sel,]
    intensities <- intensities[sel,]
  }

  #6. Remove contaminants and reverse sequences

  updateProgress(progress=progress, detail="Removing contaminants and/or reverse sequences", n=8, shiny=shiny, print=isTRUE(printProgress & (length(filter)==0)))

  if(!is.null(filter)){
    filterdata <- df[,filter, drop=FALSE]
    #Sometimes, there are no contaminants or no reverse sequences, R then reads these empty columns as "NA"
    filterdata[is.na(filterdata)] <- ""
    sel <- rowSums(filterdata!= filter_symbol)==length(filter)
    df <- df[sel,]
    intensities <- intensities[sel,]
  }

  updateProgress(progress=progress, detail="Removing less useful properties", n=8, shiny=shiny, print=isTRUE(printProgress))

  #Retain only those properties in the fData slot that are useful (or might be useful) for our further analysis:
  #This always includes the accession (protein) as well as the peptide identifier (almost always)
  #If the accession is not present, add it
  #"quant_cols" is not necessary here, as it is part of "intensities"
  if(!(accession %in% useful_properties)){useful_properties <- c(accession,useful_properties)}
  if(!(aggr_by %in% useful_properties)){useful_properties <- c(aggr_by,useful_properties)}
  df <- df[,useful_properties]

  #7. How many times shoud a peptide be identified?
  #We require by default at least 2 identifications of a peptide sequence, as with 1 identification, the model will be perfectly confounded

  updateProgress(progress=progress, detail=paste0("Removing peptides identified less than ", minIdentified," times"), n=8, shiny=shiny, print=isTRUE(printProgress))

  keepers <- rowSums(!is.na(intensities))>=minIdentified
  df <- df[keepers,]
  intensities <- intensities[keepers,]

  #Put everything back together!!!
  df <- cbind(df, intensities)

  #8. Add experiment annotation

  updateProgress(progress=progress, detail="Adding experimental annotation", n=8, shiny=shiny, print=isTRUE(printProgress & !is.null(exp_annotation)))

  if(!is.null(exp_annotation)){

    pData <- makeAnnotation(exp_annotation, type_annot=type_annot, colClasses=colClasses)
    attr(df, "MSqRob_exp_annotation") <- exp_annotation

    # pData <- makeAnnotation(exp_annotation=exp_annotation, type_annot=type_annot, colClasses=colClasses)
    #
    # #MERGE!!! -> of weglaten: beter toevoegen bij protdata!!!!
    #
    # #Check which column of the given exp_annotation (pData) contains exactly the same elements as the mass spec run names in the data
    # annotation_run <- getAnnotationRun(pData, colnames(intensities))
    #
    # #Error checks
    # check_expAnn(pData, annotation_run)
    #
    # #Sort columns in intensities by annotation_run column
    # intensities <- intensities[,match(as.character(pData[,annotation_run]), colnames(intensities))]
    # rownames(pData) <- colnames(intensities)


  }

  return(df)
}
