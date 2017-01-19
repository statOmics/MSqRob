#' Data frame to protein data
#'
#' @description Converts data frame \code{df} to a a \code{\link[=protdata-class]{protdata}} object, which can be used for further analysis.
#' @param df A data frame that needs to be converted to a \code{\link[=protdata-class]{protdata}} object.
#' @param acc_col A character string or numeric index indicating the column in data frame \code{df} that contains the identifiers by which the data should be grouped in the resulting \code{\link[=protdata-class]{protdata}} object.
#' @param quant_cols A vector of character strings or numeric indices indicating which column(s) in the data frame \code{df} contain the quantitative values of interest (mostly peptide intensities or peptide areas under the curve). If \code{quant_cols} contains only one element, \code{df2protdata} will assume that data frame \code{df} is in "long" format. If \code{quant_cols} contains more than one element, \code{df2 protdata} will assume the data to be in "wide" format.
#' @param quant_name A character string indicating the name that will be given to the column that will contain the quantitative values of interest (mostly peptide intensities or peptide areas under the curve). Defaults to \code{"quant_value"}.
#' @param run_name If quant_cols contains more than one element (i.e. the data is in "wide" format), this should contain a freely chosen character string indicating the name that will be given to the column containing the mass spec run names. If no name is chosen (\code{run_name=NULL}), this will default to \code{"run"}.
#' If quant_cols contains only one element (i.e. the data is in "long" format), \code{run_name} should contain the name of the column that contains the run names.
#' @param annotations A vector of character strings or numeric indices indicating the columns in the data frame \code{df} that contain additional information on the accessions (typically protein names, gene names, gene ontologies,...) that should be added in a separate \code{annotation} slot. In case multiple values of the same annotation column would exist for a unique accession, these values are pasted together. Defaults to \code{NULL}, in which case no annotations will be added.
#' @return A \code{\link[=protdata-class]{protdata}} object.
#' @examples #This example will convert df object peptides into a protdata object proteins.
#' #Import the data as a df object
#' pepdf <- read.table(system.file("extdata/CPTAC", "peptides.txt", package = "MSqRob"), sep="\t", header=TRUE)
#' #To save time, we only take the first 50 peptides as an example
#' pepdf <- pepdf[1:50,]
#' #Determine columns that contain the intensity values
#' quant_cols <- colnames(pepdf)[which(grepl("Intensity.",colnames(pepdf)))]
#' #Log2 transform data and change -Inf to NA
#' pepdf[,quant_cols] <- log2(pepdf[,quant_cols])
#' tmp_ints <- pepdf[,quant_cols]
#' tmp_ints <- as.matrix(tmp_ints)
#' tmp_ints[is.infinite(tmp_ints)] <- NA
#' tmp_ints <- as.data.frame(tmp_ints)
#' pepdf[,quant_cols] <- tmp_ints
#' #Keep only columns of interest
#' pepdf <- pepdf[,c(quant_cols,"Proteins","Sequence","PEP")]
#' #Dermine the column that contains the protein names
#' acc_col <- "Proteins"
#' proteins <- df2protdata(pepdf, acc_col, quant_cols)
#' @include protdata.R
#' @export
df2protdata <- function(df, acc_col, quant_cols, quant_name="quant_value", run_name=NULL, annotations=NULL){

  if(is.null(run_name) & length(quant_cols)==1){stop("Please provide the name of the column containing the mass spec runs as input in run_name.")}
  if(is.null(run_name) & length(quant_cols)>1){run_name <- "run"}

  #If there is an experimental annotation attribute, store it in data frame pData
  pData <- attr(df, "MSqRob_exp_annotation")

  #Checken of functie de output geeft die ze zou moeten geven, dat er geen dingen omgewisseld zijn enzo
  #Toevoegen: if datatype invoer niet correct
  df <- droplevels(df)
  proteins <- factor(unique(df[,acc_col]))
  datalist <- replicate(length(proteins), data.frame())
  annotation_matrix <- matrix(nrow=length(proteins), ncol=length(annotations), dimnames=list(proteins,annotations))

  #Automatically determine the column name in the pData that corresponds to run_name in the df

  if(length(quant_cols)==1){
  run_names <- unique(df[,run_name])
  } else{
    run_names <- colnames(df[,quant_cols])
    #If there are no column names, we make them ourselves:
    if(is.null(run_names)){run_names <- paste0("run_",1:length(quant_cols))}
  }

  #If pData is NULL, turn it into a data frame with one a column for run
  if(is.null(pData)){pData <- data.frame(run=run_names)
  colnames(pData) <- run_name}

  #Check which column of the given exp_annotation (pData) contains exactly the same elements as the mass spec run names in the data
  run_col_annot <- getAnnotationRun(pData=pData, run_names=run_names)

  #Error checks
  check_expAnn(pData=pData, annotation_run=run_col_annot)

  for(i in 1:length(proteins)){

    #Select the rows that correspond to protein i
    sel <- df[, acc_col] == proteins[i]
    #Select the intensities based on the user input
    intensities <- df[sel, quant_cols, drop=FALSE]
    #All that is not intensities nor proteins is properties
    properties <- df[sel,, drop=FALSE]
    annotation_matrix[i,] <- apply(properties[,annotations, drop=FALSE], 2, function(x){paste0(unique(x))})
    properties <- properties[,-which(colnames(properties) %in% c(colnames(df[,acc_col,drop=FALSE]),colnames(df[,annotations,drop=FALSE]),colnames(df[,quant_cols,drop=FALSE]))), drop=FALSE]

    #If there are multiple intensity columns, put them under each other.
    value <- as.vector(as.matrix(intensities))
    framelist <- replicate(ncol(intensities), properties, simplify = FALSE)
    frame <- plyr::rbind.fill(framelist)

    #If wide format (more than one intensity column)
    if(length(quant_cols)>1){
    #Add intensity values in the first column and run names in the last column
    frame2 <- data.frame(value, frame, rep(colnames(intensities), each=nrow(intensities)))
    #Name the newly created run column.
    colnames(frame2)[ncol(frame2)] <- run_name

    #If long format (only one intensity column)
    } else{
      frame2 <- data.frame(value, frame)
    }

    frame2 <- pasteAnnotation(df=frame2, run_col=run_name, run_col_annot=run_col_annot, pData=pData)
    #Name the newly create intensity column.
    colnames(frame2)[1] <- quant_name

    #Removing NA
    frame3 <- frame2[complete.cases(frame2[,1]),, drop=FALSE]

    datalist[[i]] <- frame3
  }

  protdata <- new("protdata", accession=proteins, data=datalist, annotation=annotation_matrix)
  return(protdata)
}
