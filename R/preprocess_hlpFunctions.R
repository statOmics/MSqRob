#' Initialize an annotation Excel file based on a MaxQuant peptides.txt file
#'
#' @description Creates an Excel file with one column containing the run names of an experiment based on the column names of a MaxQuant peptides.txt file. This function might be a start in constructing the first column of the annatation Excel file, other colunms need to be added manually.
#' @param file The name of a MaxQuant peptides.txt file. For more details about how this argument can be specified, see \code{\link[utils]{read.table}}.
#' @param savepath The file path to the directory where you want the output Excel annotation file to be saved. If set to \code{NULL} (the default), the Excel file will be saved in the working directory.
#' @param output_name The name of the output Excel annotation file. Defaults to \code{"experimental_annotation"}.
#' @param col_name The name of the only column created in the Excel file. Defaults to \code{"run"}.
#' @param pattern A character string containing a regular expression that will be matched to the file's header. The columns matching the expression should be the columns containing the peptide intensity values. Defaults to "Intensity.".
#' @param remove_pattern A logical indicating whether the expression in "pattern" should be removed from the column names in the resulting Excel file. Defaults to \code{TRUE}.
#' @return An Excel file with one column containing the different run names present in the MaxQuant peptides.txt file and column name \code{col_name}.
#' @export
init_ann_MQ_Excel <- function(file, savepath=NULL, output_name="experimental_annotation", col_name="run", pattern="Intensity.", remove_pattern=TRUE){

  if(is.null(savepath)){savepath <- getwd()}

  colInt <- MSnbase::grepEcols(file, pattern=pattern, split = "\t")
  runs <- read.table(file, header=FALSE, nrow=1, sep = "\t", quote = "",
                     stringsAsFactors = FALSE, comment.char = "")[colInt]
  runs <- make.names(runs, unique = TRUE)

  if(isTRUE(remove_pattern)){
    #Remove pattern from colnames of exprs
    runs <- gsub(pattern,"",runs)
  }

  run_column <- data.frame(run = runs, stringsAsFactors = TRUE)
  colnames(run_column) <- col_name

  #Remove potential extension in output_name to avoid double extension
  output_name <- gsub(".xlsx","",output_name)

  openxlsx::write.xlsx(run_column, file = file.path(savepath, paste0(output_name,".xlsx")), colNames = TRUE, rowNames = FALSE)

}

#' Create an annotation data frame
#'
#' @description Creates an annotation data frame based on either an existing annotation data frame or a path to a file which contains the experiment annotation.  Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}. As an error protection measurement, leading and trailing spaces in each column are trimmed off.
#' @param exp_annotation Either the path to the file which contains the experiment annotation or a data frame containing the experiment annotation. Exactly one colum in the experiment annotation should contain the mass spec run names.  Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}.
#' @param run_names A vector containing the names of the samples. Each name must be unique.
#' @param type_annot If \code{exp_annotation} is a path to a file, the type of file. \code{type_annot} is mostly obsolete as supported files will be automatically recognized by their extension. Currently only \code{"tab-delim"} (tab-delimited file), \code{"xlsx"} (Office Open XML Spreadsheet file) and \code{NULL} (file type decided based on the extension) are supported. If the extension is not recognized, the file will be assumed to be a tab-delimited file. Defaults to \code{NULL}.
#' @param colClasses character. Only used when the \code{exp_annotation} argument is a filepath. A vector of classes to be assumed for the columns. Recycled if necessary. If named and shorter than required, names are matched to the column names with unspecified values are taken to be NA.
#' Possible values are \code{"keep"} (the default, when the colClasses are unchanged for data frames and \code{type.convert} is used for files),  \code{NA} (when \code{type.convert} is always used), \code{NULL} (when the column is skipped), one of the atomic vector classes (\code{"logical"}, \code{"integer"}, \code{"numeric"}, \code{"complex"}, \code{"character"}, \code{"raw"}), or \code{"factor"}, \code{"Date"} or \code{"POSIXct"}. Otherwise there needs to be an as method (from package \code{methods}) for conversion from \code{"character"} to the specified formal class.
#' @return A data frame containing the experimental annotation. Possible leading and trailing spaces are trimmed.
#' @export
makeAnnotation <- function(exp_annotation, run_names, type_annot=NULL, colClasses="keep"){

  #If exp_annotation is not a data frame, it should be imported from a file.
  if(!is.data.frame(exp_annotation)){

    if(isTRUE(as.logical(grep(".xlsx[/\\]*$",exp_annotation))) || isTRUE(type_annot=="xlsx")){
      ann_frame <- openxlsx::read.xlsx(exp_annotation)
      #Convert characters to factors: this is the default for files and will be put into your colClasses if colClasses="keep"!
      ann_frame <- as.data.frame(unclass(ann_frame))
    } else{ #if(type_annot=="tab-delim")
      ann_frame <- read.table(exp_annotation, sep="\t", header=TRUE, row.names = NULL, stringsAsFactors = TRUE)
    }

    #If exp_annotation is a data frame, just put the data frame in the pData slot (lapply and gsub are just to remove leading and trailing spaces)
  } else{
    ann_frame <- exp_annotation
  }

  if(!is.na(colClasses) && colClasses=="keep"){colClasses <- vapply(ann_frame,class,"")}

  #First: everything to factor, gsub is to remove leading and trailing spaces
  pData_makeAllNames <- make_pData_MSqRob(ann_frame)

  #Check which column of the given exp_annotation (pData) contains exactly the same elements as the mass spec run names in the data
  annotation_run <- getAnnotationRun(pData=pData_makeAllNames, run_names=run_names)

  #Now only "make.names" on the annotation column
  pData <- make_pData_MSqRob(ann_frame, annotation_run=annotation_run)

  pData <- addColClasses(pData, colClasses)

  #Error checks
  check_expAnn(pData, annotation_run)

  return(pData)
}


#MSqRob "make names"
#If annotation_run==NULL => "make.names" on everything, else only on the annotation_run
make_pData_MSqRob <- function(pDataObject, annotation_run=NULL){

  ann_run_index <- which(colnames(pDataObject)==annotation_run)
  if(length(ann_run_index)==0){ann_run_index <- NULL}

  counter <- 0

  pDataObject <- data.frame(lapply(pDataObject, function (x) {

    counter <<- counter + 1

    #Remove leading and trailing spaces for every column
    x <- gsub("^\\s+|\\s+$", "", x)

    #If it's the unique column, apply make.names so that its names are equal to the colnames of the exprs slot of the MSnSet object (read in via read.table)
    if(counter==ann_run_index || is.null(ann_run_index)){
      x <- make.names(x, unique = TRUE)}

    return(factor(x, levels=unique(x))) #levels are specified to give factor levels the order given in the annotation file
  }
  ))
  return(pDataObject)
}



#addColClasses

addColClasses <- function(df, colClasses){

  #If the type is specified by a named vector
  names <- names(colClasses)
  if(!is.null(names)){
    for(i in 1:length(colClasses)){
      #"as" doesn't work for factors => need tryCatch
      df[,names[i]] <- tryCatch(as(df[,names[i]], colClasses[i]), error=function(e){
        return(df[,names[i]])})

      #Else, the type is specified by the order
    }} else if(!all(is.na(colClasses))){
      colClasses <- rep(colClasses, length=ncol(df))
      for(i in 1:ncol(df)){
        #"as" doesn't work for factors => need tryCatch
        df[,i] <- tryCatch(as(df[,i], colClasses[i]), error=function(e){
          return(df[,i])})

      #Else, colClasses is NA, so we return to the default behaviour of "type.convert"
      }} else{
        for(i in 1:ncol(df)){
          #"as" doesn't work for factors => need tryCatch
          df[,i] <- type.convert(df[,i])
        }}

  #The real function return:
  return(df)
}

#' Get the column with the mass spec run names
#'
#' @description Get the name of the column of an experimental annotation data frame \code{pData}. This column should contain exactly the same elements as the elements in the \code{run_names} vector.
#' This function is borderline internal.
#' @param pData A dataframe containing experimental annotation information, of which one column contains the mass spec run names.
#' @param run_names A vector with the mass spec run names.
#' @return The name of the column in \code{pData} that contains all elements in the vector \code{run_names}.
#' @export
getAnnotationRun <- function(pData, run_names){
  annotation_run <- names(which(vapply(pData, function(x) return(identical(sort(as.character(run_names)),sort(as.character(x)))), FUN.VALUE = TRUE)))
  return(annotation_run)
}


#check_expAnn

check_expAnn <- function(pData, annotation_run){
  #Error checks
  if(length(annotation_run)>1){stop(paste0("Your experiment annotation has multiple columns with elements equal to the mass spec run names in the data. Please remove one of the following columns: ",paste0(annotation_run, collapse=", "),"."))
  } else if(length(annotation_run)==0){
    uniquelvls <- names(which(apply(pData,2, function(x) length(unique(x))==nrow(pData))))
    if(length(uniquelvls)>0){stop(paste0("Make sure that exactly one column in the experiment annotation has elements equal to the mass spec run names in the data. Maybe one of the following columns has an error in at least one of its elements: ",paste0(uniquelvls, collapse=", "),"."))
    } else{stop("Make sure that exactly one column in the experiment annotation has elements equal to the mass spec run names in the data.")}
  } else if(!is.factor(pData[,annotation_run])){
    stop("Please make sure that the column in your experiment annotation that corresponds to the mass spec run names is a factor variable.")
  }
}

#' Add annotations to a data frame in "long" format
#'
#' @description This function adds extra columns to a data frame in "long" format (i.e. the data frame has one observation row per measurement (thus, multiple rows per subject)) based on either a path to an annotation file or an annotation data frame.
#' @param df A data frame in "long" format to which extra annotation columns need to be added.
#' @param run_col A character indicating the column in the data frame \code{df} that contains a different identifier for each mass spec run. All elements in this column should be present exactly once in exactly one column of the experimental annotation data frame \code{pData} and be equal to the elements of the \code{run_col_annot} column in the \code{pData} data frame.
#' @param run_col_annot A character indicating the column in the annotation data frame \code{pData} that contains a different identifier for each mass spec run. All elements in this column should be present exactly once in this column of the experimental annotation data frame \code{pData} and be equal to the elements of the \code{run_col} column in the data frame \code{df}.
#' @param pData A data frame containing the experiment annotation. Exactly one colum in the experiment annotation should contain the mass spec run names.
#' @return A data frame to which the annotations are left-joined.
#' @export
pasteAnnotation <- function(df, run_col, run_col_annot, pData){

  #Merge in experimental annotation: almost instantly done!
  df <- merge(df, pData, by.x=run_col, by.y=run_col_annot, all.x=TRUE, all.y=FALSE)

  #Make mass spec run column the last column (so that intensity remains the first column)
  df <- data.frame(df[,-1,drop=FALSE],df[,1,drop=FALSE])

  return(df)
}



