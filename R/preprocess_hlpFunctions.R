#' Create an annotation data frame
#'
#' @description Creates an annotation data frame based on either an existing annotation data frame or a path to a file which contains the experiment annotation.  Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}.  As an error protection measurement, leading and trailing spaces in each column are trimmed off.
#' @param exp_annotation Either the path to the file which contains the experiment annotation or a data frame containing the experiment annotation. Exactly one colum in the experiment annotation should contain the mass spec run names.  Annotation in a file can be both a tab-delimited text document or an Excel file. For more details, see \code{\link[utils]{read.table}} and \code{\link[openxlsx]{read.xlsx}}.
#' @param type_annot If \code{exp_annotation} is a path to a file, the type of file. \code{type_annot} is mostly obsolete as supported files will be automatically recognized by their extension. Currently only \code{"tab-delim"} (tab-delimited file), \code{"xlsx"} (Office Open XML Spreadsheet file) and \code{NULL} (file type decided based on the extension) are supported. If the extension is not recognized, the file will be assumed to be a tab-delimited file. Defaults to \code{NULL}.
#' @param colClasses character. Only used when the \code{exp_annotation} argument is a filepath. A vector of classes to be assumed for the columns. Recycled if necessary. If named and shorter than required, names are matched to the column names with unspecified values are taken to be NA.
#' Possible values are \code{NA} (the default, when \code{type.convert} is used), \code{NULL} (when the column is skipped), one of the atomic vector classes (\code{logical}, \code{integer}, \code{numeric}, \code{complex}, \code{character}, \code{raw}), or \code{factor}, \code{Date} or \code{POSIXct}. Otherwise there needs to be an as method (from package \code{methods}) for conversion from \code{character} to the specified formal class.
#' @return A data frame containing the experimental annotation. Possible leading and trailing spaces are trimmed.
#' @export
makeAnnotation <- function(exp_annotation, type_annot=NULL, colClasses=NA){

  #If exp_annotation is not a data frame, it should be imported from a file.
  if(!is.data.frame(exp_annotation)){

    if(isTRUE(as.logical(grep(".xlsx[/\\]*$",exp_annotation))) || type_annot=="xlsx"){

      if(Sys.info()['sysname']=="Windows"){

        pData <- xlsx::read.xlsx(file=exp_annotation, sheetIndex = 1)
        #Remove columns with all NA
        pData <- pData[,!colSums(apply(pData, 2, is.na))==nrow(pData),drop=FALSE]
        return(pData)

      } else{pData <- openxlsx::read.xlsx(exp_annotation)}

    } else{ #if(type_annot=="tab-delim")
      pData <- read.table(exp_annotation, sep="\t", header=TRUE, row.names = NULL)
    }

    #Default: everything to factor, gsub is to remove leading and trailing spaces
    pData <- data.frame(lapply(pData, function(x){return(factor(gsub("^\\s+|\\s+$", "", x), levels=unique(gsub("^\\s+|\\s+$", "", x))))})) #levels are specified to give factor levels the order given in the annotation file
    #pData <- data.frame(lapply(pData, function(x){return(as.factor(x))}))

    pData <- addColClasses(pData, colClasses)

    #If exp_annotation is a data frame, just put the data frame in the pData slot (lapply and gsub are just to remove leading and trailing spaces)
  } else{pData <- data.frame(lapply(exp_annotation, function (x) factor(gsub("^\\s+|\\s+$", "", x))))}

  return(pData)
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
      colClasses <- rep(colClasses, length=length(colClasses))
      for(i in 1:length(colClasses)){
        #"as" doesn't work for factors => need tryCatch
        df[,i] <- tryCatch(as(df[,i], colClasses[i]), error=function(e){
          return(df[,i])})
      }}

  #The real function return:
  return(df)
}

#getAnnotationRun

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



