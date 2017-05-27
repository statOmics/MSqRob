#' Import data independent acquisition files
#'
#' @description This function generates a peptide data frame from output files originating from either DIAumpire, OpenSWATH, PeakView, Skyline or Spectronaut.
#' It is largely based on the \code{\link[LFQbench]{FSWE.generateReports}} function from the \code{LFQbench} package and heavily depends on this package.
#' It is also possible to define a new software format using the \code{\link[LFQbench]{FSWE.addSoftwareConfiguration}} function from the \code{LFQbench} package (Navarro et al., 2016).
#' @param experimentFile The name of a file. For more details about how this argument can be specified, see \code{\link[utils]{read.table}}.
#' @param aggr_by A character indicating the column by which the data should be aggregated. The default \code{"sequence.var"} will aggregate the data over different charge states and modification statuses. If you only want to aggregate over charge states, set \code{aggr_by} to \code{"sequence.mod.var"}. If no aggregation at all is desired, set \code{aggr_by} to \code{"none"}. Data will never be aggregated over different \code{filename.var}.
#' @param softwareSource The default, \code{"guess"}, tries to guess from which software the file originates.
#' @param remove_decoys Should decoys be removed? Defaults to \code{FALSE} as removing of decoys is typically done during our preprocessing step. Can be set to \code{TRUE} to remove decoys already upfront.
#' @references Navarro P., Kuharev J., Gillet, L. C., Bernhardt, O. M., MacLean, B., R\"ost, H. L., Tate, S. A., Tsou, C., Reiter, L., Distler, U., Rosenberger, G., Perez-Riverol, Y., Nesvizhskii, A. I., Aebersold, R., and Tenzer, S. (2016) A multicenter study benchmarks software tools for label-free proteome quantification. Nature Biotechnology.
#' @export
importDIAData <- function (experimentFile, aggr_by = "sequence.var", softwareSource = "guess", remove_decoys = FALSE)
{

  if(!(aggr_by %in% c("none","sequence.var","sequence.mod.var"))){stop("aggr_by should be one of \"none\", \"sequence.var\" or \"sequence.mod.var\".")}

  #Test if experimentFile is an URL (if URL: read.xlsx doesn't work => download ourselves)
  con.url <- try(url(experimentFile, open='rb'))
  try.error <- inherits(con.url, "try-error") #If false: URL, if TRUE, no url

  #If URL: read.xlsx doesn't work => download ourselves
  if(!isTRUE(try.error)){
  tmp = paste0(tempdir(),"/",basename(experimentFile))
  download.file(url = experimentFile, destfile = tmp, mode="wb")
  experimentFile=tmp
  }

  LFQbench:::FSWE.init.defaultSoftwares()
  LFQbench:::FSWE.init.defaultModifications()
  LFQbench:::LFQbench.initConfiguration()

  ###dataSets = as.list(FSWE.dataSets)

  #LFQbench.setDataRootFolder(working_dir, createSubfolders = FALSE)
  if (softwareSource == "guess") {
    softwareSource <- LFQbench:::guessSoftwareSource(basename(experimentFile),
                                          FSWE.softwareNames)
  }

  #Make variables like qvalue.var, sequence.mod.var,... depending on the software source
  #and put them in the FSWE.Config object + error control: is softwareSource one of the defined software sources.
  LFQbench:::FSWE.switchSoftwareConfiguration(softwareSource)

  #Assign the elements of the FSWE.Config object we just made to the current environment
  nix = sapply(names(FSWE.Config), function(pn) assign(pn,
                                                       FSWE.Config[[pn]], envir = parent.env(environment())))

  sumquant <- paste0("sum(", quantitative.var, ")")
  medianquant <- paste0("median(", quantitative.var, ")")
  qvalue.filtered = FALSE
  #original.experimentFile <- experimentFile
  #experimentFile <- file.path(working_dir, experimentFile)
  cat(paste0("Preprocessing ", experimentFile,
             "\n"))

  #Import data
  if (grepl(".xls", input.extension)) {
    df <- readxl::read_excel(experimentFile, sheet = sheet.data,
                     col_names = TRUE)
    # df <- openxlsx::read.xlsx(experimentFile, sheet = sheet.data,
    #                           colNames = TRUE)
    names(df) <- gsub(" ", ".", names(df))
    if (!is.na(q_filter_threshold)) {
      df$sequence_z <- paste(df[[sequence.mod.var]], df[[charge.var]],
                             sep = "_")
      df$sequence_z <- as.character(df$sequence_z)
      df.fdr <- readxl::read_excel(experimentFile, sheet = sheet.fdr,
                           col_names = TRUE)
      names(df.fdr) <- gsub(" ", ".", names(df.fdr))
      df.fdr$sequence_z <- paste(df.fdr[[sequence.mod.var]],
                                 df.fdr[[charge.var]], sep = "_")
      df <- df %>% filter(!duplicated(sequence_z))
      df.fdr <- df.fdr %>% filter(!duplicated(sequence_z))
      tmp1 <- data.frame(sequence_z = as.character(df.fdr$sequence_z))
      tmp1$sequence_z <- as.character(tmp1$sequence_z)
      valid_peptides <- df.fdr[, grepl(fdr.var.tag, colnames(df.fdr),
                                       ignore.case = TRUE)] <= q_filter_threshold
      tmp1 <- cbind(tmp1, valid_peptides)
      df.fdr <- tmp1
      rm(tmp1)
      df2 <- left_join(df, df.fdr, by = "sequence_z")
      df2 <- unique(df2)
      df2.flags <- df2[, grepl(fdr.var.tag, colnames(df2),
                               ignore.case = TRUE)]
      df2.values <- df2[, grepl(quantitative.var.tag, colnames(df2),
                                ignore.case = TRUE) & !grepl(fdr.var.tag, colnames(df2),
                                                          ignore.case = TRUE)]
      df2.values[!df2.flags] <- NA
      rm(df2)
      rm(df2.flags)
      rm(valid_peptides)
      rm(df.fdr)
      df <- df[, !grepl(quantitative.var.tag, colnames(df),
                        ignore.case = TRUE)]
      df <- cbind(df, df2.values)
      qvalue.filtered = TRUE
    }
  }
  else {
    if (LFQbench:::guessSep(experimentFile) == ",") {
      df <- read_csv(experimentFile, na = nastrings)
    }
    if (LFQbench:::guessSep(experimentFile) == "\t") {
      df <- read_tsv(experimentFile, na = nastrings)
    }
  }

  names(df) <- make.names(names(df))
  theEnvir = environment()
  matchColumnNames <- function(paramName = "protein.var") {
    paramValues = get0(paramName)
    if (is.null(paramValues) || is.na(paramValues))
      return(F)
    if (length(paramValues) > 1) {
      columnNames = names(df)
      matchingValues = paramValues %in% columnNames
      if (!any(matchingValues)) {
        stop(paste0("software configuration does not match the input file!\n",
                    "file: ", experimentFile, "\n", paramName,
                    ": ", paste(paramValues, collapse = ", "),
                    "\n", "column names: ", paste(columnNames,
                                                  collapse = ", "), "\n"))
      }
      assign(paramName, paramValues[matchingValues][1],
             envir = theEnvir)
    }
  }

  #Add sequence var using own function setSequenceVar
  df_FSWE.Config <- setSequenceVar(df, FSWE.Config, softwareSource)
  df <- df_FSWE.Config[["df"]]
  FSWE.Config <- df_FSWE.Config[["FSWE.Config"]]

  #Again assign the elements of the FSWE.Config object to the current environment (so that we also have sequence.var)
  nix = sapply(names(FSWE.Config), function(pn) assign(pn,
                                                       FSWE.Config[[pn]], envir = parent.env(environment())))

  if (!is.na(q_filter_threshold) & !qvalue.filtered) {
    df <- eval(substitute(filter(df, var < q_filter_threshold),
                          list(var = as.name(qvalue.var))))
    qvalue.filtered = TRUE
  }

  #Remove decoys
  if(isTRUE(remove_decoys)){
  for (dectag in decoy.tags) {
    df <- df[!grepl(dectag, df[[protein.var]], ignore.case = TRUE),]
  }
  if (!is.na(decoy.var)) {
    is_decoytype_int = typeof(as.matrix(df[, decoy.var])) ==
      "integer"
    if (is_decoytype_int) {
      df <- df[tolower(df[[decoy.var]]) == 0, ]
    }
    else {
      df <- df[tolower(df[[decoy.var]]) == "false", ]
    }
  }
  }

  df <- df %>% rowwise()

  #Remove empty and multiple species
  #df <- filter(df, species != "NA", species != "multiple")
  experiment <- NA
  if (input_format == "wide") {
    # experiment <- which(sapply(dataSets, LFQbench:::guessExperiment_wide,
    #                            colnames(df)))
    tmp1 <- df[, protein.var]
    tmp1 <- cbind(tmp1, df[, sequence.var])
    tmp1 <- cbind(tmp1, df[, sequence.mod.var])
    tmp1 <- cbind(tmp1, df[, charge.var])
    quant.columns <- grepl(quantitative.var.tag, colnames(df),
                           ignore.case = TRUE)
    # injection.columns <- Reduce(`|`, lapply(dataSets[[experiment]],
    #                                         grepl, colnames(df)))
    tmp1 <- cbind(tmp1, df[, quant.columns]) # & injection.columns
    df <- tmp1
    rm(tmp1)
    quantvar.range <- which(grepl(quantitative.var.tag, colnames(df),
                                  ignore.case = TRUE))
    quantvar.colnames <- names(df)[quantvar.range]
    if (is.na(filename.var))
      filename.var <- "filename.var"
    df <- df %>% gather_(filename.var, quantitative.var,
                         quantvar.colnames) %>% arrange_(protein.var, get(aggr_by))
  }
  else if (input_format == "long") {
    df[[filename.var]] <- basename(file_path_sans_ext(df[[filename.var]]))
    df[[filename.var]] <- basename(file_path_sans_ext(df[[filename.var]]))
  }
  data <- df %>% distinct_(filename.var, sequence.var, sequence.mod.var,
                           charge.var, .keep_all = TRUE)
  data[[quantitative.var]][data[[quantitative.var]] == 0] <- NA
  if(aggr_by!="none"){
  data <- data %>% ungroup() %>% group_by_(filename.var, get(aggr_by),
                                           protein.var) %>% summarise_(quant_value = sumquant)
  }
  return(as.data.frame(data))
}



setSequenceVar <- function(df, FSWE.Config, softwareSource){

  #PViewBuiltinProteins doesn't have a sequence.mod.var nor a sequence.var => cannot be used!

  index <- which(c("Spectronaut", "DIAumpire", "Skyline", "PeakView","PViewNoFilter", "OpenSWATH", "DIAumpBuiltinProteins")==softwareSource)

  #For PeakView and PViewNoFilter, add a new colunm Sequence to data frame df
  if(index %in% c(4,5)){df$Sequence <- gsub("\\[.*\\]", "", df[[FSWE.Config$sequence.mod.var]])}

  sequence.vars <- c("EG.StrippedSequence", "Sequence", "PeptideSequence", gsub("\\[.*\\]", "", "Sequence"),
                    "Sequence", "Sequence", "Sequence")

  FSWE.Config$sequence.var <- sequence.vars[index]

  df_FSWE.Config <- list(df=df, FSWE.Config=FSWE.Config)

  return(df_FSWE.Config)
}

