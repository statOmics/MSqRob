#' Save the variables of a data.frame in distinct binary files
#'
#' @description \code{saves_MSqRob} is almost a pure copy of the \code{saves} function from the \code{saves} package by Dar\'oczi (2013) with some minor code tweaks to make it work for MSqRob. It saves dataframe(s) or list(s) to disk in a special, binary format. This binary format consists of distinct binary files of all separate variables of a dataframe/list merged into an uncompressed tar archive. This is done via a loop, which saves each variable/column to an an external representation of the R objects via save in a temporary directory. Theese 'RData' files are archived to an 'RDatas' tar archive, uncompressed for better speed.
#' @param ... R objects: the names of the objects to be saved (as symbols or character strings)
#' @param list character vector: the name(s) of the data frame(s) or list(s) to save
#' @param file character vector: the (RDatas) filename(s) in which to save the variables in the current working directory
#' @param overwrite boolean: if TRUE, existing files will be deleted before saving. Default set to FALSE, which will report error on conflicting file names.
#' @param ultra.fast boolean: if TRUE, ultra fast (...) processing is done without any check to parameters, also no archiving or compression is done. Be sure if using this setting, as many uncompressed files could be generated in the working directory's subdirectory named to df. Only recommended for servers dealing with lot of R objects' saves and loads in a monitored environment.
#' @return The saved filename(s) (invisible).
#' @family \code{\link{loads_MSqRob}} to load R objects from RDatas binary format.
#' @family \code{\link{inspect_loads_MSqRob}} to inspect the content of an RDatas binary object.
#' @examples ## Not run:
#' ## Saving the demo dataset to evs.2000.hun.RDatas in current working directory.
#' data(evs.2000.hun)
#' saves_MSqRob(evs.2000.hun)
#' ## Saving both the demo dataset and mtcars to current working directory
#' saves_MSqRob(evs.2000.hun, mtcars)
#' saves_MSqRob(list=c('evs.2000.hun', 'mtcars'))
#' ## Saving all kind of cars :)
#' saves_MSqRob(cars, mtcars, overwrite = TRUE)
#' saves_MSqRob(list=c('cars', 'mtcars'), overwrite = TRUE)
#'
#' ## End(Not run)
#' @references Dar\'oczi, G. (2013). saves: Fast load variables. R package version 0.5, URL http://cran.r-project.org/package=saves
#' @export
saves_MSqRob <- function (..., envir = environment(), list = character(), file = NULL, overwrite = FALSE,
                          ultra.fast = FALSE, shiny=FALSE, printProgress=FALSE, message=NULL)
{

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  names <- as.character(substitute(list(...)))[-1L]
  list <- c(list, names)
  if (ultra.fast == TRUE) {
    df <- list[1]
    data <- get(df, envir=envir)
    dir.create(df)
    e <- as.environment(data)
    for (i in 1:length(data)) {
      save(list = names(data)[i], file = paste(df, "/",
                                               names(data)[i], ".RData", sep = ""), compress = FALSE,
           precheck = FALSE, envir = e)
    }
    return(invisible(df))
  }
  if (is.null(file))
    file <- paste(list, ".RDatas", sep = "")
  if (length(list) != length(file))
    stop("Bad number of files given!")
  for (i in 1:length(list)) {
    if (inherits(try(data <- get(list[i], envir=envir), silent = TRUE),
                 "try-error"))
      stop(paste("No dataframe/list given or `", list[i],
                 "` is not a dataframe/list!"))
    if (file.exists(file[i])) {
      if (overwrite == TRUE) {
        file.remove(file[i])
      }
      else {
        stop(paste("Destination filename `", file[i],
                   "` already exists! Use other filename or use paramater `overwrite` set to TRUE."))
      }
    }
    if (!is.data.frame(data) & !is.list(data))
      stop(paste("No dataframe/list given or `", list[i],
                 "` is not a dataframe/list!"))
    tmp <- tempfile("saves.dir-")
    dir.create(tmp)
    e <- as.environment(data)

    count <- 0
    lapply(names(data), function(x) {
      count <<- count+1
      updateProgress(progress=progress, detail=paste0("Saving ",x,"."), n=length(data), shiny=shiny, print=isTRUE(printProgress))
      save(list = x, file = paste(tmp,"/", x, ".RData", sep = ""), envir = e)
    })

    w <- getwd()

    #Added to make it possible to save to other directories
    if(dirname(file[i])=="."){savedir <- w
    } else{savedir <- dirname(file[i])}

    setwd(tmp)
    tar(paste(savedir, "/", basename(file[i]), sep = ""), ".", compression = "none")

    setwd(w)
    unlink(tmp, recursive = TRUE)
  }
  invisible(file)
}

#' Loading only given variables of a data.frame from binary file
#'
#' @description \code{loads_MSqRob} is almost a pure copy of the \code{loads} function from the \code{saves} package by Dar\'oczi (2013) with some minor code tweaks to make it work for MSqRob. It loads data from a special binary file format (RDatas) made up by the \code{\link{MSqRob_saves}} function. This special, uncompressed tar archive inlcudes several separate RData files (saved by \code{\link{MSqRob_saves}} function) as being columns/variables of a data frame.
#' @param file character string: the (RDatas) filename from which to load the variables. If using \code{ultra.fast = TRUE} option, specify the directory holding the uncompressed R objects (saved via \code{MSqRob_saves(..., ultra.fast = TRUE))}.
#' @param variables Optional: a character vector containing the variable names to load. If not specified, all variables will be loaded.
#' @param to.data.frame boolean: the default behavior of loads is to concatenate the variables to a list. This could be overriden with \code{TRUE} argument specified at to.data.frame parameter, which will return a dataframe instead of list. Only do this if all your variables have the same number of cases!
#' @param ultra.fast boolean: if \code{TRUE}, ultra fast (...) processing is done without any check to parameters or file existence/permissions. Be sure if using this setting as no debugging is done! Only recommended for servers dealing with lot of R objects' saves and loads in a monitored environment. Also, for performance gain, it is advised not to convert the list to data frame (\code{to.data.frame = FALSE}).
#' @details The purpose of this function is to be able only a few variables of a data.frame really fast. It is done by reading and writing datas in binary format without any transformations, and combining the speed of only reading the needed part of an archive.
#'
#' Some minor experiments shows a huge performance gain against using SQLite/MySQL backends or loading whole binary data, but be conscious always choosing the aprropriate method to write and read data.
#'
#' The author of the \code{saves} package (Dar\'oczi) emphasizes: this package could be useful only in few cases!
#' @return Loaded data.frame.
#' @family \code{\link{saves_MSqRob}} to save R objects to RDatas binary format
#' @family \code{\link{inspect_loads_MSqRob}} to inspect the content of an RDatas binary object.
#' @examples ## Not run:
#' # Loading the 'v1' and 'v5' variables of the demo dataset.
#' data(evs.2000.hun)
#' saves(evs.2000.hun)
#' evs.filtered.list <- loads("evs.2000.hun.RDatas", c('v1', 'v5'))
#' evs.filtered.df <- loads("evs.2000.hun.RDatas", c('v1', 'v5'), to.data.frame=TRUE)
#'
#' ## End(Not run)
#' @references Dar\'oczi, G. (2013). saves: Fast load variables. R package version 0.5, URL http://cran.r-project.org/package=saves
#' @export
loads_MSqRob <- function (file = NULL, variables = NULL,
                          ultra.fast = FALSE, printProgress=FALSE, shiny=FALSE, message=NULL)
{

  progress <- NULL
  if(isTRUE(shiny)){
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = message, value = 0)
  }

  tmp <- tempfile("saves.dir-")
  dir.create(tmp)
  untar(file, exdir = tmp)
  if(is.null(variables)){variables <- gsub(".RData","",dir(tmp))}

  data <- list()
  if (ultra.fast == TRUE) {
    for (i in 1:length(variables)) {
      updateProgress(progress=progress, detail=paste0("Loading object ", i," of ",length(variables),": ",variables[i],"."), n=length(variables), shiny=shiny, print=isTRUE(printProgress))
      f <- paste(file, "/", variables[i], ".RData", sep = "")
      data[[paste(variables[i])]] <- local(get(load(f)))
    }
    names(data) <- variables
    return(data)
  }
  if (is.null(file)) {
    stop("Argument missing! Specify a filename to load.")
  }
  if (!file.exists(file)) {
    stop("Archive not found!")
  }

  for (i in 1:length(variables)) {
    updateProgress(progress=progress, detail=paste0("Loading object ", i," of ",length(variables),": ",variables[i],"."), n=length(variables), shiny=shiny, print=isTRUE(printProgress))
    f <- paste(tmp, "/", variables[i], ".RData", sep = "")
    if (!file.exists(f)) {
      stop(paste("Variable: <<", variables[i], ">> not found!"))
    }
    data[[paste(variables[i])]] <- local(get(load(f)))
  }
  #names(data) <- variables
  unlink(tmp, recursive = TRUE)
  return(data)
}

#'Inspect the variables present in a data.frame from binary file
#'
#' @description \code{inspect_loads_MSqRob} returns the names of the variables present in a .RDatas file, some or all of which can be used in the \code{variables} argument of the \code{\link{loads_MSqRob}} function if one prefers to load only some variables.
#' @param file character string: the (RDatas) filename from which to inspect the variable names.
#' @return The variable names present in the .RDatas file.
#' @family \code{\link{saves_MSqRob}} to save R objects to RDatas binary format
#' @family \code{\link{loads_MSqRob}} to load R objects from RDatas binary format.
#' @export
inspect_loads_MSqRob <- function (file = NULL)
{
  tmp <- tempfile("saves.dir-")
  dir.create(tmp)
  untar(file, exdir = tmp)
  variables <- gsub(".RData","",dir(tmp))
  unlink(tmp, recursive = TRUE)
  return(variables)
}


