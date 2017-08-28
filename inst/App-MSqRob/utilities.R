
#A function to prompt a folder chooser under Mac OS X,
#normal choose.dir sometimes gives just NA because system has no access to finder
choose.dir2 <- function(default = NA, caption = NA) {
  command = 'osascript'

  #Dit was het:
  #args = "-e 'tell application \"Finder\"' -e 'activate' -e 'POSIX path of (choose folder{{prompt}}{{default}})' -e 'end tell'"
  #'-e "POSIX path of (choose folder{{prompt}}{{default}})"'

  #osascript -e 'tell application "Terminal" to return POSIX path of (choose file)'

  #Find App that is in front
  args1 = "-e 'tell application \"System Events\"' -e 'set frontApp to name of first application process whose frontmost is true' -e 'end tell'"
  suppressWarnings({
    frontApp = system2(command, args = args1, stderr = TRUE)
  })

  #The application that is in front should open the  file choose window
  args = paste0("-e 'tell application \"",frontApp,"\" to return POSIX path of (choose folder{{prompt}}{{default}})'")

  #system2('osascript', args = "-e 'tell application \"Finder\"' -e 'POSIX path of (choose folder{{prompt}}{{default}})' -e 'end tell'", stderr = TRUE)

  if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
    prompt = sprintf(' with prompt \\"%s\\"', caption)
  } else {
    prompt = ''
  }
  args = sub('{{prompt}}', prompt, args, fixed = TRUE)

  if (!is.null(default) && !is.na(default) && nzchar(default)) {
    default = sprintf(' default location \\"%s\\"', path.expand(default))
  } else {
    default = ''
  }
  args = sub('{{default}}', default, args, fixed = TRUE)

  suppressWarnings({
    path = system2(command, args = args, stderr = TRUE)
  })
  if (!is.null(attr(path, 'status')) && attr(path, 'status')) {
    # user canceled
    path = NA
  }

  return(path)
}

choose.dir_Linux <- function(default = NA, caption = NA) {

  command = "zenity"

  args1 = "--file-selection --directory"

  if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
    prompt = sprintf(' with prompt \\"%s\\"', caption)
  } else {
    prompt = ''
  }

  if (!is.null(default) && !is.na(default) && nzchar(default)) {
    default = path.expand(default) #sprintf(' default location \\"%s\\"', path.expand(default))
  } else {
    default = "/dev/null"
  }

  suppressWarnings({
    path = system2(command, args = args1, stderr = default, stdout=TRUE)
  })

  if (!is.null(attr(path, 'status')) && attr(path, 'status')) {
    # user canceled
    path = NA
  }

  return(path)
}


#Function to convert data paths if on windows
getDataPath <- function(datapath){
  if(Sys.info()['sysname']=="Windows"){
    datapath <- gsub("\\","/",datapath, fixed=TRUE)
  }
  return(datapath)
}


#Function to plot densities
plotDens=function(eset, densAll, xlim=NULL, ylim=NULL, colors=1, las=1, frame.plot=FALSE, ...)
{
      plot(densAll[[1]],col=colors[1],xlim=xlim,ylim=ylim, las=las, frame.plot=frame.plot, ...)
      if (length(colors)>1) for (i in 2:ncol(eset)) lines(densAll[[i]],col=colors[i])
      else for (i in 2:ncol(eset)) lines(densAll[[i]],col=colors)
}

#Function to check and process input
processInput <- function(input){

    if(isTRUE(input$onlysite) && is.null(input$proteingroups$datapath)){stop("Please provide a protein groups file or untick the box \"Remove proteins that are only identified by modified peptides\".")}

    if(input$save==2 && is.null(input$loadmodel$datapath)){stop("Please provide a saved RData file or don't choose the option \"Upload an existing model\" under \"Save/load options\".")}

    type_annot <- NULL

    if(isTRUE(as.logical(grep(".xlsx[/\\]*$",input$annotation$name)))){type_annot <- "xlsx"}

    proteins <- input$proteins
    annotations <- input$annotations
    filter <- input$filter

    if(!is.null(proteins)){proteins <- gsub(" ",".",proteins)}
    if(!is.null(annotations)){annotations <- gsub(" ",".",annotations)}
    if(!is.null(filter)){filter <- gsub(" ",".",filter)}

    peptides <- input$peptides

    peptides$datapath <- getDataPath(as.character(peptides$datapath))

    processedvals <- list("proteins"=proteins, "peptides"=peptides, "annotations"=annotations,"filter"=filter, "type_annot"=type_annot)
    return(processedvals)

  }


#' Folder Upload Control
#'
#' Create a folder upload control that can be used to upload one or more filepaths pointing to folders. Strongly based on Shiny's File Upload Control.
#'
#' Whenever a folder upload completes, the corresponding input variable is set
#' to a character path.
#'
#' @family input elements
#'
#' @param inputId	The \code{input} slot that will be used to access the value.
#' @param label	Display label for the control, or \code{NULL} for no label.
#' @param value	Initial value.
#' @param width	The width of the input, e.g. \code{'400px'}, or \code{'100%'}; see \code{\link{validateCssUnit}}.
#' @param placeholder	A character string giving the user a hint as to what can be entered into the control. Internet Explorer 8 and 9 do not support this option.
#' @param multiple Whether the user should be allowed to select and upload
#'   multiple folders at once. \bold{Does not work on older browsers, including
#'   Internet Explorer 9 and earlier.}
#' @param accept A character vector of MIME types; gives the browser a hint of
#'   what kind of folders the server is expecting.
#' @param style The style attribute for the HTML tag. Used to hide/unhide the progress bar.
#'
#' @examples
#' ## Only run examples in interactive R sessions
#' if (interactive()) {
#'
#' ui <- fluidPage(
#'   sidebarLayout(
#'     sidebarPanel(
#'       fileInput("file1", "Choose CSV File",
#'         accept = c(
#'           "text/csv",
#'           "text/comma-separated-values,text/plain",
#'           ".csv")
#'         ),
#'       tags$hr(),
#'       checkboxInput("header", "Header", TRUE)
#'     ),
#'     mainPanel(
#'       tableOutput("contents")
#'     )
#'   )
#' )
#'
#' server <- function(input, output) {
#'   output$contents <- renderTable({
#'     # input$file1 will be NULL initially. After the user selects
#'     # and uploads a file, it will be a data frame with 'name',
#'     # 'size', 'type', and 'datapath' columns. The 'datapath'
#'     # column will contain the local filenames where the data can
#'     # be found.
#'     inFile <- input$file1
#'
#'     if (is.null(inFile))
#'       return(NULL)
#'
#'     read.csv(inFile$datapath, header = input$header)
#'   })
#' }
#'
#' shinyApp(ui, server)
#' }
#' @export
folderInput <- function(inputId, label, value = NA, multiple = FALSE, accept = NULL,
                        width = NULL, style="") {

  restoredValue <- restoreInput(id = inputId, default = NULL)

  # Catch potential edge case - ensure that it's either NULL or a data frame.
  if (!is.null(restoredValue) && !dir.exists(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }

  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }

  inputTag <- tags$input(
    id = inputId,
    name = inputId,
    type = "button",
    style = "display: none;",
    `data-restore` = restoredValue,
    class = "btn action-button"
    # webkitdirectory = NA,
    # directory = NA
  )

  if (multiple)
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse=',')


  div(class = "form-group shiny-input-container",
      style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"),
      shiny:::`%AND%`(label, tags$label(label)),

      div(class = "input-group",
          tags$label(class = "input-group-btn",
                     span(id = paste(inputId, "_label", sep = ""), class = "btn btn-default btn-file",
                          "Browse...",
                          inputTag
                     )
          ),
          tags$input(type = "text", class = "form-control", value=value,
                     placeholder = "No folder selected", readonly = "readonly"
          )
      ),

      tags$div(
        id=paste(inputId, "_progress", sep=""),
        class="progress progress-striped active shiny-file-input-progress", style=style, #"visibility: visible;"
        tags$div(class="progress-bar", style="width: 100%;","Folder selected")
      )
  )
}





#'@export
fileInput <- function (inputId, label, multiple = FALSE, accept = NULL, width = NULL)
{
  restoredValue <- restoreInput(id = inputId, default = NULL)
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  inputTag <- tags$input(id = inputId, name = inputId, type = "file",
                         style = "display: none;", `data-restore` = restoredValue)
  if (multiple)
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse = ",")
  div(class = "form-group shiny-input-container", style = if (!is.null(width))
    paste0("width: ", validateCssUnit(width), ";"),
    shiny:::`%AND%`(label, tags$label(label)),
    div(class = "input-group", tags$label(class = "input-group-btn",
                                                               span(id = paste(inputId, "_label", sep = ""), class = "btn btn-default btn-file", "Browse...",
                                                                    inputTag)), tags$input(type = "text", class = "form-control",
                                                                                           placeholder = "No file selected", readonly = "readonly")),
    tags$div(id = paste(inputId, "_progress", sep = ""),
             class = "progress progress-striped active shiny-file-input-progress",
             tags$div(class = "progress-bar")))
}
