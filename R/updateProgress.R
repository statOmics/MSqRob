
#' @export
updateProgress <- function(progress=NULL, detail=NULL, n=NULL, shiny=FALSE, print=TRUE){
  if(isTRUE(shiny)){

    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = detail)

  } else if(isTRUE(print)){
    print(detail)
  }
}
