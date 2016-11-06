#source('directoryInput.R')
library(shiny)
library(DT)
library(shinyjs)
library(lme4)
library(MSqRob)
library(grDevices)
library(limma)
### Plotfunction for normalisation tab
plotnorm = function(p,method = 'quantiles',logp = TRUE,mdscol=1){
  if (logp) p = log(p,base = 2)
  ##Change -Inf values in the peptide intensities to NA
  exprs = Biobase::exprs(p)
  exprs[is.infinite(exprs)] = NA
  Biobase::exprs(p) = exprs

  unnorm = as.data.frame(p)
  p = normalise(p,method)
  norm = as.data.frame(p)
  
  dens = density(unlist(unnorm[1,]),na.rm = TRUE)
  densrest = lapply(2:nrow(unnorm),function(i){
    density(unlist(unnorm[i,]),bw = dens$bw,na.rm = TRUE)
  })
  densn = density(unlist(norm[1,]),na.rm = TRUE)
  densrestn = lapply(2:nrow(norm),function(i){
    density(unlist(norm[i,]),bw = densn$bw,na.rm = TRUE)
  })
  
  ###PLOT
  par(mfrow = c(2,2))
  ymax = max(dens$y,sapply(densrest,function(d){d$y}))
  xlim = range(unnorm,na.rm = TRUE)#,
  plot(dens,type = 'l',ylim = c(0,ymax),xlim = xlim, main = 'Normalisation method: None',col=mdscol[1])
  t = sapply(2:nrow(norm),function(i){
if (length(mdscol)>1) lines(densrest[[i-1]],col = mdscol[i]) else lines(densrest[[i-1]],col = i)
  })
  ymax = max(densn$y,sapply(densrestn,function(d){d$y}))
  xlim = range(norm,na.rm = TRUE)#,
  plot(densn,type = 'l',ylim = c(0,ymax),xlim = xlim, main = paste('Normalisation method:',method))
  t = sapply(2:nrow(norm),function(i){
if (length(mdscol)>1) lines(densrestn[[i-1]],col = mdscol[i]) else lines(densrestn[[i-1]],col = i)

  })
  plotMDS(t(unnorm),col=mdscol)  
  plotMDS(exprs(p),col=mdscol)
  
}


#Max file size: 500 MB
options(shiny.maxRequestSize=500*1024^2)


#source('directoryInput.R')
library(shiny)
library(DT)
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  #Function to check and process input

  processInput <- function(input, fileProtGrps, ann_name){

    if(isTRUE(input$onlysite) && is.null(fileProtGrps$datapath)){stop("Please provide a protein groups file or untick the box \"Remove proteins that are only identified by modified peptides\".")}

    if(input$save==2 && is.null(input$loadmodel$datapath)){stop("Please provide a saved RData file or don't choose the option \"Upload an existing model\" under \"Save/load options\".")}

    type_annot <- NULL

    if(isTRUE(as.logical(grep(".xlsx[/\\]*$",ann_name)))){type_annot <- "xlsx"}

    proteins <- gsub(" ",".",input$proteins)
    annotations <- gsub(" ",".",input$annotations)
    filter <- gsub(" ",".",input$filter)

    processedvals <- list("proteins"=proteins, "annotations"=annotations,"filter"=filter, "type_annot"=type_annot)
    return(processedvals)

  }

  #Needed for the input of a directory
  # observeEvent(
  #   ignoreNULL = TRUE,
  #   eventExpr = {
  #     input$directory
  #   },
  #   handlerExpr = {
  #     if (input$directory > 0) {
  #       # condition prevents handler execution on initial app launch
  #
  #       # launch the directory selection dialog with initial path read from the widget
  #       path = choose.dir(default = readDirectoryInput(session, 'directory'))
  #
  #       # update the widget value
  #       updateDirectoryInput(session, 'directory', value = path)
  #     }
  #   }
  # )

  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should
  #     re-execute automatically when inputs change
  #  2) Its output type is a plot

  filePep <- reactive({input$peptides})

  filterOptions <- reactive({
    if(is.null(filePep())){NULL
      } else{as.vector(as.matrix(read.table(filePep()$datapath, nrows=1, sep="\t", quote="")))}
  })

  selectedFilter <- reactive({
  if(!any(c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions())) {
    NULL
  }else{c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant")[c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions()]}
  })

  output$selectFilters <- renderUI({
    selectInput("filter", "Also filter based on these columns", filterOptions(), multiple=TRUE, selected=selectedFilter() )
  })


  fileAnn <- reactive({input$annotation})

  fixedOptions2 <- reactive({
    if(is.null(fileAnn()$name)){
      NULL
    } else{
    if(isTRUE(as.logical(grep(".xlsx[/\\]*$",fileAnn()$name)))){
      openxlsx::read.xlsx(fileAnn()$datapath)
    } else{
      read.table(fileAnn()$datapath, sep="\t", header=TRUE, row.names = NULL, quote="")
    }
  }
  })

  fixedOptions <- reactive({
    c(as.vector(colnames(fixedOptions2())),filterOptions())
  })

  levelOptions <- reactive({
    if(is.null(fixedOptions2()) || is.null(input$fixed)){
      NULL
    } else{
      optionsFixedSelected <- fixedOptions2()[,input$fixed,drop=FALSE]
      levelOptions <- unique(as.vector(as.matrix(sapply(colnames(optionsFixedSelected),function(name){ paste0(name,optionsFixedSelected[,name])}))))
      return(levelOptions)
      }
  })

 fileProtGrps <- reactive({input$proteingroups})

  #Without tabs:
  # output$selectLevels <- renderUI({
  #   #if(!exists("nContr")){nContr <- reactive({1})}
  #   if(!is.null(levelOptions())){
  #     lapply(1:input$nContr, function(j) {
  #     lapply(1:length(levelOptions()), function(i) {
  #         numericInput(paste0("contrast ",i,"_",j), levelOptions()[i], value=0, min = NA, max = NA, step = NA, width = NULL)
  #     })
  #     })
  #   }
  # })

  #With tabs:
  output$selectLevels = renderUI({
    nTabs = input$nContr
    myTabs = lapply(paste0('Contrast ', 1:nTabs), function(x){return(tabPanel(x,

                         if(!is.null(levelOptions())){
                         lapply(1:length(levelOptions()), function(i) {
                         numericInput(paste0("contrast ",i,"_",x), levelOptions()[i], value=0, min = NA, max = NA, step = NA, width = NULL)
                         })
                         }

                         ))})
    do.call(tabsetPanel, myTabs)
  })


  # fixedOptions <- reactive({
  #   if(is.null(fileAnn()$name)){
  #     NULL
  #   } else{
  #     if(isTRUE(as.logical(grep(".xlsx[/\\]*$",fileAnn()$name)))){
  #       c(as.vector(as.matrix(openxlsx::read.xlsx(fileAnn()$datapath, colNames=FALSE, rows=0:1))),filterOptions())
  #     } else{
  #       c(as.vector(as.matrix(read.table(fileAnn()$datapath, sep="\t", header=TRUE, row.names = NULL, nrows=1, quote=""))),filterOptions())
  #     }
  #   }
  # })

  nmsFixedOptions <- reactive({names(fixedOptions2())})


  output$selectFixed <- renderUI({
    selectInput("fixed", "Select fixed effects", nmsFixedOptions(), multiple=TRUE ) #"fixedOptions()" if you want to allow everything as a fixed effect, but then more difficult for contrasts
  })

  selectedRandom <- reactive({
    if(!c("Sequence") %in% filterOptions()) {
      NULL
    }else{"Sequence"}
  })

  output$selectRandom <- renderUI({
    selectInput("random", "Select random effects", fixedOptions(), multiple=TRUE, selected=selectedRandom() )
  })

  selectedProteins <- reactive({
    if(!c("Proteins") %in% filterOptions()) {
      NULL
    }else{"Proteins"}
  })

  output$selectProteins <- renderUI({
    selectInput("proteins", "Select the grouping factor (mostly the \"Proteins\" column)", filterOptions(), multiple=FALSE, selected=selectedProteins() )
  })

  selectedAnnotations <- reactive({
    if(!any(c("Gene names", "Protein names","Gene.names", "Protein.names") %in% filterOptions())) {
      NULL
    }else{c("Gene names", "Protein names","Gene.names", "Protein.names")[c("Gene names", "Protein names","Gene.names","Protein.names") %in% filterOptions()]}
  })

  output$selectAnnotations <- renderUI({
    selectInput("annotations", "Select additional annotation columns you want to keep", filterOptions(), multiple=TRUE, selected=selectedAnnotations() )
  })

  output$download_button <- renderUI({
    if(!is.null(outputlist())){
      downloadButton("downloadData", "Download")
    }
  })

  #Here comes what happens when we activate the go button, here are the real calculations
  #Maybe with progress bar...

  outputlist <- eventReactive(input$go, {

    outputlist=list(RData=list(proteins=NULL, models=NULL), test=NULL, results=NULL)

    if(input$save==2){load(input$load_model$datapath)
      outputlist$RData$proteins <- proteins
      outputlist$RData$models <- models
    }else{

      fs <- list()
      fs_type <- NULL

      processedvals = processInput(input, fileProtGrps(), fileAnn()$name)
      peptides = read_MaxQuant(filePep()$datapath, pattern="Intensity.")
      useful_properties = unique(c(processedvals[["proteins"]],processedvals[["annotations"]],input$fixed,input$random)[c(processedvals[["proteins"]],processedvals[["annotations"]],input$fixed,input$random) %in% colnames(Biobase::fData(peptides))])
      peptides2 = preprocess_MaxQuant(peptides, accession=processedvals[["proteins"]], annotation=fileAnn()$datapath, type_annot=processedvals[["type_annot"]], logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, useful_properties=useful_properties, filter=processedvals[["filter"]], remove_only_site=input$onlysite, file_proteinGroups=fileProtGrps()$datapath,  filter_symbol="+", minIdentified=input$minIdentified)

      #save(peptides2, file)


      #save(peptides2,processedvals, file="/Users/lgoeminn/MSqRob/test.RData")

      #levels(Biobase::pData(peptides2)$treat)
      #head(Biobase::exprs(peptides2))[2,]
      #levels(Biobase::fData(peptides2)$Proteins)

      #peptides2 <- peptides2[1:30]
      Biobase::fData(peptides2) <- droplevels(Biobase::fData(peptides2))
      proteins = MSnSet2protdata(peptides2, accession=processedvals[["proteins"]], annotations=processedvals[["annotations"]])

      #save(proteins, file="/Users/lgoeminn/MSqRob/proteins_test.RData")

      models <- fit.model(protdata=proteins, response="value", fixed=input$fixed, random=input$random)

      outputlist$RData$proteins <- proteins
      outputlist$RData$models <- models

      #If no save
      if(input$save==3){
        outputlist$RData$proteins <- NULL
        outputlist$RData$models <- NULL
      }
    }

    nTabs <- input$nContr
    L <- matrix(0, nrow=length(levelOptions()), ncol=nTabs)

    colnames(L) <- paste0('Contrast ', 1:nTabs)
    rownames(L) <- levelOptions()
    for(x in 1:nTabs){
                  if(!is.null(levelOptions())){
                  for(i in 1:length(levelOptions())){
                  L[i,x] <- input[[paste0("contrast ",i,"_Contrast ",x)]]
                  }
                  }
    }

    #If standard
    if(input$analysis_type=="standard"){
    RidgeSqM <- test.contrast_adjust(models, L, simplify=FALSE)

    #If stagewise
    } else if(input$analysis_type=="stagewise"){
    RidgeSqM <- test.contrast_stagewise(models, L, simplify=FALSE)

    #If ANOVA
    } else if(input$analysis_type=="ANOVA"){
    RidgeSqM <- test.contrast_adjust(models, L, simplify=FALSE, anova=TRUE, anova.na.ignore=FALSE)
    names(RidgeSqM) <- "ANOVA"
    }

    outputlist$results <- RidgeSqM
    outputlist$test <- "DONE!"

    return(outputlist)

  })

  estimate <- reactive({
    if(input$analysis_type %in% c("standard","stagewise")){
      estimate <- "estimate"
    } else if(input$analysis_type=="ANOVA"){
      estimate <- "AveExpr"
    }
    return(estimate)
  })


  output$downloadData <- downloadHandler(filename = function() { paste0(input$project_name,".zip") }, #default name
                                         content = function(file){

                                           #Works, but all folders on top are also included in the zip:

                                           #!!!Temporary folder will be removed, careful when changing this!!!
                                           temppath <- file.path(getwd(), "temp")
                                           dir.create(temppath)

                                           proteins <- outputlist()$RData$proteins
                                           models <- outputlist()$RData$models
                                           save(proteins,models, file=file.path(temppath, "models.RData"))

                                           results <- outputlist()$results
                                           
                                           if(Sys.info()['sysname']=="Windows"){
                                             names_res_xlsx <- character()
                                             for(i in 1:length(results)){
                                               names_res_xlsx[i] <- paste0("results",i,".xlsx")
                                           xlsx::write.xlsx(results[[i]], file = file.path(temppath, names_res_xlsx[i]), col.names = TRUE, row.names = TRUE)
                                             }
                                             files <- file.path(temppath, "models.RData")
                                             for(i in 1:length(results)){
                                               files[(i+1)] <- file.path(temppath, names_res_xlsx[i])
                                             }
                                             zip(zipfile=file, files=files)
                                           }else{
                                           openxlsx::write.xlsx(results, file = file.path(temppath, "results.xlsx"), colNames = TRUE, rowNames = TRUE)
                                            zip(zipfile=file, files=c(file.path(temppath, "models.RData"),file.path(temppath, "results.xlsx")))
                                             }
                                           
 
                                           #!!!Removes temporary folder, careful when changing this!!!!
                                           unlink(temppath, recursive = TRUE)

                                         },
                                         contentType = "application/zip"
  )

  ranges <- reactiveValues(x = NULL, y = NULL)

  #Een dropdown om het contrast te kiezen voor de plots

  contrastOptions <- reactive({
    nTabs = input$nContr
    paste0('Contrast ', 1:nTabs)
  })

  output$plot_contrast <- renderUI({
    if(input$analysis_type!="ANOVA"){
    selectInput("plot_contrast", "Select the contrast you want to visualize", contrastOptions())
    }
  })

  dataset <- reactive({
    if(input$analysis_type %in% c("standard","stagewise")){
      dataset <- outputlist()$results[[input$plot_contrast]]
    } else if(input$analysis_type=="ANOVA"){
      dataset <- outputlist()$results[["ANOVA"]]
    }
    #!!! "as.numeric:"  Quick fix voor ANOVA waarbij alles NA is (e.g. data Emmy, treatKO-treatWT en treatKO_LPS_1h-treatWT_LPS_1h), verder verfijnen!!!!:
    dataset$minus_log10_p <- -log10(as.numeric(dataset$pval)) #Moet er zijn omdat anders yvar niet gevonden kan worden!
    dataset <- data.frame(Accessions=rownames(dataset), dataset)
    rownames(dataset) <- NULL
    return(dataset)
  })


 
  output$plot1 <- renderPlot({

    #!!!Quick fix voor ANOVA waarbij alles NA is (e.g. data Emmy, treatKO-treatWT en treatKO_LPS_1h-treatWT_LPS_1h), verder verfijnen!!!!:
    if(!all(is.na(dataset()$minus_log10_p))){
    colBool <- dataset()$qval<0.05#subset(dataset, dataset$qval<0.05)
    colors <- rep(NA,length(dataset()$qval))
    colors[colBool] <- "red"
    colors[!colBool] <- "black"

    if(input$analysis_type %in% c("standard","stagewise")){
      xlab <- "estimate"
    } else if(input$analysis_type=="ANOVA"){
      xlab <- "average expression"
    }

    plot(dataset()[[estimate()]], dataset()$minus_log10_p, main="Volcano plot MSqRob", xlab=xlab, ylab="-log10(p)", xlim = ranges$x, ylim = ranges$y, las=1, col=colors, frame=FALSE)
    #points(sign_MSqRob$estimate, sign_MSqRob$minus_log10_p, col="red")

    s = input$table_rows_selected
    #Door het selecteren verandert de plot...

    if (length(s)) {
      subdataset <- clickInfo()[s, , drop = FALSE]

      colBool2 <- subdataset$qval<0.05
      colors2 <- rep(NA,length(subdataset$qval))
      colors2[colBool2] <- "red"
      colors2[!colBool2] <- "black"

      points(subdataset[[estimate()]], subdataset$minus_log10_p, pch = 19, cex = 2, col=colors2)
    }
    }
  })

   #When a double-click happens, check if there's a brush on the plot.
   #If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)

    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }

    #Set selection to zero: happens already if ranges change, but should also happen on normal double click
    proxy %>% DT::selectRows(NULL)

  })


  clickInfo <- reactive({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
    # if(nrow(brushedPoints(dataset(), input$plot1_brush, xvar="estimate", yvar="minus_log10_p"))==0){
      if(!is.null(ranges$x) && !is.null(ranges$y)){clickInfo <- subset(dataset(), (dataset()$estimate>ranges$x[1] & dataset()$estimate<ranges$x[2] & dataset()$minus_log10_p>ranges$y[1] & dataset()$minus_log10_p<ranges$y[2]))
      } else if(is.null(ranges$x) && is.null(ranges$y)){clickInfo <- dataset()}
    # } else{
    #   clickInfo <- brushedPoints(dataset(), input$plot1_brush, xvar="estimate", yvar="minus_log10_p")
    #   #Raar soort bug: enkel bij brushedPoints en als je DT::renderDataTable gebruikt, 2 lege rijen onderaan toegevoegd in de print...
    #   #Verwijder 2 laatste rijen als ze niet voorkomen in de dataset
    #   if(!(clickInfo[nrow(clickInfo),1] %in% dataset()[,1]) && !(clickInfo[(nrow(clickInfo)-1),1] %in% dataset()[,1])){clickInfo <- clickInfo[-c(nrow(clickInfo),(nrow(clickInfo)-1)),]}
    # }
    return(clickInfo)
  })

data <- reactive(
{
  data <- clickInfo()
  #oldnames <- c("se","df","Tval","pval","qval","signif","pvalS1","qvalS1","signifS1","AveExpr","df_num","df_den","Fval")
  #newnames <- c("standard error","degrees of freedom","T value","p value", "false discovery rate","significant","p value stage 1","false discovery rate stage 1","significant stage 1","average expression","degrees of freedom numerator","degrees of freedom denominator","F value")

  #for(i in 1:length(oldnames)){
  #  colnames(data)[colnames(data)==oldnames[i]] <- newnames[i]
  #}
  return(data)
  }
)



output$table<-DT::renderDataTable(data(),selection="single")


#output$table <- DT::renderDataTable( #DT version
 #   data()  ,selection='single'  ,  options = list(
 #                                      pageLength = 10,
 #                                      initComplete = I("function(settings, json) {alert('Done.');}")
 #                                    ))






  proxy = dataTableProxy('table')

  #clickInfo <- nearPoints(dataset, input$plot1_click, addDist = TRUE, xvar="estimate", yvar="minus_log10_p")

  #Add and remove points by clicking in the plot window
  observeEvent(input$plot1_click, {

    selected <- nearPoints(clickInfo(), input$plot1_click, addDist = TRUE,maxpoints=1, xvar=estimate(), yvar="minus_log10_p")
    sel_rows <- which(clickInfo()$Accessions %in% selected$Accessions)
    #Rows which were selected and selected again are removed, rows which were already selected but not selected again are retained
    #Don't sort this! Otherwise reacalculated.
    new_rows <- c(sel_rows[!sel_rows%in%input$table_rows_selected], input$table_rows_selected[!input$table_rows_selected%in%sel_rows])

    proxy %>% DT::selectRows(new_rows)
    })

#IKweg  #Enable or disable add brush to selection and remove brush from selection buttons
  observe({
    if (is.null(input$plot1_brush)) {
     shinyjs::disable("add_area_selection")
     shinyjs::disable("remove_area_selection")
    } else {
      shinyjs::enable("add_area_selection")
      shinyjs::enable("remove_area_selection")
    }
  })

  observeEvent(input$add_area_selection, {

    selected <- brushedPoints(clickInfo(), input$plot1_brush, xvar=estimate(), yvar="minus_log10_p")
    sel_rows <- which(clickInfo()$Accessions %in% selected$Accessions)
    #Rows which were selected and selected again are removed, rows which were already selected but not selected again are retained
    #Don't sort this! Otherwise reacalculated.
    new_rows <- unique(c(input$table_rows_selected,sel_rows))

    proxy %>% DT::selectRows(new_rows)
  })

  observeEvent(input$remove_area_selection, {

    selected <- brushedPoints(clickInfo(), input$plot1_brush, xvar=estimate(), yvar="minus_log10_p")
    sel_rows <- which(clickInfo()$Accessions %in% selected$Accessions)
    #Rows which were selected and selected again are removed, rows which were already selected but not selected again are retained
    #Don't sort this! Otherwise reacalculated.
    new_rows <- input$table_rows_selected[!(input$table_rows_selected %in% sel_rows)]

    proxy %>% DT::selectRows(new_rows)
  })

  #Drop down menu for plot 2
  plot2DependentVars <- reactive({
    as.list(c(input$fixed, input$random))
  })

  plot2OtherVars <- reactive({
    as.list(c("none",input$fixed, input$random))
  })

  output$selectPlot2 <- renderUI({
    selectInput("selPlot2", "Select independent variable", plot2DependentVars())
    })

  output$selectColPlot2 <- renderUI({
    selectInput("selColPlot2", "Select color variable", plot2OtherVars())
  })

  #Plot 2

  acc_plot2 <- reactive({as.character(clickInfo()[input$table_rows_selected,"Accessions"])}) #[input$table_rows_selected,"Accessions"]

  indep_var_plot2 <- reactive({input$selPlot2})
  color_var_plot2 <- reactive({input$selColPlot2})


  output$plot2 <- renderPlot({

    accessions <- acc_plot2() #Geeft "NA"
    proteins <- outputlist()$RData$proteins
    indep_var <- indep_var_plot2()
    color_var <- color_var_plot2()
    if(color_var=="none"){colors <- 1} else{
    colordata <- tryCatch(getData(proteins[accessions])[[color_var]], error=function(e){
      return(NULL)
    })

    #http://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(length(unique(colordata)))
    colors <- colors[as.numeric(as.factor(colordata))]
    }

    if(length(accessions)==1){

      #proteins[acc_plot2()]$value
      boxplot(getData(proteins[accessions])$value~getData(proteins[accessions])[[indep_var]], outline=FALSE, ylim=c(min(getData(proteins[accessions])$value)-0.2,max(getData(proteins[accessions])$value)+0.2), ylab="log2(peptide intensity)", xlab="", main=getAccessions(proteins[accessions]), las=2, frame.plot=FALSE, frame=FALSE, col="grey", pars=list(boxcol="white")) #, cex.main=2, cex.lab=2, cex.axis=2, cex=2
      points(jitter((as.numeric(getData(proteins[accessions])[[indep_var]])), factor=2),getData(proteins[accessions])$value, col=colors) #,cex=2, lwd=2, col=c(1,2,3,4,"cyan2",6)
      # title(ylab="Log2(Intensity)", line=5, cex.lab=2, family="Calibri Light")

    } else{NULL}

  })


  # Nieuwe select button maken voor brushed area
  # observeEvent(input$select1, {
  #   proxy %>% selectRows(as.numeric(input$rows))
  # })

  output$nText <- renderText({
    outputlist()$test
  })

 #######################################################################
  ## Normalisation tab 
  
#   v <- reactiveValues(peps = NULL,normmethod = 'quantiles',normlog = TRUE)
#  observeEvent(input$norm1, {
#    v$normmethod = input$norm1
#  })
#
#  observeEvent(input$peptides, {
#    v$peps = read_MaxQuant(input$peptides$datapath)
#    })

#  output$plotnorm1 = renderPlot({
#     if(!is.null(v$peps)){
#       plotnorm(v$peps,v$normmethod,v$normlog)
#   }},
#   width = 700, height = 800)

#})

  #Drop down menu for plot normalization Plot
  plotNorm1DependentVars <- reactive({
    as.list(c("none",colnames(fixedOptions2())))
  })

   output$selectColPlotNorm1 <- renderUI({
    selectInput("selColPlotNorm1", "Select variable for visualisation",  plotNorm1DependentVars())
  })


   v <- reactiveValues(peps = NULL,normmethod = 'quantiles',normlog = TRUE)
  observeEvent(input$norm1, {
    v$normmethod = input$norm1
  })

  observeEvent(input$peptides, {
    v$peps = read_MaxQuant(input$peptides$datapath)
    })
  
  color_var_plotnorm1 <- reactive({input$selColPlotNorm1})

  output$plotnorm1<- renderPlot(
{
    color_var <- color_var_plotnorm1()
    colors <- 1
    try(
    {colordata <- fixedOptions2()[,color_var]
    #http://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
    #colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
    #                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(length(unique(colordata)))
    colors<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(length(unique(colordata)))
    colors <- colors[as.numeric(as.factor(colordata))]
    },silent=TRUE) 
    if(!is.null(v$peps)){
       plotnorm(v$peps,v$normmethod,v$normlog,mdscol=colors)}}
   ,
   width = 700, height = 800)
#plot(1,1))
})
