#source('directoryInput.R')
library(shiny)
library(DT)
library(shinyjs)
library(lme4)
library(MSqRob)
library(MSnbase)
library(grDevices)
library(limma)
source("utilities.R")

#Max file size: 500 MB
options(shiny.maxRequestSize=500*1024^2)
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  ###########################################
  #Input Tab
  ###########################################

  saveFolder <- reactiveValues(folder = NULL)

  #Function to check if save folder is set
  check_save_folder <- function(save_folder){
    if(is.null(save_folder) | length(save_folder)==0){stop("No output folder selected!")
    } else if(is.na(save_folder)){
      stop("No output folder selected!")
    } else if(!dir.exists(save_folder)){
      stop("Output folder does not exist!")
    }
  }

  output$folderError <- renderText({
    check_save_folder(saveFolder$folder)
  })

  observeEvent(input$outputFolder,{

    if(Sys.info()['sysname']=="Darwin"){
      saveFolder$folder <- tryCatch(
        choose.dir2()
        , error=function(e){
          library(svDialogs)
          svDialogs::dlgDir()$res
        })
    } else{
      saveFolder$folder <- tryCatch(
        choose.dir()
        , error=function(e){
          library(svDialogs)
          svDialogs::dlgDir()$res
        })
    }
  })

  output$outputFolderOut <- renderUI({

    if(is.null(saveFolder$folder) | length(saveFolder$folder)==0){
      style=""
      value = NA
    } else if(is.na(saveFolder$folder)){
      style=""
      value = NA
    } else if(!dir.exists(saveFolder$folder)){
      style=""
      value = NA
    } else{
      style="visibility: visible;"
      value = saveFolder$folder
    }

    folderInput(inputId="outputFolder", label="Specify the location where your output will be saved", value = value, multiple = FALSE, accept = NULL, width = NULL, style=style)
  })

  observe({

    #If no save folder or no peptides.txt file specified, do not allow to try to create an annotation data frame.
    res <- try(check_save_folder(saveFolder$folder),silent = TRUE)

    if(is.null(input$peptides) | class(res) == "try-error"){
      shinyjs::disable("create_annot")
    } else{
      shinyjs::enable("create_annot")
    }
  })

  #Button to initialize experimental annotation file
  newExpAnnText <- eventReactive(input$create_annot, {
    #Check if save folder is set
    check_save_folder(saveFolder$folder)
    init_ann_MQ_Excel(peptidesDatapath(), savepath=saveFolder$folder, output_name=paste0(input$project_name,"_experimental_annotation"), col_name="run", pattern="Intensity.", remove_pattern=TRUE)
    newExpAnnText <- paste0("Annotation file initialized. Check ",saveFolder$folder,"/",input$project_name,"_experimental_annotation.xlsx. \n Adjust this file according to your experimental settings and upload it as your experimental annotation file.")
    return(newExpAnnText)
  })

  output$newExpAnnText <- renderText({
    newExpAnnText()
  })

  ########################################################
  #Clear datapaths of backslashes (Needed on Windows only)
  ########################################################

  peptidesDatapath <- reactive({getDataPath(input$peptides$datapath)})
  annotationDatapath <- reactive({getDataPath(input$annotation$datapath)})
  modelDatapath <- reactive({getDataPath(input$load_model$datapath)})
  proteinGroupsDatapath <- reactive({getDataPath(input$proteingroups$datapath)})

  ########################################
  #Set Filter option
  ########################################
  filterOptions <- reactive({
    if(is.null(input$peptides)){
      NULL
    } else{
      # req(input$peptides)
      as.vector(as.matrix(read.table(peptidesDatapath(), nrows=1, sep="\t", quote="")))
    }
  })
  selectedFilter <- reactive({
    if(!any(c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions())) {
      NULL
    } else{c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant")[c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions()]}
  })

  output$selectFilters <- renderUI({
    selectInput("filter", "Also filter based on these columns", filterOptions(), multiple=TRUE, selected=selectedFilter())})


  ########################################
  #Generate options for fixed effect variables
  ########################################
  #from annotation file
  fixedOptions2 <- reactive({
    if(is.null(input$annotation$name)){
      NULL
    } else{
      if(isTRUE(as.logical(grep(".xlsx[/\\]*$",input$annotation$name)))){

        #Dirty fix for Windows 7, 64 bits problem with rJava
        # if(Sys.info()['sysname']=="Windows"){
        #   test <- xlsx::read.xlsx(file=annotationDatapath(), sheetIndex = 1)
        #   #Remove columns with all NA
        #   test <- test[,!colSums(apply(test, 2, is.na))==nrow(test),drop=FALSE]
        #   return(test)
        # } else{
        openxlsx::read.xlsx(annotationDatapath())
        # }

      } else{
        read.table(annotationDatapath(), sep="\t", header=TRUE, row.names = NULL, quote="")
      }
    }
  })
  #from peptides file (filterOptions)
  fixedOptions <- reactive({
    c(as.vector(colnames(fixedOptions2())),filterOptions())
  })
  #Generate option of factor levels
  levelOptions <- reactive({
    if((is.null(fixedOptions2()) | is.null(input$fixed)) & (input$save!=2 | is.null(input$load_model$datapath))){
      NULL
    } else if(input$save==2 & !is.null(input$load_model$datapath)){
      #Load models
      # progressSave <- NULL
      # # Create a Progress object
      # progressSave <- shiny::Progress$new()
      #
      # # Make sure it closes when we exit this reactive, even if there's an error
      # on.exit(progressSave$close())
      # progressSave$set(message = "Loading settings...", value = 0)

      #Load only levelOptions (much faster!)
      RData <- try(loads_MSqRob(file=modelDatapath(), variables="levelOptions"), silent=TRUE)
      if(inherits(RData, 'try-error')){stop("Loading of model file failed. Please provide a valid RDatas model file.")}
      levelOptions <- RData$levelOptions
    } else{
      optionsFixedSelected <- fixedOptions2()[,input$fixed,drop=FALSE]
      levelOptions <- unique(as.vector(as.matrix(sapply(colnames(optionsFixedSelected),function(name){ paste0(name,optionsFixedSelected[,name])}))))
      return(levelOptions)
    }
  })

  ###########################################
  #Functionalities for Quantification Tab
  ###########################################
  nmsFixedOptions <- reactive({names(fixedOptions2())})

  ####select Fixed effects, random effects, Proteins and store options ####
  output$selectFixed <- renderUI({
    selectInput("fixed", "Select fixed effects", nmsFixedOptions(), multiple=TRUE )
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




  ####set contrasts###
  #With tabs:
  output$selectLevels = renderUI({
    nTabs = input$nContr
    myTabs = lapply(paste0('Contrast ', 1:nTabs), function(x){return(tabPanel(x,

                                                                              if(!is.null(levelOptions())){
                                                                                lapply(1:length(levelOptions()), function(i) {

                                                                                  #Preserve values of input :)
                                                                                  isolate(
                                                                                    if(!is.null(input[[paste0("contrast ",i,"_",x)]])){value <- input[[paste0("contrast ",i,"_",x)]]
                                                                                    } else{value <- 0}
                                                                                  )

                                                                                  textInput(paste0("contrast ",i,"_",x), levelOptions()[i], value=value, width = NULL) #, min = NA, max = NA, step = NA,
                                                                                })
                                                                              }

    ))})
    do.call(tabsetPanel, myTabs)
  })

  ###Select analysis type
  estimate <- reactive({
    if(input$analysis_type %in% c("standard","stagewise")){
      estimate <- "estimate"
    } else if(input$analysis_type=="ANOVA"){
      estimate <- "AveExpr"
    }
    return(estimate)
  })


  ###Generation of output button###
  # output$download_button <- renderUI({
  #   if(!is.null(outputlist())){
  #     downloadButton("downloadData", "Download")
  #   }
  # })
  ###Function invoked when output button is pushed###
  #Here comes what happens when we activate the go button, here are the real calculations
  #Maybe with progress bar...

  observe({

    #Observe input$filter so that when going immediately to tab 3, select load model and then go to tab 2 leads also to disabled input$filter field.
    input$filter

    # Change on selection of tab?
    # input$Input
    # input$Quantification
    # input$Preprocessing

    if(input$save!=2){
      shinyjs::enable("peptides")
      shinyjs::enable("peptides_label")
      shinyjs::enable("annotation")
      shinyjs::enable("annotation_label")
      shinyjs::enable("onlysite")
      shinyjs::enable("proteingroups")
      shinyjs::enable("proteingroups_label")
      shinyjs::enable("smallestUniqueGroups")
      shinyjs::enable("minIdentified")
      shinyjs::enable("filter")
      shinyjs::enable("logtransform")
      shinyjs::enable("log_base")
      shinyjs::enable("normalisation")
      shinyjs::enable("proteins")
      shinyjs::enable("annotations")
      shinyjs::enable("fixed")
      shinyjs::enable("random")

      shinyjs::enable("evalnorm")
      shinyjs::enable("selColPlotNorm1")
      shinyjs::enable("plotMDSPoints")
      shinyjs::enable("plotMDSLabels")

      shinyjs::enable("borrowFixed")
      shinyjs::enable("borrowRandom")

    } else if(input$save==2){
      shinyjs::disable("peptides")
      shinyjs::disable("peptides_label")
      shinyjs::disable("annotation")
      shinyjs::disable("annotation_label")
      shinyjs::disable("onlysite")
      shinyjs::disable("proteingroups")
      shinyjs::disable("proteingroups_label")
      shinyjs::disable("smallestUniqueGroups")
      shinyjs::disable("minIdentified")
      shinyjs::disable("filter")
      shinyjs::disable("logtransform")
      shinyjs::disable("log_base")
      shinyjs::disable("normalisation")
      shinyjs::disable("proteins")
      shinyjs::disable("annotations")
      shinyjs::disable("fixed")
      shinyjs::disable("random")

      shinyjs::disable("evalnorm")
      shinyjs::disable("selColPlotNorm1")
      shinyjs::disable("plotMDSPoints")
      shinyjs::disable("plotMDSLabels")

      shinyjs::disable("borrowFixed")
      shinyjs::disable("borrowRandom")
    }
  })

  outputlist <- eventReactive(input$go, {

    validate(
      need((!is.null(saveFolder$folder) & length(saveFolder$folder)!=0), "No output folder selected!")
    )

    validate(
      need((input$save==2 | !is.null(input$fixed)), "Please select at least one fixed effect!")
    )

    nTabs <- input$nContr
    L <- matrix(0, nrow=length(levelOptions()), ncol=nTabs)

    colnames(L) <- paste0('Contrast ', 1:nTabs)
    rownames(L) <- levelOptions()
    for(x in 1:nTabs){
      if(!is.null(levelOptions())){
        for(i in 1:length(levelOptions())){
          validate(
            #1. only numbers and mathematical operators
            need((grep("^[-0-9\\*\\+\\/\\.\\(\\)]*$", input[[paste0("contrast ",i,"_Contrast ",x)]])==1), "All contrast input should be numeric!"),
            #2. the mathematical expression needs to be valid!
            need(is.numeric(try(eval(parse(text=input[[paste0("contrast ",i,"_Contrast ",x)]])))), "All contrast input should be numeric!")
          )
          L[i,x] <- eval(parse(text=input[[paste0("contrast ",i,"_Contrast ",x)]]))
        }
      }
    }

    #Check if saveFolder is correctly specified!
    check_save_folder(saveFolder$folder)

    outputlist=list(RData=list(proteins=NULL,
                               models=NULL),
                    test=NULL,
                    results=NULL)

    if(input$save==2){

      #Load models: levelOptions and plot2DependentVars are loaded in their respective reactives!
      RData <- try(loads_MSqRob(file=modelDatapath(), variables=c("proteins","random","models"), printProgress=TRUE, shiny=TRUE, message="Loading models..."), silent=TRUE)
      if(class(RData)=="try-error"){stop("Loading of model file failed. Please provide a valid RDatas model file.")}
      outputlist$RData$proteins <- RData$proteins
      outputlist$RData$models <- RData$models
      fixed <- RData$fixed
      random <- RData$random
      #levelOptions is loaded in "levelOptions" reactive, needs to work before "go" button is pressed!
    }else{

      fixed <- input$fixed
      random <- input$random

      fs <- list()
      fs_type <- NULL

      processedvals = processInput(input)
      peptides = read_MaxQuant(peptidesDatapath(), pattern="Intensity.", shiny=TRUE, message="Importing data...")

      useful_properties = unique(c(processedvals[["proteins"]],processedvals[["annotations"]],fixed,random)[c(processedvals[["proteins"]],processedvals[["annotations"]],fixed,random) %in% colnames(Biobase::fData(peptides))])
      peptides2 = preprocess_MaxQuant(peptides, accession=processedvals[["proteins"]], exp_annotation=annotationDatapath(), type_annot=processedvals[["type_annot"]], logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, useful_properties=useful_properties, filter=processedvals[["filter"]], remove_only_site=input$onlysite, file_proteinGroups=proteinGroupsDatapath(),  filter_symbol="+", minIdentified=input$minIdentified, shiny=TRUE, printProgress=TRUE, message="Preprocessing data...")
      Biobase::fData(peptides2) <- droplevels(Biobase::fData(peptides2))
      proteins = MSnSet2protdata(peptides2, accession=processedvals[["proteins"]], annotations=processedvals[["annotations"]], printProgress=TRUE, shiny=TRUE, message="Converting data...")

      par_squeeze <- NULL

      if(isTRUE(input$borrowRandom)){par_squeeze <- c(par_squeeze, random)}
      if(isTRUE(input$borrowFixed)){par_squeeze <- c(par_squeeze,"ridgeGroup.1")}

      models <- fit.model(protdata=proteins, response="quant_value", fixed=fixed, random=random, par_squeeze=par_squeeze, printProgress=TRUE, shiny=TRUE, message_fitting="Fitting models...", message_thetas="Extracting variances...", message_squeeze="Squeezing variances...", message_update="Updating models...")

      #We save the squeezed models!

      outputlist$RData$proteins <- proteins
      outputlist$RData$models <- models

    }

    outputlist$L <- L

    #If standard
    if(input$analysis_type=="standard"){
      RidgeSqM <- test.contrast_adjust(outputlist$RData$models, L, simplify=FALSE, printProgress=TRUE, shiny=TRUE, message_extract="Calculating contrasts...", message_test="Testing contrasts...")

      #If stagewise
    } else if(input$analysis_type=="stagewise"){
      RidgeSqM <- test.contrast_stagewise(outputlist$RData$models, L, simplify=FALSE, printProgress=TRUE, shiny=TRUE, message_extractS1="Calculating contrasts stage 1...", message_testS1="Testing contrasts stage 1...", message_extractS2="Calculating contrasts stage 2...", message_testS2="Testing contrasts stage 1...")

      #If ANOVA
    } else if(input$analysis_type=="ANOVA"){
      RidgeSqM <- test.contrast_adjust(outputlist$RData$models, L, simplify=FALSE, anova=TRUE, anova.na.ignore=FALSE, printProgress=TRUE, shiny=TRUE, message_extract="Calculating contrasts...", message_test="Testing contrasts...")
      names(RidgeSqM) <- "ANOVA"
    }

    outputlist$results <- RidgeSqM
    outputlist$test <-  "DONE!"

    ###Save output (unless no save: input$save==3)###
    if(input$save!=3){

      proteins <- outputlist$RData$proteins
      models <- outputlist$RData$models
      results <- outputlist$results

      savepath <- getDataPath(saveFolder$folder)
      savepath <- gsub("//","/",file.path(savepath, paste0(input$project_name,"_",gsub(" |:","_",Sys.time()))))
      dir.create(savepath)
      RData <- outputlist$RData #2 slots: "proteins" and "models"
      RData$levelOptions <- levelOptions()
      RData$fixed <- fixed
      RData$random <- random
      #RData$plot2DependentVars <- plot2DependentVars()

      RData_env <- new.env()
      assign("RData", RData, RData_env)
      wd_old <- getwd()
      setwd(savepath)
      saves_MSqRob(RData, envir=RData_env, file=paste0(input$project_name,"_","models.RDatas"), shiny=TRUE, printProgress=TRUE, message="Saving models")
      setwd(wd_old)

      progressSaveExcel <- NULL
      # Create a Progress object
      progressSaveExcel <- shiny::Progress$new()

      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progressSaveExcel$close())
      progressSaveExcel$set(message = "Saving results...", value = 0)

      #save(RData, file=file.path(savepath, paste0(input$project_name,"_","models.RDatas")))
      openxlsx::write.xlsx(results, file = file.path(savepath, paste0(input$project_name,"_","results.xlsx")), colNames = TRUE, rowNames = TRUE)

      updateProgress(progress=progressSaveExcel, detail="Saving to Excel file", n=1, shiny=TRUE, print=TRUE)

      # #Bold header in Excel file:
      # headerStyle <- openxlsx::createStyle(textDecoration = "bold")
      # openxlsx::addStyle(wb, sheet = 1:length(results), headerStyle, rows = 1, cols = 1:ncol(results[[1]])) #, gridExpand = TRUE
    }
    ######

    return(outputlist)

  })

  # output$downloadData <- downloadHandler(filename = function() { paste0(input$project_name,"_",gsub(" |:","_",Sys.time()),".zip") }, #default name
  #                                        content = function(file){
  #
  #                                          file <- getDataPath(file)
  #
  #                                          #Works, but all folders on top are also included in the zip:
  #
  #                                          #!!!Temporary folder will be removed, careful when changing this!!!
  #                                          temppath <- file.path(getwd(), "temp")
  #                                          dir.create(temppath)
  #
  #                                          proteins <- outputlist()$RData$proteins
  #                                          models <- outputlist()$RData$models
  #                                          save(proteins,models, file=file.path(temppath, "models.RData"))
  #
  #                                          results <- outputlist()$results
  #
  #                                          # if(Sys.info()['sysname']=="Windows"){
  #                                          #   names_res_xlsx <- character()
  #                                          #   for(i in 1:length(results)){
  #                                          #     names_res_xlsx[i] <- paste0("results",i,".xlsx")
  #                                          # xlsx::write.xlsx(results[[i]], file = file.path(temppath, names_res_xlsx[i]), col.names = TRUE, row.names = TRUE)
  #                                          #   }
  #                                          #   files <- file.path(temppath, "models.RData")
  #                                          #   for(i in 1:length(results)){
  #                                          #     files[(i+1)] <- file.path(temppath, names_res_xlsx[i])
  #                                          #   }
  #                                          #   zip(zipfile=file, files=files)
  #                                          # }else{
  #                                          openxlsx::write.xlsx(results, file = file.path(temppath, "results.xlsx"), colNames = TRUE, rowNames = TRUE)
  #                                           zip(zipfile=file, files=c(file.path(temppath, "models.RData"),file.path(temppath, "results.xlsx")))
  #                                            # }
  #
  #
  #                                          #!!!Removes temporary folder, careful when changing this!!!!
  #                                          unlink(temppath, recursive = TRUE)
  #
  #                                        },
  #                                        contentType = "application/zip"
  # )

  output$contrastL <- renderPrint({
    #The contrast corresponding to the selected contrast
    L <- outputlist()$L[,input$plot_contrast,drop=FALSE]
    L <- L[L!=0,,drop=FALSE]
    return(L)
  })


  ###########################################
  #Plots for Quantification Tab
  ###########################################
  ###ranges for zooming in our out in plot###
  ranges <- reactiveValues(x = NULL, y = NULL)

  ###Choice of contrast to be visualized ###
  contrastOptions <- reactive({
    nTabs = input$nContr
    paste0('Contrast ', 1:nTabs)
  })
  output$plot_contrast <- renderUI({
    if(input$analysis_type!="ANOVA"){
      selectInput("plot_contrast", "Select the contrast you want to visualize", contrastOptions())
    }
  })

  ###Generation of all data for output###
  dataset <- reactive({
    if(!is.null(outputlist())){
      if(input$analysis_type %in% c("standard","stagewise")){
        dataset <- outputlist()$results[[input$plot_contrast]]
      } else if(input$analysis_type=="ANOVA"){
        dataset <- outputlist()$results[["ANOVA"]]
      }
      #!!! "as.numeric:"  Quick fix voor ANOVA waarbij alles NA is (e.g. data Emmy, treatKO-treatWT en treatKO_LPS_1h-treatWT_LPS_1h), verder verfijnen!!!!:
      dataset$minus_log10_p <- -log10(as.numeric(dataset$pval)) #Necessary to select data in table according to the zoom in the plot
      dataset <- data.frame(Accessions=rownames(dataset), dataset)
      rownames(dataset) <- NULL
    } else{dataset <- NULL}
    return(dataset)
  })



  ###Volcanoplot###
  output$plot1 <- renderPlot({
    #!!!Quick fix voor ANOVA waarbij alles NA is (e.g. data Emmy, treatKO-treatWT en treatKO_LPS_1h-treatWT_LPS_1h), verder verfijnen!!!!:
    if(!all(is.na(dataset()$minus_log10_p))){
      colBool <- dataset()$qval<input$alpha
      colors <- rep(NA,length(dataset()$qval))
      colors[colBool] <- "red"
      colors[!colBool] <- "black"

      if(input$analysis_type %in% c("standard","stagewise")){
        xlab <- "estimate"
      } else if(input$analysis_type=="ANOVA"){
        xlab <- "average expression"
      }

      plot(dataset()[[estimate()]], dataset()$minus_log10_p, main="Volcano plot MSqRob", xlab=xlab, ylab="-log10(p)", xlim = ranges$x, ylim = ranges$y, las=1, col=colors, bty="n")
      #points(sign_MSqRob$estimate, sign_MSqRob$minus_log10_p, col="red")

      s = input$table_rows_selected
      #Door het selecteren verandert de plot...

      if (length(s)) {
        subdataset <- clickInfo()[s, , drop = FALSE]

        colBool2 <- subdataset$qval<input$alpha
        colors2 <- rep(NA,length(subdataset$qval))
        colors2[colBool2] <- "purple"
        colors2[!colBool2] <- "darkgrey"

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

  #for zooming
  clickInfo <- reactive({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
    if(!is.null(ranges$x) && !is.null(ranges$y)){clickInfo <- subset(dataset(), (dataset()$estimate>ranges$x[1] & dataset()$estimate<ranges$x[2] & dataset()$minus_log10_p>ranges$y[1] & dataset()$minus_log10_p<ranges$y[2]))
    } else if(is.null(ranges$x) && is.null(ranges$y)){clickInfo <- dataset()}
    return(clickInfo)
  })

  ###Generation of datatable for output###
  #Data table
  data <- reactive(
    {
      data <- clickInfo()
      oldnames <- c("se","df","Tval","pval","qval","signif","pvalS1","qvalS1","signifS1","AveExpr","df_num","df_den","Fval")
      data$signif=data$qval<input$alpha

      newnames <- c("standard error","degrees of freedom","T value","p value", "false discovery rate","significant","p value stage 1","false discovery rate stage 1","significant stage 1","average expression","degrees of freedom numerator","degrees of freedom denominator","F value")
      for(i in 1:length(oldnames)){
        colnames(data)[colnames(data)==oldnames[i]] <- newnames[i]
      }
      return(data)
    }
  )

  output$table<-DT::renderDataTable(
    data()
  )

  #Set table Proxy so as to reduce the table according to the zoom in the plot and to highlight points
  proxy = dataTableProxy('table')
  #Add and remove points by clicking in the plot window
  observeEvent(input$plot1_click, {

    selected <- nearPoints(clickInfo(), input$plot1_click, addDist = TRUE,maxpoints=1, xvar=estimate(), yvar="minus_log10_p")
    sel_rows <- which(clickInfo()$Accessions %in% selected$Accessions)
    #Rows which were selected and selected again are removed, rows which were already selected but not selected again are retained
    #Don't sort this! Otherwise reacalculated.
    new_rows <- c(sel_rows[!sel_rows%in%input$table_rows_selected], input$table_rows_selected[!input$table_rows_selected%in%sel_rows])

    proxy %>% DT::selectRows(new_rows)
  })

  #Enable or disable add brush to selection and remove brush from selection buttons
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

  observeEvent(input$remove_all_selection, {
    proxy %>% DT::selectRows(NULL)
  })


  ###Detail Plot###
  #Drop down menu for plot 2
  plot2DependentVars <- reactive({
    if(input$save==2 & !is.null(input$load_model$datapath)){
      RData <- try(loads_MSqRob(file=modelDatapath(), variables=c("fixed","random")), silent=TRUE)
      if(inherits(RData, 'try-error')){stop("Loading of model file failed. Please provide a valid RDatas model file.")}
      plot2DependentVars <- as.list(c(RData$fixed, RData$random))
    } else{
      plot2DependentVars <- as.list(c(input$fixed, input$random))
    }
    return(plot2DependentVars)
  })

  plot2OtherVars <- reactive({
    c("none",plot2DependentVars())
  })

  plot2MainVars <- reactiveValues(
    values=NULL
  )

  observeEvent(input$go, {
    plot2MainVars$values <- names(data())
  })

  output$selectMainPlot2 <- renderUI({
    selectInput("selMainPlot2", "Select title variable", plot2MainVars$values)
  })

  output$selectPlot2 <- renderUI({
    selectInput("selPlot2", "Select independent variable", plot2DependentVars())
  })

  output$selectColPlot2 <- renderUI({
    selectInput("selColPlot2", "Select color variable", plot2OtherVars())
  })

  output$selectPchPlot2 <- renderUI({
    selectInput("selPchPlot2", "Select shape variable", plot2OtherVars())
  })

  acc_plot2 <- reactive({as.character(clickInfo()[input$table_rows_selected,"Accessions"])})
  #indep_var_plot2 <- reactive({input$selPlot2})
  #color_var_plot2 <- reactive({input$selColPlot2})


  colorsPlot2 <- reactive({
    accessions <- acc_plot2()
    proteins <- outputlist()$RData$proteins
    #indep_var <- indep_var_plot2()
    #color_var <- color_var_plot2()
    if(input$selColPlot2=="none"){colors <- 1} else{
      colordata <- tryCatch(getData(proteins[accessions])[,input$selColPlot2], error=function(e){
        return(NULL)
      })
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(length(unique(colordata)))
      colors <- colors[as.numeric(droplevels(as.factor(colordata)))]
    }
    return(colors)
  })


  pchPlot2 <- reactive({
    accessions <- acc_plot2()
    proteins <- outputlist()$RData$proteins
    #indep_var <- indep_var_plot2()
    #color_var <- color_var_plot2()
    if(input$selPchPlot2=="none"){pch_vals <- 1} else{
      pchdata <- tryCatch(getData(proteins[accessions])[,input$selPchPlot2], error=function(e){
        return(NULL)
      })
      pch_vals <- c(0:25,32:127)
      points <- as.numeric(droplevels(as.factor(pchdata)))
      #Repeat pch_vals if there would be more than 122 unique levels
      pch_vals <- rep_len(pch_vals, length(unique(points)))
      pch_vals <- pch_vals[points]
    }
    return(pch_vals)
  })


  #plot 2
  output$plot2 <- renderPlot({
    accessions <- acc_plot2() #Geeft "NA"
    proteins <- outputlist()$RData$proteins
    #indep_var <- indep_var_plot2()
    #color_var <- color_var_plot2()

    if(length(accessions)==1){

      #Needed for "main"
      s = input$table_rows_selected
      subdataset <- data()[s, , drop = FALSE]
      main <- subdataset[,input$selMainPlot2]

      if (is.factor(getData(proteins[accessions])[[input$selPlot2]])){
        boxplot(getData(proteins[accessions])$quant_value~getData(proteins[accessions])[[input$selPlot2]], outline=FALSE, ylim=c(min(getData(proteins[accessions])$quant_value)-0.2,max(getData(proteins[accessions])$quant_value)+0.2), ylab="preprocessed peptide intensity", xlab="", main=main, las=2, frame.plot=FALSE, frame=FALSE, col="grey", pars=list(boxcol="white")) #, cex.main=2, cex.lab=2, cex.axis=2, cex=2, getAccessions(proteins[accessions])
        points(jitter((as.numeric(getData(proteins[accessions])[[input$selPlot2]])), factor=2),getData(proteins[accessions])$quant_value, col=colorsPlot2(), pch=pchPlot2()) #,cex=2, lwd=2, col=c(1,2,3,4,"cyan2",6)
        # title(ylab="Log2(Intensity)", line=5, cex.lab=2, family="Calibri Light")
      } else plot(getData(proteins[accessions])$quant_value~getData(proteins[accessions])[[input$selPlot2]],ylim=c(min(getData(proteins[accessions])$quant_value)-0.2,max(getData(proteins[accessions])$quant_value)+0.2), ylab="preprocessed peptide intensity", xlab="", main=getAccessions(proteins[accessions]), las=2, bty="n")
    } else{NULL}
  })

  output$nText <- renderText({
    outputlist()$test
  })

  ##############################################
  #Normalization tab
  #############################################
  ###Function plotDens see utilities.R

  ###Drop down menu for plot normalization Plot###
  plotNorm1DependentVars <- reactive({
    as.list(c("none",colnames(fixedOptions2())))
  })

  output$selectColPlotNorm1 <- renderUI({
    selectInput("selColPlotNorm1", "Select color variable",  plotNorm1DependentVars())
  })

  output$npeptidesNormalized = renderText(NULL)

  ####Raw peptide density plot####
  colorsNorm <- reactive({
    colors <- 1
    try(
      {colordata <- fixedOptions2()[,input$selColPlotNorm1]
      colors<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(length(unique(colordata)))
      colors <- colors[as.numeric(droplevels(as.factor(colordata)))]
      },silent=TRUE)
    return(colors)
  })





  ####NEW######

  peps <- reactive({
    if(!is.null(peptidesDatapath())){
      read_MaxQuant(peptidesDatapath(), shiny=TRUE, message="Importing data...")
    } else{NULL}
  })

  esetN <- reactive({

    if(!is.null(peps()) & isTRUE(input$evalnorm)){
      #If remove only identified by site==TRUE and fileProteinGroups is NULL, this also throws an error
      pepsN <- preprocess_MaxQuant(peps(), logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, filter=gsub(" ",".",input$filter), remove_only_site=input$onlysite, file_proteinGroups=proteinGroupsDatapath(), filter_symbol="+", minIdentified=input$minIdentified, shiny=TRUE, printProgress=TRUE, message="Preprocessing data...")
      esetN <- exprs(pepsN)
    }

  })

  getDensXlimYim <- function(eset){
    densAll=apply(eset,2,density,na.rm=TRUE)
    ymax=max(sapply(densAll,function(d) max(d$y)))
    xlim=range(eset,na.rm=TRUE)
    ylim=c(0,ymax)
    return(list(densAll=densAll, xlim=xlim, ylim=ylim))
  }

  eset <- reactive({

    if(!is.null(peptidesDatapath())){

      if(isTRUE(input$logtransform)) {eset <- log(exprs(peps()),base=input$log_base)

      } else {eset <- exprs(peps())}
      eset[is.infinite(eset)] <- NA
    } else{eset <- NULL}

    return(eset)
  })

  output$plotRaw<- renderPlot({
    if(isTRUE(input$onlysite) && is.null(input$proteingroups)){stop("Please provide a protein groups file or untick the box \"Remove proteins that are only identified by modified peptides\".")}
    if(!is.null(eset())){
      densXlimYim <- getDensXlimYim(eset())
      plotDens(eset(), densXlimYim[["densAll"]], densXlimYim[["xlim"]], densXlimYim[["ylim"]], colorsNorm(), main="")
      output$npeptidesRaw = renderText(nrow(peps()))
    }
  })


  output$plotNorm1<- renderPlot({
    if(isTRUE(input$onlysite) && is.null(input$proteingroups)){stop("Please provide a protein groups file or untick the box \"Remove proteins that are only identified by modified peptides\".")}
    if(isTRUE(input$evalnorm) & !is.null(esetN())){
      densXlimYimN <- getDensXlimYim(esetN())
      plotDens(esetN(), densXlimYimN[["densAll"]], densXlimYimN[["xlim"]], densXlimYimN[["ylim"]], colorsNorm(), main="")
      output$npeptidesNormalized = renderText(nrow(esetN()))
    }
  })

  ####NEW######







  # output$plotRaw<- renderPlot({
  #
  #  if(input$onlysite && is.null(input$proteingroups)){stop("Please provide a protein groups file or untick the box \"Remove proteins that are only identified by modified peptides\".")}
  #
  #  if(!is.null(peptidesDatapath())){
  #       peps=read_MaxQuant(peptidesDatapath(), shiny=TRUE)
  #       if(input$logtransform) eset=log(exprs(peps),base=input$log_base) else eset=exprs(peps)
  #       eset[is.infinite(eset)]=NA
  #
  #       densAll=apply(eset,2,density,na.rm=TRUE)
  #       ymax=max(sapply(densAll,function(d) max(d$y)))
  #       xlim=range(eset,na.rm=TRUE)
  #       ylim=c(0,ymax)
  #
  #        plotDens(eset, densAll, xlim, ylim, colorsNorm(), main="")
  #        output$npeptidesRaw = renderText(nrow(peps))
  #       }
  # })
  #
  #
  # ####Normalized peptide density plot####
  # output$plotNorm1<- renderPlot({
  #
  # if(input$onlysite && is.null(input$proteingroups)){stop("Please provide a protein groups file or untick the box \"Remove proteins that are only identified by modified peptides\".")}
  #
  #
  #       if(!is.null(peptidesDatapath())){
  #       peps=read_MaxQuant(peptidesDatapath(), shiny=TRUE)
  #       if (input$evalnorm)
  #       error=try({
  #
  #           filter <- gsub(" ",".",input$filter)
  #
  #           pepsN=preprocess_MaxQuant(peps, logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, filter=filter, remove_only_site=input$onlysite, file_proteinGroups=proteinGroupsDatapath(), filter_symbol="+", minIdentified=input$minIdentified, shiny=TRUE, printProgress=TRUE, message="Preprocessing data...")
  #
  #           esetN=exprs(pepsN)
  #
  #           densAllN=apply(esetN,2,density,na.rm=TRUE)
  #           ymax=max(sapply(densAllN,function(d) max(d$y)))
  #           xlim=range(esetN,na.rm=TRUE)
  #           ylim=c(0,ymax)
  #
  #        plotDens(esetN, densAllN, xlim, ylim, colorsNorm(), main="")
  #        output$npeptidesNormalized = renderText(nrow(pepsN))
  #
  #       },silent=TRUE)
  #       }
  # })

  ranges2 <- reactiveValues(x = NULL, y = NULL)

  ####MDS plot with zoom####
  observeEvent(input$plotMDS_dblclick, {
    brush <- input$plotMDS_brush
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      ranges2$y <- c(brush$ymin, brush$ymax)

    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })

  output$plotMDS<- renderPlot(
    {
      if(isTRUE(input$onlysite) && is.null(input$proteingroups)){stop("Please provide a protein groups file or untick the box \"Remove proteins that are only identified by modified peptides\".")}

      if(isTRUE(input$evalnorm) & !is.null(esetN())){

        if(isTRUE(input$plotMDSLabels) & !isTRUE(input$plotMDSPoints)){
          #Only labels
          limma::plotMDS(esetN(), col=colorsNorm(), xlim = ranges2$x, ylim = ranges2$y, las=1, bty="n")
        } else{

          mds <- plotMDS(esetN(), plot=FALSE)

          #Dots with labels
          if(isTRUE(input$plotMDSLabels) & isTRUE(input$plotMDSPoints)){

            #Need extra space on top for labels
            if(is.null(ranges2$y)){
              yrange=c(min(mds$y),(max(mds$y)+(range(mds$y)[2]-range(mds$y)[1])/10))
            } else{
              yrange=c(ranges2$y[1],(ranges2$y[2]+(ranges2$y[2]-ranges2$y[1])/10))
            }

            plot(mds, col=colorsNorm(), xlim = ranges2$x, ylim = yrange, las=1, bty="n", xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
            text(mds, labels=colnames(esetN()), col=colorsNorm(), cex= 1, pos=3)
          }

          #Only dots
          if(!isTRUE(input$plotMDSLabels) & isTRUE(input$plotMDSPoints)){
            plot(mds, col=colorsNorm(), xlim = ranges2$x, ylim = ranges2$y, las=1, bty="n", xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
          }

          #No dots and no labels => no plot!

        }
        #pch NULL, pch NA, text
        # plotMDS(mds, las=1, bty="n",pch=NULL)
        # text(mds, labels=colnames(test), cex= 0.7, pos=3)

      }
    })

})
