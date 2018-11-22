#source('directoryInput.R')
library(shiny)
library(DT)
library(shinyjs)
library(lme4)
library(MSqRob)
library(MSnbase)
library(grDevices)
library(limma)
library(graphics)
source("utilities.R")

#Max file size: 50 TB
options(shiny.maxRequestSize=50000*1024^3)
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  ###########################################
  #Input Tab
  ###########################################

  #Change the default for protein groups: if MaxQuant: set to true, if not, set to FALSE
  observe({
    if(input$input_type == "MaxQuant"){
      value <- TRUE
    } else{
      value <- FALSE
    }

    updateCheckboxInput(session, "onlysite", value = value)
  })

  saveFolder <- reactiveValues(folder = NULL)

  #Function to check if save folder is set
  check_save_folder <- function(save_folder){
    if(is.null(save_folder) | length(save_folder)==0){stop("No output location selected!")
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
    } else if(Sys.info()["sysname"]=="Linux"){
      saveFolder$folder <- tryCatch(
        choose.dir_Linux()
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

    # folderInput(inputId="outputFolder", label=span("Specify the location where your output will be saved", popify(bsButton("q_outputFolder", label = "", icon = icon("question"), style = "info", size = "extra-small")
    #                                                , title = "Output directory",
    #                                                content = "A new folder will be created in this directory. This folder will contain all your output. This includes your resulting Excel file as well",
    #                                                placement = "right")), value = value, multiple = FALSE, accept = NULL, width = NULL, style=style)
    folderInput(inputId="outputFolder", label=NULL, value = value, multiple = FALSE, accept = NULL, width = NULL, style=style)

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
      if(input$input_type=="Progenesis"){
        test <- read.table(peptidesDatapath(), sep=",", skip=2, nrows=1, quote="", comment.char = "", stringsAsFactors = FALSE)
        a <- gsub("^\\\"","",test)
        b <- gsub("\\\"$","",a)
        make.names(b)
      } else{
        make.names(as.vector(as.matrix(read.table(peptidesDatapath(), nrows=1, sep="\t", quote="", comment.char = ""))))
      }
    }
  })
  selectedFilter <- reactive({
    if(!any(c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions())) {
      NULL
    } else{c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant")[c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions()]}
  })

  output$selectFilters <- renderUI({
    selectInput("filter", NULL, filterOptions(), multiple=TRUE, selected=selectedFilter(), width = '100%')})

  selNormalisation <- reactive({
    if(input$input_type=="Progenesis"){
      selNormalisation <- "none"
    } else{
      selNormalisation <- NULL
    }
  })

  output$selectNormalisation <- renderUI({
    selectInput("normalisation", NULL, c("loess.fast", "rlr", "quantiles", "quantiles.robust", "vsn", "center.median", "center.mean", "max", "sum", "none"), selected=selNormalisation(), width = '100%') #"loess.affy" and "loess.pairs" left out on purpose because they remove everything with at least 1 NA!
    })


  ########################################
  #Generate options for fixed effect variables
  ########################################

  colClasses <- reactive({
    if(isTRUE(input$asis_numeric)){
      colClasses <- "keep"
    } else{colClasses <- "factor"}
    return(colClasses)
  })

  #from annotation file
  exp_annotation <- reactive({
    if(is.null(input$annotation$name) || is.null(eset())){ #Needs to be "eset()" here and not "eset3()", "evaluation nested too deeply" error because exp_annotation() depends both on exp_annotation and peptide input
      exp_annotation <- NULL
    } else{
      if(isTRUE(as.logical(grep(".xlsx[/\\]*$",input$annotation$name)))){

        #Dirty fix for Windows 7, 64 bits problem with rJava
        # if(Sys.info()['sysname']=="Windows"){
        #   test <- xlsx::read.xlsx(file=annotationDatapath(), sheetIndex = 1)
        #   #Remove columns with all NA
        #   test <- test[,!colSums(apply(test, 2, is.na))==nrow(test),drop=FALSE]
        #   return(test)
        # } else{
        exp_annotation <- openxlsx::read.xlsx(annotationDatapath())
        #Convert characters to factors:
        exp_annotation <- as.data.frame(unclass(exp_annotation))
        # }

      } else{
        exp_annotation <- read.table(annotationDatapath(), sep="\t", header=TRUE, row.names = NULL, quote="", stringsAsFactors = TRUE)
      }
      exp_annotation <- makeAnnotation(exp_annotation, run_names=colnames(eset()), colClasses=colClasses())
    }
    return(exp_annotation)
  })
  #from peptides file (filterOptions)
  fixedOptions <- reactive({
    c(as.vector(colnames(exp_annotation())),filterOptions())
  })
  #Generate option of factor levels
  levelOptions <- reactive({
    if((is.null(exp_annotation()) | is.null(input$fixed)) & (input$save!=2 | is.null(input$load_model$datapath))){
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

      optionsFixedSelected <- exp_annotation()[,input$fixed,drop=FALSE]

      #Should stay "sapply" because result can sometimes be double or character
      levelOptions <- unique(unlist(lapply(colnames(optionsFixedSelected),function(name){
        if(is.numeric(optionsFixedSelected[,name])){result <- name
        } else{
          result <- paste0(name,optionsFixedSelected[,name])
        }
        return(result)
      })))

      return(levelOptions)
    }
  })

  ###########################################
  #Functionalities for Quantification Tab
  ###########################################
  nmsFixedOptions <- reactive({names(exp_annotation())})

  ####select Fixed effects, random effects, Proteins and store options ####
  output$selectFixed <- renderUI({
    selectInput("fixed", label=NULL, choices=nmsFixedOptions(), multiple=TRUE )
    })

  selectedRandom <- reactive({
    if(!any(c("Sequence","sequence","peptide") %in% filterOptions())) {
      NULL
    }else{
      c("Sequence","sequence","peptide")[which(c("Sequence","sequence","peptide") %in% filterOptions())[1]]
    }
  })

  output$selectRandom <- renderUI({
      selectInput("random", label=NULL, choices=fixedOptions(), multiple=TRUE, selected=selectedRandom() )
    })

  #   renderText({
  #   verbatimTextOutput(fixedOptions())
  # })

  #Check waarom er nog altijd "Select accession" staat en niet "prot"!!!!
  selectedProteins <- reactive({
    if(!any(c("Proteins","Protein", "prot","Accession","accession") %in% filterOptions())) {
      ""
    }else{
      c("Proteins","Protein","prot","Accession","accession")[which(c("Proteins","Protein","prot","Accession","accession") %in% filterOptions())[1]]
      }
  })

  output$selectProteins <- renderUI({
    selectInput("proteins", label=NULL, filterOptions(), multiple=FALSE, selected=selectedProteins() )
  #!Be careful with "selectizeInput" -> fixed and random all of a sudden might not work anymore...
  })


  selectedAnnotations <- reactive({
    if(!any(c("Gene names", "Protein names","Gene.names", "Protein.names", "gene") %in% filterOptions())) {
      NULL
    }else{c("Gene names", "Protein names","Gene.names", "Protein.names", "gene")[c("Gene names", "Protein names","Gene.names","Protein.names","gene") %in% filterOptions()]}
  })

  output$selectAnnotations <- renderUI({
    selectInput("annotations", label=NULL, filterOptions(), multiple=TRUE, selected=selectedAnnotations() )
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
  shinyjs::onclick("button_outputFolder",
                   shinyjs::toggle(id = "tooltip_outputFolder", anim = TRUE))
  shinyjs::onclick("button_project_name",
                   shinyjs::toggle(id = "tooltip_project_name", anim = TRUE))
  shinyjs::onclick("button_input_type",
                   shinyjs::toggle(id = "tooltip_input_type", anim = TRUE))
  shinyjs::onclick("button_peptides",
                   shinyjs::toggle(id = "tooltip_peptides", anim = TRUE))
  shinyjs::onclick("button_annotation",
                   shinyjs::toggle(id = "tooltip_annotation", anim = TRUE))
  shinyjs::onclick("button_asis_numeric",
                   shinyjs::toggle(id = "tooltip_asis_numeric", anim = TRUE))
  shinyjs::onclick("button_newExpAnnText",
                   shinyjs::toggle(id = "tooltip_newExpAnnText", anim = TRUE))
  shinyjs::onclick("button_cite",
                   shinyjs::toggle(id = "tooltip_cite", anim = TRUE))
  shinyjs::onclick("button_notinlist",
                   shinyjs::toggle(id = "tooltip_notinlist", anim = TRUE))
  shinyjs::onclick("button_logtransform",
                   shinyjs::toggle(id = "tooltip_logtransform", anim = TRUE))
  shinyjs::onclick("button_log_base",
                   shinyjs::toggle(id = "tooltip_log_base", anim = TRUE))
  shinyjs::onclick("button_normalisation",
                   shinyjs::toggle(id = "tooltip_normalisation", anim = TRUE))
  shinyjs::onclick("button_onlysite",
                   shinyjs::toggle(id = "tooltip_onlysite", anim = TRUE))
  shinyjs::onclick("button_proteingroups",
                   shinyjs::toggle(id = "tooltip_proteingroups", anim = TRUE))
  shinyjs::onclick("button_smallestUniqueGroups",
                   shinyjs::toggle(id = "tooltip_smallestUniqueGroups", anim = TRUE))
  shinyjs::onclick("button_minIdentified",
                   shinyjs::toggle(id = "tooltip_minIdentified", anim = TRUE))
  shinyjs::onclick("button_filter",
                   shinyjs::toggle(id = "tooltip_filter", anim = TRUE))
  # shinyjs::onclick("button_evalnorm",
  #                  shinyjs::toggle(id = "tooltip_evalnorm", anim = TRUE))
  shinyjs::onclick("button_selColPlotNorm",
                   shinyjs::toggle(id = "tooltip_selColPlotNorm", anim = TRUE))
  shinyjs::onclick("button_preprocessing_extension",
                   shinyjs::toggle(id = "tooltip_preprocessing_extension", anim = TRUE))
  shinyjs::onclick("button_h4_int_transformation",
                   shinyjs::toggle(id = "tooltip_h4_int_transformation", anim = TRUE))
  shinyjs::onclick("button_h4_full_preprocessing",
                   shinyjs::toggle(id = "tooltip_h4_full_preprocessing", anim = TRUE))
  shinyjs::onclick("button_h4_MDS_full_preprocessing",
                   shinyjs::toggle(id = "tooltip_h4_MDS_full_preprocessing", anim = TRUE))

  shinyjs::onclick("button_proteins",
                   shinyjs::toggle(id = "tooltip_proteins", anim = TRUE))
  shinyjs::onclick("button_annotations",
                   shinyjs::toggle(id = "tooltip_annotations", anim = TRUE))
  shinyjs::onclick("button_fixed",
                   shinyjs::toggle(id = "tooltip_fixed", anim = TRUE))
  shinyjs::onclick("button_random",
                   shinyjs::toggle(id = "tooltip_random", anim = TRUE))
  shinyjs::onclick("button_save",
                   shinyjs::toggle(id = "tooltip_save", anim = TRUE))
  shinyjs::onclick("button_load_model",
                   shinyjs::toggle(id = "tooltip_load_model", anim = TRUE))
  shinyjs::onclick("button_analysis_type",
                   shinyjs::toggle(id = "tooltip_analysis_type", anim = TRUE))
  shinyjs::onclick("button_lfc",
                   shinyjs::toggle(id = "tooltip_lfc", anim = TRUE))
  shinyjs::onclick("button_nContr",
                   shinyjs::toggle(id = "tooltip_nContr", anim = TRUE))
  shinyjs::onclick("button_plot_contrast",
                   shinyjs::toggle(id = "tooltip_plot_contrast", anim = TRUE))
  shinyjs::onclick("button_result_extension",
                   shinyjs::toggle(id = "tooltip_result_extension", anim = TRUE))
  shinyjs::onclick("button_h4_volcano_plot",
                   shinyjs::toggle(id = "tooltip_h4_volcano_plot", anim = TRUE))
  shinyjs::onclick("button_h4_detail_plot",
                   shinyjs::toggle(id = "tooltip_h4_detail_plot", anim = TRUE))
  shinyjs::onclick("button_alpha",
                   shinyjs::toggle(id = "tooltip_alpha", anim = TRUE))
  shinyjs::onclick("button_selMainPlot2",
                   shinyjs::toggle(id = "tooltip_selMainPlot2", anim = TRUE))
  shinyjs::onclick("button_selPlot2",
                   shinyjs::toggle(id = "tooltip_selPlot2", anim = TRUE))
  shinyjs::onclick("button_selColPlot2",
                   shinyjs::toggle(id = "tooltip_selColPlot2", anim = TRUE))
  shinyjs::onclick("button_selPchPlot2",
                   shinyjs::toggle(id = "tooltip_selPchPlot2", anim = TRUE))
})

  observe({

    #Observe input$filter so that when going immediately to tab 3, select load model and then go to tab 2 leads also to disabled input$filter field.
    input$filter

    # Change on selection of tab?
    # input$Input
    # input$Quantification
    # input$Preprocessing

    if(input$save!=2){
      shinyjs::enable("peptides")
      shinyjs::enable("annotation")
      shinyjs::enable("asis_numeric")
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

      # shinyjs::enable("evalnorm")
      shinyjs::enable("selColPlotNorm1")
      shinyjs::enable("preprocessing_extension")
      shinyjs::enable("plotMDSPoints")
      shinyjs::enable("plotMDSLabels")

      shinyjs::enable("borrowFixed")
      shinyjs::enable("borrowRandom")

    } else if(input$save==2){
      shinyjs::disable("peptides")
      shinyjs::disable("annotation")
      shinyjs::disable("asis_numeric")
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

      # shinyjs::disable("evalnorm")
      shinyjs::disable("selColPlotNorm1")
      shinyjs::disable("preprocessing_extension")
      shinyjs::disable("plotMDSPoints")
      shinyjs::disable("plotMDSLabels")

      shinyjs::disable("borrowFixed")
      shinyjs::disable("borrowRandom")
    }
  })

  processedvals = reactive({processInput(input)})

  useful_properties = reactive({
    if(is.null(peps)){useful_properties <- NULL
    } else{
      useful_properties <- unique(c(processedvals()[["proteins"]],processedvals()[["annotations"]],input$fixed,input$random)[c(processedvals()[["proteins"]],processedvals()[["annotations"]],input$fixed,input$random) %in% colnames(Biobase::fData(peps()))])
    }
    return(useful_properties)
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

      processedvals = isolate(processedvals())
      useful_properties = isolate(useful_properties())

      # peptides <- isolate(peps())

      # if(input$input_type=="MaxQuant"){
      # peptides = read_MaxQuant(peptidesDatapath(), shiny=TRUE, message="Importing data...")
      # } else if(input$input_type=="moFF"){
      # peptides = read_moFF(peptidesDatapath(), shiny=TRUE, message="Importing data...")
      # } else{
      # peptides = read_MaxQuant(peptidesDatapath(), shiny=TRUE, message="Importing data...")
      # }

      peptides2 <- isolate(pepsN())

      # if(input$input_type=="MaxQuant"){
      #   peptides2 = preprocess_MaxQuant(peptides, accession=processedvals[["proteins"]], exp_annotation=annotationDatapath(), type_annot=processedvals[["type_annot"]], logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, useful_properties=useful_properties, filter=processedvals[["filter"]], remove_only_site=input$onlysite, file_proteinGroups=proteinGroupsDatapath(), colClasses=colClasses(), filter_symbol="+", minIdentified=input$minIdentified, shiny=TRUE, printProgress=TRUE, message="Preprocessing data...")
      # } else if(input$input_type=="moFF"){
      #   peptides2 = preprocess_moFF(peptides, accession=processedvals[["proteins"]], exp_annotation=annotationDatapath(), type_annot=processedvals[["type_annot"]], logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, useful_properties=useful_properties, filter=processedvals[["filter"]], minIdentified=input$minIdentified, shiny=TRUE, colClasses=colClasses(), printProgress=TRUE, message="Preprocessing data...")
      # } else{
      #   peptides2 = preprocess_MaxQuant(peptides, accession=processedvals[["proteins"]], exp_annotation=annotationDatapath(), type_annot=processedvals[["type_annot"]], logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, useful_properties=useful_properties, filter=processedvals[["filter"]], remove_only_site=input$onlysite, file_proteinGroups=proteinGroupsDatapath(), colClasses=colClasses(),  filter_symbol="+", minIdentified=input$minIdentified, shiny=TRUE, printProgress=TRUE, message="Preprocessing data...")
      # }

      # Biobase::fData(peptides2) <- droplevels(Biobase::fData(peptides2)) # -> nu geïmplementeerd in preprocess_MaxQuant en preprocess_moFF!
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
    lfc <- input$lfc

    #If standard
    if(input$analysis_type=="standard"){
      RidgeSqM <- test.contrast_adjust(outputlist$RData$models, L, simplify=FALSE, lfc=lfc, printProgress=TRUE, shiny=TRUE, message_extract="Calculating contrasts...", message_test="Testing contrasts...")

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
    if(input$analysis_type!="ANOVA"){
    #The contrast corresponding to the selected contrast
    L <- outputlist()$L[,input$plot_contrast,drop=FALSE]
    L <- L[L!=0,,drop=FALSE]
    return(L)
    }
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
      div(class="MSqRob_input_container",
          list(
            tags$label("Show contrast", `for`="plot_contrast", class="MSqRob_label"),
            tags$button(id="button_plot_contrast", tags$sup("[?]"), class="MSqRob_tooltip"),
            selectInput("plot_contrast", label=NULL, contrastOptions()),
            hidden(helpText(id="tooltip_plot_contrast","Select the contrast you want to visualize."))
          )
      )
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



  ###Volcano plot###

  makeVolcanoPlot <- function(dataset,estimate,clickInfo,input){
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
  }

  output$plot1 <- renderPlot({
    makeVolcanoPlot(dataset,estimate,clickInfo,input)
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
    if(!is.null(ranges$x) && !is.null(ranges$y)){clickInfo <- subset(dataset(), (dataset()[[estimate()]]>ranges$x[1] & dataset()[[estimate()]]<ranges$x[2] & dataset()$minus_log10_p>ranges$y[1] & dataset()$minus_log10_p<ranges$y[2]))
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
    div(class="MSqRob_input_container",
        list(
          tags$label("Title variable", `for`="selMainPlot2", class="MSqRob_label"),
          tags$button(id="button_selMainPlot2", tags$sup("[?]"), class="MSqRob_tooltip"),
    selectInput("selMainPlot2", label=NULL, plot2MainVars$values),
    hidden(helpText(id="tooltip_selMainPlot2","
                    Select the title variable.
                    This is the title that will be shown on top of the detail plot.
                    "))
        )
    )
  })

  output$selectPlot2 <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Independent variable", `for`="selPlot2", class="MSqRob_label"),
          tags$button(id="button_selPlot2", tags$sup("[?]"), class="MSqRob_tooltip"),
    selectInput("selPlot2", label=NULL, plot2DependentVars()),
    hidden(helpText(id="tooltip_selPlot2","
                    Select the independent variable.
                    This is the variable that will be present on the x-axis of the plot.
                    "))
        )
    )
  })

  output$selectColPlot2 <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Color variable", `for`="selColPlot2", class="MSqRob_label"),
          tags$button(id="button_selColPlot2", tags$sup("[?]"), class="MSqRob_tooltip"),
    selectInput("selColPlot2", label=NULL, plot2OtherVars()),
    hidden(helpText(id="tooltip_selColPlot2","
                    Select the color variable.
                    This is the variable by which the individual peptides should be colored.
                    "))
        )
    )
  })

  output$selectPchPlot2 <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Shape variable", `for`="selPchPlot2", class="MSqRob_label"),
          tags$button(id="button_selPchPlot2", tags$sup("[?]"), class="MSqRob_tooltip"),
          selectInput("selPchPlot2", label=NULL, plot2OtherVars()),
          hidden(helpText(id="tooltip_selPchPlot2","
                          Select the shape variable.
                          This is the variable that is represented by the different plot symbols.
                          "))
        )
    )
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

  ###The detail plot###

  makeDetailPlot <- function(data,acc_plot2,outputlist,input){
    accessions <- acc_plot2() #Geeft "NA"
    proteins <- outputlist()$RData$proteins
    #indep_var <- indep_var_plot2()
    #color_var <- color_var_plot2()

    if(length(accessions)==1){

      #Needed for "main"
      s = input$table_rows_selected
      subdataset <- data()[s, , drop = FALSE]
      main <- subdataset[,input$selMainPlot2]

      #if(is.factor(getData(proteins[accessions])[[input$selPlot2]])){
        boxplot(getData(proteins[accessions])$quant_value~as.factor(getData(proteins[accessions])[[input$selPlot2]]), outline=FALSE, ylim=c(min(getData(proteins[accessions])$quant_value)-0.2,max(getData(proteins[accessions])$quant_value)+0.2), ylab="preprocessed peptide intensity", xlab="", main=main, las=2, frame.plot=FALSE, frame=FALSE, col="grey", pars=list(boxcol="white")) #, cex.main=2, cex.lab=2, cex.axis=2, cex=2, getAccessions(proteins[accessions])
        points(jitter((as.numeric(getData(proteins[accessions])[[input$selPlot2]])), factor=2),getData(proteins[accessions])$quant_value, col=colorsPlot2(), pch=pchPlot2()) #,cex=2, lwd=2, col=c(1,2,3,4,"cyan2",6)
        # title(ylab="Log2(Intensity)", line=5, cex.lab=2, family="Calibri Light")
      #} else plot(getData(proteins[accessions])$quant_value~getData(proteins[accessions])[[input$selPlot2]], outline=FALSE, ylim=c(min(getData(proteins[accessions])$quant_value)-0.2,max(getData(proteins[accessions])$quant_value)+0.2), ylab="preprocessed peptide intensity", xlab="", main=main, las=2, frame.plot=FALSE, frame=FALSE, bty="n")
    } else{NULL}
  }

  #plot 2
  output$plot2 <- renderPlot({
    makeDetailPlot(data,acc_plot2,outputlist,input)
  })

  output$nText <- renderText({
    outputlist()$test
  })

  ########################################################################################
  ###Save the volcano and detail plots when the downloadResultPlots button is clicked#####
  ########################################################################################

  observeEvent(input$downloadResultPlots, {

    savepathFigs <- getDataPath(saveFolder$folder)

    plotnames <- c("volcano","detail")

    for(i in 1:length(plotnames)){

      if(input$result_extension == "svg"){

        grDevices::svg(file.path(savepathFigs,paste0(plotnames[i],".svg")), width=10, height=10,onefile=FALSE) ##onefile=FALSE to prevent the plots that are made previously to appear in the svg

      } else if(input$result_extension == "png"){
        grDevices::png(file.path(savepathFigs,paste0(plotnames[i],".png")),width = 3.25,
            height    = 3.25,
            units     = "in",
            res       = 2000,
            pointsize = 5.4
        )
      } else if(input$result_extension == "pdf") {
        grDevices::cairo_pdf(file.path(savepathFigs,paste0(plotnames[i],".pdf")),width=10,height=10,onefile=FALSE) #onefile=FALSE to prevent the plots that are made previously to appear in the PDF

      } else if(input$result_extension == "tiff") {
        grDevices::tiff(file.path(savepathFigs,paste0(plotnames[i],".tiff")),width = 3.25,
                        height    = 3.25,
                        units     = "in",
                        res       = 2000,
                        pointsize = 5.4)
      } else if(input$result_extension == "jpeg") {
        grDevices::jpeg(file.path(savepathFigs,paste0(plotnames[i],".jpeg")),width = 3.25,
                        height    = 3.25,
                        units     = "in",
                        res       = 2000,
                        pointsize = 5.4)
      } else if(input$result_extension == "bmp") {
        grDevices::bmp(file.path(savepathFigs,paste0(plotnames[i],".bmp")),width = 3.25,
                        height    = 3.25,
                        units     = "in",
                        res       = 2000,
                        pointsize = 5.4)
      } else if(input$result_extension == "postscript") {
        grDevices::postscript(file.path(savepathFigs,paste0(plotnames[i],".ps")),width = 3.25,
                       height    = 3.25,
                       pointsize = 5.4)
      } else if(input$result_extension == "xfig") {
        grDevices::xfig(file.path(savepathFigs,paste0(plotnames[i],".fig")),width = 3.25,
                              height    = 3.25,
                              pointsize = 5.4)
      }

      if(i==1){makeVolcanoPlot(dataset,estimate,clickInfo,input)
      } else if(i==2){makeDetailPlot(data,acc_plot2,outputlist,input)}

      dev.off()  # turn the device off

    }

  })

  #Show a notification when the plots are saved
  observeEvent(input$downloadResultPlots, {
    showNotification("Result plots saved.", duration = 1.5)
  })

  ##############################################
  #Normalization tab
  #############################################
  ###Function plotDens see utilities.R

  ###Drop down menu for plot normalization Plot###
  plotNorm1DependentVars <- reactive({
    as.list(c("none",colnames(exp_annotation())))
  })

  output$selectColPlotNorm1 <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Color variable", `for`="selColPlotNorm", class="MSqRob_label"),
          tags$button(id="button_selColPlotNorm", tags$sup("[?]"), class="MSqRob_tooltip"),
          selectInput("selColPlotNorm1", label=NULL,  plotNorm1DependentVars()),
          hidden(helpText(id="tooltip_selColPlotNorm","Select the variable by which the densities should be colored."))
        )
    )
  })

  output$npeptidesNormalized = renderText(NULL)

  ####Raw peptide density plot####
  colorsNorm <- reactive({
    colors <- 1
    try(
      {colordata <- exp_annotation()[,input$selColPlotNorm1]
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(length(unique(colordata)))
      colors <- colors[as.numeric(droplevels(as.factor(colordata)))]
      },silent=TRUE)
    return(colors)
  })





  ####NEW######

  rawPeptides <- reactive({
    if(!is.null(peptidesDatapath())){

      rawPeptides = import2MSnSet(file=peptidesDatapath(), filetype=input$input_type, shiny=TRUE, message="Importing data...")

    } else{rawPeptides <- NULL}
    return(rawPeptides)
  })

  peps <- reactive({
    if(!is.null(rawPeptides())){

      if(input$input_type=="Progenesis"){
      peps = aggregateMSnSet(rawPeptides(), aggr_by="Sequence", aggr_function="sum", split=", ", shiny=TRUE, printProgress=TRUE, message="Aggregating peptides")
      } else if(input$input_type=="mzTab"){
      peps = aggregateMSnSet(rawPeptides(), aggr_by="sequence", aggr_function="sum", split=", ", shiny=TRUE, printProgress=TRUE, message="Aggregating peptides")
      } else {
      peps = rawPeptides()
      }

    } else{peps <- NULL}

    return(peps)
  })

  pepsN <- reactive({
    if(!is.null(peps()) && !is.null(input$proteins)){ #!is.null(input$proteins) is there to prevent this from double running, e.g. with moFF! #& isTRUE(input$evalnorm)
      #If remove only identified by site==TRUE and fileProteinGroups is NULL, this also throws an error

      if(input$input_type=="MaxQuant"){
        pepsN <- preprocess_MaxQuant(peps(), accession=processedvals()[["proteins"]], exp_annotation=annotationDatapath(), type_annot=processedvals()[["type_annot"]], logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, useful_properties=useful_properties(), filter=processedvals()[["filter"]], remove_only_site=input$onlysite, file_proteinGroups=proteinGroupsDatapath(), colClasses=colClasses(), filter_symbol="+", minIdentified=input$minIdentified, shiny=TRUE, printProgress=TRUE, message="Preprocessing data...")
      } else {
        pepsN <- preprocess_generic(peps(), MSnSetType=input$input_type, exp_annotation=annotationDatapath(), type_annot=processedvals()[["type_annot"]], logtransform=input$logtransform, base=input$log_base, normalisation=input$normalisation, smallestUniqueGroups=input$smallestUniqueGroups, useful_properties=useful_properties(), filter=processedvals()[["filter"]], minIdentified=input$minIdentified, colClasses=colClasses(), shiny=TRUE, printProgress=TRUE, message="Preprocessing data...") #Deze wordt twee keer gelopen!!!
      }
    } else{
      pepsN <- NULL
    }
    return(pepsN)
  })

  esetN <- reactive({
    if(!is.null(pepsN())){
      esetN <- Biobase::exprs(pepsN())
    } else{
      esetN <- NULL
    }
    return(esetN)
  })

  getDensXlimYlim <- function(eset){
    densAll=apply(eset,2,density,na.rm=TRUE)
    ymax=max(vapply(densAll,function(d) max(d$y),1))
    rangematrix <- vapply(densAll,function(d) range(d$x, na.rm=TRUE), c(1,1)) #no longer range(eset), but range of d$x!
    xlim=range(rangematrix,na.rm=TRUE)
    ylim=c(0,ymax)
    return(list(densAll=densAll, xlim=xlim, ylim=ylim))
  }

  eset <- reactive({

    if(!is.null(peps())){

      if(isTRUE(input$logtransform)) {eset <- log(Biobase::exprs(peps()),base=input$log_base)

      } else {eset <- Biobase::exprs(peps())}
      eset[is.infinite(eset)] <- NA

    } else{eset <- NULL}

    return(eset)
  })

  eset3 <- reactive({
    if(is.null(eset()) || is.null(exp_annotation())){
      eset <- NULL
    } else{
      #Sort columns of unprocessed data in exprs by annotation_run column: needed to make correct plots!
      pData <- exp_annotation()
      annotation_run <- getAnnotationRun(pData=pData, run_names=colnames(eset()))
      eset <- eset()[,match(as.character(pData[,annotation_run]), colnames(eset()))]
    }
    return(eset)
  })

  ###Possibilities for zooming
  rangesRaw <- reactiveValues(x = NULL, y = NULL)
  rangesNorm1 <- reactiveValues(x = NULL, y = NULL)
  rangesMDS <- reactiveValues(x = NULL, y = NULL)

  ####Raw density plot with zoom####
  observeEvent(input$plotRaw_dblclick, {
    brush <- input$plotRaw_brush
    if (!is.null(brush)) {
      rangesRaw$x <- c(brush$xmin, brush$xmax)
      rangesRaw$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesRaw$x <- NULL
      rangesRaw$y <- NULL
    }
  })

  ####Normal density plot with zoom####
  observeEvent(input$plotNorm1_dblclick, {
    brush <- input$plotNorm1_brush
    if (!is.null(brush)) {
      rangesNorm1$x <- c(brush$xmin, brush$xmax)
      rangesNorm1$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesNorm1$x <- NULL
      rangesNorm1$y <- NULL
    }
  })

  ####MDS plot with zoom####
  observeEvent(input$plotMDS_dblclick, {
    brush <- input$plotMDS_brush
    if (!is.null(brush)) {
      rangesMDS$x <- c(brush$xmin, brush$xmax)
      rangesMDS$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesMDS$x <- NULL
      rangesMDS$y <- NULL
    }
  })

  makePlotRaw <- function(input,eset,colorsNorm,rangesRaw){

    if(isTRUE(input$onlysite) && is.null(input$proteingroups)){stop("Please provide a proteinGroups.txt file or untick the box \"Remove only identified by site\".")}
    if(!is.null(eset3())){

      densXlimYlim <- getDensXlimYlim(eset3())

      #Allow for zooming
      if(is.null(rangesRaw$y) & is.null(rangesRaw$x)){
        xlim=densXlimYlim[["xlim"]]
        ylim=densXlimYlim[["ylim"]]
      } else{
        xlim=rangesRaw$x
        ylim=rangesRaw$y
      }

      plotDens(eset3(), densXlimYlim[["densAll"]], xlim, ylim, colorsNorm(), main="")
      output$npeptidesRaw = renderText(nrow(peps()))
    }
  }

  makePlotNorm1 <- function(input,esetN,colorsNorm,rangesNorm){
    if(isTRUE(input$onlysite) && is.null(input$proteingroups)){stop("Please provide a proteinGroups.txt file or untick the box \"Remove only identified by site\".")}
    if(!is.null(esetN())){ #isTRUE(input$evalnorm) &

      densXlimYlimN <- getDensXlimYlim(esetN())

      #Allow for zooming
      if(is.null(rangesNorm1$y) & is.null(rangesNorm1$x)){
        xlim=densXlimYlimN[["xlim"]]
        ylim=densXlimYlimN[["ylim"]]
      } else{
        xlim=rangesNorm1$x
        ylim=rangesNorm1$y
      }

      plotDens(esetN(), densXlimYlimN[["densAll"]], xlim, ylim, colorsNorm(), main="")
      output$npeptidesNormalized = renderText(nrow(esetN()))
    }
  }

  makeMDSPlot <- function(input,esetN,colorsNorm,rangesMDS){
    if(isTRUE(input$onlysite) && is.null(input$proteingroups)){stop("Please provide a proteinGroups.txt file or untick the box \"Remove only identified by site\".")}

    if(!is.null(esetN()) & (isTRUE(input$plotMDSLabels) | isTRUE(input$plotMDSPoints))){ #isTRUE(input$evalnorm) & #Last condition: at least one of the boxes labels or points must be ticked, otherwise no plot!

      mds <- plotMDS(esetN(), plot=FALSE)
      labels_mds <- names(mds$x) #Doesn't matter whether you take mds$x or mds$y here

      #Only dots: plot is always made in order to set the ranges correctly: strwidth and strheight are calculated based on the previous plot!
      plot(mds, col=colorsNorm(), xlim = rangesMDS$x, ylim = rangesMDS$y, las=1, bty="n", xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")

      #Calculate maximal increase in plot size due to labels
      max_increase <- (par("usr")[2]-par("usr")[1]+max(graphics::strwidth(labels_mds)))/(par("usr")[2]-par("usr")[1])

      #Need extra space on x axis for labels
      if(is.null(rangesMDS$x)){
        xrange=c((min(mds$x)-max(graphics::strwidth(labels_mds))*max_increase/2),(max(mds$x)+max(graphics::strwidth(labels_mds))*max_increase/2)) #min(mds$x), max(mds$x)
      } else{
        xrange=c(rangesMDS$x[1],rangesMDS$x[2])
      }

      #Only labels
      if(isTRUE(input$plotMDSLabels) & !isTRUE(input$plotMDSPoints)){
        limma::plotMDS(esetN(), col=colorsNorm(), xlim = xrange, ylim = rangesMDS$y, las=1, bty="n")
      }

      #Labels and dots
      if(isTRUE(input$plotMDSLabels) & isTRUE(input$plotMDSPoints)){

        #Need extra space on top for labels
        if(is.null(rangesMDS$y)){
          yrange=c(min(mds$y),(max(mds$y)+max(graphics::strheight(labels_mds)))) #c(min(mds$y),(max(mds$y)+(range(mds$y)[2]-range(mds$y)[1])/10))
        } else{
          yrange=c(rangesMDS$y[1],rangesMDS$y[2]) #c(rangesMDS$y[1],(rangesMDS$y[2]+(rangesMDS$y[2]-rangesMDS$y[1])/10))
        }

        plot(mds, col=colorsNorm(), xlim = xrange, ylim = yrange, las=1, bty="n", xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
        text(mds, labels=colnames(esetN()), col=colorsNorm(), cex= 1, pos=3)
      }

      #pch NULL, pch NA, text
      # plotMDS(mds, las=1, bty="n",pch=NULL)
      # text(mds, labels=colnames(test), cex= 0.7, pos=3)

    }
  }

  #Render the preprocessing plots
  output$plotRaw<- renderPlot({
    makePlotRaw(input,eset,colorsNorm,rangesRaw)
  })

  output$plotNorm1<- renderPlot({
    makePlotNorm1(input,esetN,colorsNorm,rangesNorm)
  })

  output$plotMDS <- renderPlot({
    makeMDSPlot(input,esetN,colorsNorm,rangesMDS)
  })

  #Show a notification when the plots are saved
  observeEvent(input$downloadPreprocessingPlots, {
    showNotification("Diagnostic plots saved.", duration = 1.5)
  })

  #Save the preprocessing plots when the downloadPreprocessingPlots button is clicked
  observeEvent(input$downloadPreprocessingPlots, {

    savepathFigs <- getDataPath(saveFolder$folder)

    plotnames <- c("rawDensity","preprocessedDensity","MDSPlot")

    for(i in 1:length(plotnames)){

    if(input$preprocessing_extension == "svg"){

      grDevices::svg(file.path(savepathFigs,paste0(plotnames[i],".svg")), width=10, height=10) #filenames,i

    } else if(input$preprocessing_extension == "png"){
      grDevices::png(file.path(savepathFigs,paste0(plotnames[i],".png")),width = 3.25,
         height    = 3.25,
         units     = "in",
         res       = 2000,
         pointsize = 5.4
      )
    }
    else if(input$preprocessing_extension == "pdf") {
      grDevices::cairo_pdf(file.path(savepathFigs,paste0(plotnames[i],".pdf")),width=10,height=10,onefile=FALSE) #onefile=FALSE to prevent the plots that are made previously to appear in the PDF
    } else if(input$preprocessing_extension == "tiff") {
      grDevices::tiff(file.path(savepathFigs,paste0(plotnames[i],".tiff")),width = 3.25,
                      height    = 3.25,
                      units     = "in",
                      res       = 2000,
                      pointsize = 5.4)
    } else if(input$preprocessing_extension == "jpeg") {
      grDevices::jpeg(file.path(savepathFigs,paste0(plotnames[i],".jpeg")),width = 3.25,
                      height    = 3.25,
                      units     = "in",
                      res       = 2000,
                      pointsize = 5.4)
    } else if(input$preprocessing_extension == "bmp") {
      grDevices::bmp(file.path(savepathFigs,paste0(plotnames[i],".bmp")),width = 3.25,
                      height    = 3.25,
                      units     = "in",
                      res       = 2000,
                      pointsize = 5.4)
    } else if(input$preprocessing_extension == "postscript") {
      grDevices::postscript(file.path(savepathFigs,paste0(plotnames[i],".ps")),width = 3.25,
                     height    = 3.25,
                     pointsize = 5.4)
    } else if(input$preprocessing_extension == "xfig") {
      grDevices::xfig(file.path(savepathFigs,paste0(plotnames[i],".fig")),width = 3.25,
                            height    = 3.25,
                            pointsize = 5.4)
    }

    if(i==1){makePlotRaw(input,eset,colorsNorm,rangesRaw)
    } else if(i==2){makePlotNorm1(input,esetN,colorsNorm,rangesNorm)
    } else if (i==3){makeMDSPlot(input,esetN,colorsNorm,rangesMDS)}

    dev.off()  # turn the device off

  }

  })

  observe({

    #If no save folder specified or the preprocessing plots are not made, do not allow to try to download the preprocessing plots.
    res <- try(check_save_folder(saveFolder$folder),silent = TRUE)

    if(class(res) == "try-error" | (isTRUE(input$onlysite) && is.null(input$proteingroups))){ #!isTRUE(input$evalnorm) |
      shinyjs::disable("downloadPreprocessingPlots")
    } else{
      shinyjs::enable("downloadPreprocessingPlots")
    }
  })

  #Stop the App when closing the browser or ending the session
  session$onSessionEnded(stopApp)
})
