library(shiny)
library(DT)
library(shinyjs)
#source('directoryInput.R')

#runApp("/Users/lgoeminn/Documents/phD/App-MSqRob")

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  #Use shinyjs package
  shinyjs::useShinyjs(),

  # Application title
 # titlePanel("MSqRob for MaxQuant data"),
   navbarPage("MSqRob for MaxQuant data",tabPanel('Input',
    sidebarLayout(
    sidebarPanel(

  #mainPanel(
      # sliderInput("bins",
      #             "Number of bins:",
      #             min = 5,
      #             max = 50,
      #             value = 30),

      textInput("project_name", "Project Name", value = "project", width = NULL, placeholder = NULL),

      fileInput(inputId="annotation", label="Specify the location of your experimental annotation file", multiple = FALSE, accept = NULL, width = NULL),
      fileInput(inputId="peptides", label="Specify the location of your peptides.txt file", multiple = FALSE, accept = NULL, width = NULL),


      #Goed!
      conditionalPanel(
        condition = "input.onlysite == true",
        fileInput(inputId="proteingroups", label="Specify the location of your proteinGroups.txt file", multiple = FALSE, accept = NULL, width = NULL)
      )
),
      mainPanel(
      # conditionalPanel(
      #   condition = "output.rows",
      #   selectInput("filter", "Filter out these columns",
      #               "output.rows", multiple=TRUE)),
      h4("Filtering"),
      checkboxInput("onlysite", "Remove proteins that are only identified by modified peptides", value=TRUE),

      numericInput("minIdentified", "Minimal number of times a peptide sequence should be identified", value=2, min = 1, max = NA, step = 1, width = NULL),
      htmlOutput("selectFilters")

      #!!!Normalisation "none" implementeren!!!

	)
	)
)

############################
#Preprocessing
###########################
,tabPanel('Preprocessing',
        sidebarLayout(
        sidebarPanel(
 	checkboxInput("logtransform", "Log-transform data", value=TRUE),

      	conditionalPanel(
        condition = "input.logtransform == true",
        numericInput("log_base", "Base", value=2, min = 0, max = NA, step = NA, width = NULL)
      	),
        h4(" \n"),
     	 selectInput("normalisation", "Normalisation", c("quantiles", "quantiles.robust", "vsn", "center.median", "center.mean", "max", "sum", "none")),
                htmlOutput("selectColPlotNorm1")
	),
        mainPanel(width = 5,h4("Evaluate normalization"),
        plotOutput('plotnorm1')
        )
       )
)

###########################
#Quantification
###########################
,tabPanel('Quantification',

  # Sidebar for model specification
  sidebarLayout(
    sidebarPanel(

      htmlOutput("selectProteins"),

      htmlOutput("selectAnnotations"),

      htmlOutput("selectFixed"),

      htmlOutput("selectRandom"),

      radioButtons("save", "Save/load options:",
                   c("Download the model" = 1,  #Conditioneel stukje toevoegen bij downloadHandler dat hij ook model.RData downloadt
                     "Upload an existing model" = 2#,
                     #"Don't save the model" = 3
                      )),

      # conditionalPanel(
      #   condition = "input.save == 1",
      #   directoryInput('directory', label = 'Select the directory where the model will be saved', value = '~')
      # ),

      conditionalPanel(
        condition = "input.save == 2",
        fileInput(inputId="load_model", label="Specify the location of your saved model file", multiple = FALSE, accept = NULL, width = NULL)
      ),

      selectInput("analysis_type", "Select the type of analysis", c("standard", "stagewise", "ANOVA")),

      numericInput("nContr", "Number of contrasts you want to test", value=1, min = 1, max = NA, step = 1, width = NULL),

      htmlOutput("selectLevels"),

      #uiOutput("new"),

      #textInput("contrasts", "Specify contrasts of interest", value = "", width = NULL, placeholder = NULL),
      #htmlOutput("specifyContrasts"),


      # Partial example
      # selectInput("dataset", "Dataset", c("diamonds", "rock", "pressure", "cars")),
      # conditionalPanel(
      #   condition = "output.nrows",
      #   checkboxInput("headonly", "Only use first 1000 rows")),
      actionButton("go", "Go"),
      htmlOutput("download_button")

      #downloadButton('downloadData', 'Download')
    ),

    # Show a plot of the generated distribution
     mainPanel(

       verbatimTextOutput("nText"),

       htmlOutput("plot_contrast"),
       #
       #fluidRow(class="well",
       #h4("MDS-Plot"),
       #plotOutput("mdsplot", height = 300)
       #),
       
       #Volcano plot
       fluidRow(
         column(width = 6,  #6 out of 12 => half the screen!
                h4("Volcano plot"),
 #            p("Select and deselect points by clicking on them either in the volcano plot or in the results table.#"),

                plotOutput("plot1", height = 300,
                           click = "plot1_click",
                           dblclick = "plot1_dblclick",
                           brush = brushOpts(
                             id = "plot1_brush",
                             resetOnNew = TRUE
                           )
                ),
               p("Select and deselect points by clicking on them either in the volcano plot or in the results table.
Brush and double-click on the selected area to zoom in. Double click outside the selected area to zoom out.
                  "),
               actionButton("add_area_selection", "Add selected area to selection"),
                actionButton("remove_area_selection", "Remove selected area from selection")
         ),
         column(width = 6, #6 out of 12 => half the screen!
                h4("Detail plot"),
                plotOutput("plot2", height = 300,
                           click = "plot2_click",
                           dblclick = "plot2_dblclick",
                           brush = brushOpts(
                             id = "plot2_brush",
                             resetOnNew = TRUE
                           )
                ),
 p("If only one data point is selected in the results table, this plot shows the individual log2 peptide intensities."),
                htmlOutput("selectPlot2"),
                htmlOutput("selectColPlot2")
         )
         ),

       fluidRow(
         column(width = 12,
                h4("Results table"),
                #verbatimTextOutput("click_info")
                DT::dataTableOutput('table')
         )
       )

       #textOutput("nText")

    #   plotOutput("distPlot"),
    #   verbatimTextOutput("nrows"),
    #   verbatimTextOutput("filterOptions")
     )

  )
)
)))


#source('directoryInput.R')

#runApp("/Users/lgoeminn/Documents/phD/App-MSqRob")

# Define UI for application that draws a histogram
#shinyUI(
#fluidPage(

  #Use shinyjs package
#  shinyjs::useShinyjs(),


  # Application title
#  titlePanel("MSqRob for MaxQuant data"),
#  sidebarLayout(
#    sidebarPanel(
      # sliderInput("bins",
      #             "Number of bins:",
      #             min = 5,
      #             max = 50,
      #             value = 30),

     #mainPanel(

     #  fluidRow(
      #   column(width = 12,
      #          h4("Results table"),
      #          #verbatimTextOutput("click_info")
      #          DT::dataTableOutput('table')
      #   )
      # )
#)
#))
