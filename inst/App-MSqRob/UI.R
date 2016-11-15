library(shiny)
library(DT)
library(shinyjs)
source("utilities.R")

#runApp("pathToApp/App-MsqRob-devel")

####################################
###User Interface
####################################

shinyUI(fluidPage(

  #Use shinyjs package
  shinyjs::useShinyjs(),

  #CSS styles
  # tags$head(
  #   tags$style(HTML("
  # [disabled] {
  # color: pink;
  # }
  #   "))
  # ),

############################################################################
#Navigation bar with 3 panel:Input, preprocessing, quantification
############################################################################
   navbarPage("MSqRob for MaxQuant data",


    ####################################
    #input tab
    ####################################

    tabPanel('Input',
    #sidebar
    sidebarLayout(
    #Project name
  	sidebarPanel(textInput("project_name", "Project Name", value = "project", width = NULL, placeholder = NULL),
  	#Folder where everything will be saved
  	htmlOutput("outputFolderOut"),
  	#folderInput(inputId="outputFolder", label="Specify the location where your output will be saved", placeholder = "No folder selected", multiple = FALSE, accept = NULL, width = NULL),
  	#Annotation file
		fileInput(inputId="annotation", label="Specify the location of your experimental annotation file", multiple = FALSE, accept = NULL, width = NULL),
		#Peptides.txt file
		fileInput(inputId="peptides", label="Specify the location of your peptides.txt file", multiple = FALSE, accept = NULL, width = NULL)
     	),
     	#mainpanel
     	mainPanel()
     ),
		#Main panel with number of output and plots
		mainPanel(width = 5,
		          htmlOutput("folderError")
		          )
    )

    ############################
    #Preprocessing tab
    ###########################
    ,tabPanel('Preprocessing',
    sidebarLayout(
	#Sidebar with input
        sidebarPanel(
      	h4("Filtering"),
      	#Filter on peptides only modified by site
        checkboxInput("onlysite", "Remove proteins that are only identified by modified peptides", value=TRUE),
      	conditionalPanel(
        	condition = "input.onlysite == true",
        	fileInput(inputId="proteingroups", label="Specify the location of your proteinGroups.txt file", multiple = FALSE, accept = NULL, width = NULL)
      	),
      	#Filter on peptides number of occurances
      	numericInput("minIdentified", "Minimal number of times a peptide sequence should be identified", value=2, min = 1, max = NA, step = 1, width = NULL),
      htmlOutput("selectFilters"),

	h4("Normalization"),
 	checkboxInput("logtransform", "Log-transform data", value=TRUE),
	conditionalPanel(
        	condition = "input.logtransform == true",
        	numericInput("log_base", "Base", value=2, min = 0, max = NA, step = NA, width = NULL)
      		),
  selectInput("normalisation", "Normalisation", c("quantiles", "quantiles.robust", "vsn", "center.median", "center.mean", "max", "sum", "none"))
	),
	#Main panel with number of output and plots
        mainPanel(width = 5,
 	      checkboxInput("evalnorm", "Evaluate Normalization", value=FALSE),
        strong('Number of peptides before normalization:'),textOutput('npeptidesRaw',container = span),div(),
        strong('Number of peptides after preprocessing:'),textOutput('npeptidesNormalized',container = span),div(),
        plotOutput('plotRaw'),
 	      htmlOutput("selectColPlotNorm1"),
        h4("Normalized intensities"),
        plotOutput('plotNorm1'),
        h4("MDS-plot after normalization"),
        checkboxInput("plotMDSPoints", "Plot MDS points", value=FALSE),
        checkboxInput("plotMDSLabels", "Plot MDS labels", value=TRUE),
        plotOutput('plotMDS',
				click = "plotMDS_click",
                                 dblclick = "plotMDS_dblclick",
                                 brush = brushOpts(
                                  id = "plotMDS_brush",
                                   resetOnNew = TRUE
                                 )
        ),
      	p("Brush and double-click on the selected area to zoom in. Double click outside the selected area to zoom out.")
        )
       )
)

    ###########################
    #Quantification tab
    ###########################
    ,tabPanel('Quantification',

     # Sidebar with model specification
     sidebarLayout(
     	sidebarPanel(
	htmlOutput("selectProteins"),
	htmlOutput("selectAnnotations"),
	htmlOutput("selectFixed"),
	htmlOutput("selectRandom"),
	checkboxInput("borrowFixed", "Borrow information across fixed effects", value = FALSE, width = NULL),
	checkboxInput("borrowRandom", "Borrow information across random effects", value = FALSE, width = NULL),
	radioButtons("save", "Save/load options:",
                   c("Save the models" = 1,  #Conditioneel stukje toevoegen bij downloadHandler dat hij ook model.RData downloadt
                     "Load existing models" = 2,
                     "Don't save the models" = 3
                      )),
	#load saved model
        conditionalPanel(
        	condition = "input.save == 2",
        	fileInput(inputId="load_model", label="Specify the location of your saved model file", multiple = FALSE, accept = NULL, width = NULL)
      	),
	#Type of analysis
     	#selectInput("analysis_type", "Select the type of analysis", c("standard", "stagewise", "ANOVA")),
	selectInput("analysis_type", "Select the type of analysis", c("standard")),

      	numericInput("nContr", "Number of contrasts you want to test", value=1, min = 1, max = NA, step = 1, width = NULL),
	#Specification of contrasts
      	htmlOutput("selectLevels"),
	#Run button
	actionButton("go", "Go") #,
	#Download button
      	# htmlOutput("download_button")
    ),

    #Main panel with results and plots
	mainPanel(
	verbatimTextOutput("nText"),
	htmlOutput("plot_contrast"),
	verbatimTextOutput("contrastL"),

	#Volcano plot
       	fluidRow(
        column(width = 6,  #6 out of 12 => half the screen!
        	h4("Volcano plot"),
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
Choose significance level to visualize features with an FDR level below alpha.
                  "),
		actionButton("add_area_selection", "Add selected area to selection"),
                actionButton("remove_area_selection", "Remove selected area from selection"),
              	numericInput("alpha", "Significance level (alpha)", value=.05, min = 0, max = 1, step = 0.01, width = NULL)
	),
        #Detail plot
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
		            htmlOutput("selectMainPlot2"),
		            htmlOutput("selectPlot2"),
                htmlOutput("selectColPlot2"),
		            htmlOutput("selectPchPlot2")
         	)
         ),

	fluidRow(column(width = 12, h4("Results table"),DT::dataTableOutput('table')))
     	)
    )
)
#close navbar, page, etc.
)))


