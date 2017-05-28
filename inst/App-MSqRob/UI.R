library(shiny)
library(DT)
library(shinyjs)
library(shinythemes)
source("utilities.R")

#runApp("pathToApp/App-MsqRob-devel")

####################################
###User Interface
####################################

shinyUI(fluidPage(theme = "MSqRob.css",

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
   navbarPage("MSqRob for MaxQuant data v 0.7.0", inverse=TRUE,


    ####################################
    #input tab
    ####################################

    tabPanel('Input',
    #sidebar
    sidebarLayout(
    #Project name
  	sidebarPanel(
  	  
  	  h3("Settings", class="MSqRob_topheader"),
  	  
  	  div(class="MSqRob_input_container",
  	      list(
  	  tags$label("Project name", `for`="project_name", class="MSqRob_label"),
  	  tags$button(id="button_project_name", tags$sup("[?]"), class="MSqRob_tooltip"),
  	  textInput("project_name", NULL, value = "project", width = '100%', placeholder = NULL),
  	hidden(helpText(id="tooltip_project_name",
  	       "Give your project a meaningful name.
  	       This name will be given to your results folder and the files therein.
  	       A time stamp will be automatically appended to the folder name."))
  	      )
  	  ),

  	#Folder where everything will be saved
  	div(class="MSqRob_input_container",
  	    list(
  	tags$label("Output location", `for`="outputFolder", class="MSqRob_label"),
  	tags$button(id="button_outputFolder", tags$sup("[?]"), class="MSqRob_tooltip"),
  	htmlOutput("outputFolderOut"),
  	#folderInput(inputId="outputFolder", label="Specify the location where your output will be saved", placeholder = "No folder selected", multiple = FALSE, accept = NULL, width = NULL),
  	hidden(helpText(id="tooltip_outputFolder","Specify the folder where your output will be saved."))
  	    )
  	),
  	
  	#Peptides.txt file
  	div(class="MSqRob_input_container",
  	    list(
  	tags$label("peptides.txt file", `for`="peptides", class="MSqRob_label"),
  	tags$button(id="button_peptides", tags$sup("[?]"), class="MSqRob_tooltip"),
  	fileInput(inputId="peptides", label=NULL, multiple = FALSE, accept = NULL, width = NULL),
  	hidden(helpText(id="tooltip_peptides","Specify the location of your peptides.txt file."))
  	    )
  	),
  	
  	#Annotation file
  	div(class="MSqRob_input_container",
  	    list(
  	tags$label("Annotation file", `for`="annotation", class="MSqRob_label"),
  	tags$button(id="button_annotation", tags$sup("[?]"), class="MSqRob_tooltip"),
		fileInput(inputId="annotation", label=NULL, multiple = FALSE, accept = NULL, width = NULL),
		hidden(helpText(id="tooltip_annotation","Specify the location of your experimental annotation file."))
		)
		)
  	),
		
		#Main panel with number of output and plots
		mainPanel(width = 5,
		          h3("Frequently asked questions", class="MSqRob_topheader"),
		          htmlOutput("folderError"),
		          div(class="MSqRob_h4_container",
		          list(
		          h4("What is an annotation file?"),
		          tags$button(id="button_newExpAnnText",tags$sup("[show]"), class="MSqRob_tooltip")
		          )
		          ),
              hidden(helpText(id="tooltip_newExpAnnText",
                "An experimental annotation file contains the description of your experiment. 
                Indeed, each mass spec run corresponds to e.g. a certain treatment, biological repeat, etc.
                This should be told to MSqRob via an Excel file or a tab delimited file wherein the first column contains all run names 
                and the other columns contain all predictors of interest.
                Examples of experimental annotation files for the Francisella and CPTAC experiments can be found ",
                a("here", href="https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/Francisella/label-free_Francisella_annotation.xlsx"),
                "and",
                a("here.", href="https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/label-free_CPTAC_annotation.xlsx"),
                "Click the button to initialize an Excel file with a \"run\" column (works only if peptides.txt is already uploaded!).
                The annotation file will be saved in the output location.
                You still need to add other relevant columns (treatments, biological repeats, technical repeat, etc.) manually!")),
		          actionButton("create_annot", "Create annotation file", class="MSqRob_button_space"),
		          verbatimTextOutput("newExpAnnText"),
		          div(class="MSqRob_h4_container",
		          list(
		          h4("How do I cite MSqRob?"),
		          tags$button(id="button_cite",tags$sup("[show]"), class="MSqRob_tooltip")
		          )
		          ),
		          hidden(helpText(id="tooltip_cite",
		          "MSqRob is free for you to use and completely open source. 
		          When making use of MSqRob, we would appreciate it if you could cite our two published articles.",
		          br(),
		          span("(1) The MSqRob algorithm: ", class="bold"),
		          br(),
		          "
		          Goeminne, L. J. E., Gevaert, K., and Clement, L. (2016) Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity,
		          and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics 15(2), pp 657-668.",
		          br(),
		          span("(2) The MSqRob GUI tutorial article:", class="bold"),
		          br(),
		          "
		          Goeminne, L. J. E., Gevaert, K. and Clement, L. (2017). 
		          Experimental design and data-analysis in label-free quantitative LC/MS proteomics: 
    		      A tutorial with MSqRob. Journal of Proteomics (in press).")),
		          div(class="MSqRob_h4_container",
		          list(
		          h4("My question is not in this list!"),
		          tags$button(id="button_notinlist",tags$sup("[show]"), class="MSqRob_tooltip")
		          )
		          ),
		          hidden(helpText(id="tooltip_notinlist",
		          "We are always ready to help you with any kind of issue that you might encounter!
		          If for some reason, using MSqRob is hard or counter-intuitive to you,
		          or you encounter some weird and unexpected results,
		          please do not hesitate to contact us at", a("ludger.goeminne@vib-ugent.be.", href="mailto:ludger.goeminne@vib-ugent.be"), 
		          strong("User feedback is very important to us in order to improve MSqRob's user-friendliness.")
		          ))
		)
     )
		# ,
		# #Main panel with number of output and plots
		# mainPanel(width = 5,
		#           htmlOutput("folderError")
		#           )
    )

    ############################
    #Preprocessing tab
    ###########################
    ,tabPanel('Preprocessing',
    sidebarLayout(
	#Sidebar with input
        sidebarPanel(
          
        h3("Settings", class="MSqRob_topheader"),
          
        h4("Transformation", class=c("MSqRob_sidebar")),
        
        div(class="MSqRob_input_container",
            list(
        checkboxInput("logtransform", label="Log-transform data", value=TRUE),
        tags$button(id="button_logtransform", tags$sup("[?]"), class="MSqRob_tooltip"),
        hidden(helpText(id="tooltip_logtransform",
                        "Leave this box ticked to log-transform the data.
                        Log-transformation is almost always performed to make the data less skewed.
                        Only when the data has already been log-transformed, this box can be unticked."
        ))
            )
        ),
        
        div(class="MSqRob_input_container",
            list(
          conditionalPanel(
            condition = "input.logtransform == true",
            tags$label("Base", `for`="log_base", class="MSqRob_label"),
            tags$button(id="button_log_base", tags$sup("[?]"), class="MSqRob_tooltip"),
            numericInput("log_base", NULL, value=2, min = 1, max = NA, step = NA, width = '100%'),
            hidden(helpText(id="tooltip_log_base",
                            "The base of the logarithm when log-transformation is performed.
                            Often, a base of 2 is chosen, because the results will then have interpretations
                            in terms of log2 fold changes.
                            "))
          )
          )
        ),

        h4("Normalization", class="MSqRob_sidebar"),
        
        div(class="MSqRob_input_container",
            list(
        tags$label("Normalization", `for`="normalisation", class="MSqRob_label"),
        tags$button(id="button_normalisation", tags$sup("[?]"), class="MSqRob_tooltip"),
        selectInput("normalisation", NULL, c("quantiles", "quantiles.robust", "vsn", "center.median", "center.mean", "max", "sum", "none"), width = '100%'),
        hidden(helpText(id="tooltip_normalisation",
                        "Select the type of normalization from the dropdown menu.
                        By default, we suggest to use quantile normalization (\"quantiles\").
                        Choose \"none\" if no normalization should be performed 
                        or if the data has already been normalised.
                        "))
            )
        ),

      	h4("Filtering", class="MSqRob_sidebar"),
        
      	#Filter on peptides only modified by site
        div(class="MSqRob_input_container",
            list(
        checkboxInput("onlysite", "Remove only identified by site", value=TRUE),
        tags$button(id="button_onlysite", tags$sup("[?]"), class="MSqRob_tooltip"),
        hidden(helpText(id="tooltip_onlysite",
                        "Keep this box ticked to remove proteins that are exclusively identified
                        by peptides that carry a modification site.
                        This filtering step that is by default performed in a standard Perseus analysis
                        requires you to upload the proteinGroups.txt file.
                        "))
            )
        ),
        
        div(class="MSqRob_input_container",
            list(
      	conditionalPanel(
        	condition = "input.onlysite == true",
        	tags$label("proteinGroups.txt file", `for`="proteingroups", class="MSqRob_label"),
        	tags$button(id="button_proteingroups", tags$sup("[?]"), class="MSqRob_tooltip"),
        	fileInput(inputId="proteingroups", label=NULL, multiple = FALSE, accept = NULL, width = '100%'),
        	hidden(helpText(id="tooltip_proteingroups",
        	             "Specify the location of your proteinGroups.txt file.
        	             This is only needed when you want to remove proteins that are only identified by modified peptides.
        	             "))
      	)
            )
        ),

        div(class="MSqRob_input_container",
            list(
        checkboxInput("smallestUniqueGroups", "Remove comprising protein groups", value=TRUE),
        tags$button(id="button_smallestUniqueGroups", tags$sup("[?]"), class="MSqRob_tooltip"),
        hidden(helpText(id="tooltip_smallestUniqueGroups",
                        "Remove protein groups for which any of its member proteins is present in a smaller protein group.
                        This might be done to remove any overlap of proteins in different protein groups.
                        "))
            )
        ),
        
        #Filter on peptides number of occurances
        div(class="MSqRob_input_container",
            list(
              tags$label("Minimum number of peptides", `for`="minIdentified", class="MSqRob_label"),
              tags$button(id="button_minIdentified", tags$sup("[?]"), class="MSqRob_tooltip"),
      	numericInput("minIdentified", label=NULL, value=2, min = 1, max = NA, step = 1, width = '100%'),
      	hidden(helpText(id="tooltip_minIdentified","
      	                The minimal number of times a peptide sequence should be identified over all samples.
      	                Peptide sequences that are identified less than this number will be removed from the dataset.
      	                The default of 2 has the rationale that it is impossible to discern between the peptide-specific effect and 
      	                any other effects for a peptide that has only been identified once.
      	                "))
      	)
        ),
        
        div(class="MSqRob_input_container",
            list(
              tags$label("Filter columns", `for`="filter", class="MSqRob_label"),
              tags$button(id="button_filter", tags$sup("[?]"), class="MSqRob_tooltip"),
      htmlOutput("selectFilters"),
      hidden(helpText(id="tooltip_filter","
                      Indicate the columns on which filtering should be done.
                      Peptides for which a \"+\" is present in these columns will be removed from the dataset.
                      This kind of filtering is typically done for common contaminants (e.g. operator's keratin)
                      and reversed sequences from the identificiation step that are still present in the data.
                      "))
            )
        )
	),
	#Main panel with number of output and plots
        mainPanel(width = 5,
                  
                  h3("Diagnostic plots", class="MSqRob_topheader"),
                  
                  div(class="MSqRob_input_container MSqRob_h4_checkbox",
                      list(
                        checkboxInput("evalnorm", "Evaluate preprocessing", value=TRUE),
                        tags$button(id="button_evalnorm", tags$sup("[?]"), class="MSqRob_tooltip"),
                        hidden(helpText(id="tooltip_evalnorm","Tick the box to get diagnostic plots for the preprocessing.
                                        Untick the box to skip the making of these diagnostic plots."))
                      )
                  ),
        strong('Number of peptides before preprocessing:'),textOutput('npeptidesRaw',container = span),div(),
        strong('Number of peptides after preprocessing:'),textOutput('npeptidesNormalized',container = span),div(),
        htmlOutput("selectColPlotNorm1"),
        
        
        div(class="MSqRob_h4_container",
            list(
              h4("Intensities after transformation"),
              tags$button(id="button_h4_int_transformation",tags$sup("[?]"), class="MSqRob_tooltip")
            )
        ),
        hidden(helpText(id="tooltip_h4_int_transformation","
                        A density plot showing the distribution of the peptide intensities when only 
                        the transformation is executed. 
                        Transformation is included because a density plot of untransformed intensities is often uninformative 
                        due to a strong skew to the right.
                        Brush and double-click on a selected area to zoom in. 
                        Double click outside a selected area to zoom out.")),
        
        plotOutput('plotRaw',
                   click = "plotRaw_click",
                   dblclick = "plotRaw_dblclick",
                   brush = brushOpts(
                     id = "plotRaw_brush",
                     resetOnNew = TRUE
                   )
                   ),
        
        div(class="MSqRob_h4_container",
            list(
              h4("Intensities after full preprocessing"),
              tags$button(id="button_h4_full_preprocessing",tags$sup("[?]"), class="MSqRob_tooltip")
            )
        ),
        hidden(helpText(id="tooltip_h4_full_preprocessing","
                        A density plot showing the distribution of the peptide intensities 
                        after execution of all preprocessing steps. 
                        This allows you to evaluate the effect of the preprocessing.
                        Brush and double-click on a selected area to zoom in. 
                        Double click outside a selected area to zoom out.")),
        
        plotOutput('plotNorm1',
                   click = "plotNorm1_click",
                   dblclick = "plotNorm1_dblclick",
                   brush = brushOpts(
                     id = "plotNorm1_brush",
                     resetOnNew = TRUE
                   )
                   ),
        
        div(class="MSqRob_h4_container",
            list(
              h4("MDS plot after full preprocessing"),
              tags$button(id="button_h4_MDS_full_preprocessing",tags$sup("[?]"), class="MSqRob_tooltip")
            )
        ),
        hidden(helpText(id="tooltip_h4_MDS_full_preprocessing","A multidimensional scaling plot. This plot shows a two-dimensional scatterplot 
                      so that distances on the plot approximate the typical log2 fold changes between the samples based on a pairwise comparison 
                        of the 500 most different peptides.
                        Brush and double-click on a selected area to zoom in. 
                        Double click outside a selected area to zoom out.")),
        
        
        div(checkboxInput("plotMDSPoints", "Plot MDS points", value=FALSE)),
        div(checkboxInput("plotMDSLabels", "Plot MDS labels", value=TRUE)),
        plotOutput('plotMDS',
				click = "plotMDS_click",
                                 dblclick = "plotMDS_dblclick",
                                 brush = brushOpts(
                                  id = "plotMDS_brush",
                                   resetOnNew = TRUE
                                 )
        )
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
     	  
     	  h3("Settings", class="MSqRob_topheader"),
     	  
  div(class="MSqRob_input_container",
  list(
  tags$label("Group by", `for`="proteins", class="MSqRob_label"),
  tags$button(id="button_proteins", tags$sup("[?]"), class="MSqRob_tooltip"),
	htmlOutput("selectProteins"),
	hidden(helpText(id="tooltip_proteins","
	                Select the level on which the data should be grouped.
	                This is mostly the \"Proteins\" column, as for a traditional shotgun experiment, one is mostly interested in which proteins are differentially abundant.
	                However, sometimes, one would for example like to do inference on the peptides.
	                In these more advanced cases, select the appropriate grouping level.
	                "))
  )
  ),
	
  div(class="MSqRob_input_container",
  list(
  tags$label("Annotation columns", `for`="annotations", class="MSqRob_label"),
  tags$button(id="button_annotations", tags$sup("[?]"), class="MSqRob_tooltip"),
	htmlOutput("selectAnnotations"),
	hidden(helpText(id="tooltip_annotations","
	                Some input files contain additional information about a protein, such as a gene name, a full protein name, a pathway annotation, etc.
	                If you whish to see these columns in the output, select these additional annotation	columns here.
	                "))
  )
  ),
	
  div(class="MSqRob_input_container",
  list(
  tags$label("Fixed effects", `for`="fixed", class="MSqRob_label"),
  tags$button(id="button_fixed", tags$sup("[?]"), class="MSqRob_tooltip"),
	htmlOutput("selectFixed"),
	hidden(helpText(id="tooltip_fixed","
	                Select the fixed effects.
	                Fixed effects are effects for which all levels of interest are included in the experiment.
	                They are typically controlled for by the experimenter and the levels would have exactly the same interpretation
	                if the experiment were to be repeated. Examples include all effects of interest and effects with a limited number of levels, such as \"gender\".
	                When in doubt which effects should be included as fixed, please contact us!"))
  )
  ),
  
  div(class="MSqRob_input_container",
  list(
  tags$label("Random effects", `for`="random", class="MSqRob_label"),
  tags$button(id="button_random", tags$sup("[?]"), class="MSqRob_tooltip"),
	htmlOutput("selectRandom"),
  hidden(helpText(id="tooltip_random","
                  Select the random effects.
	                Random effects are effects that are random draws from a (theoretically) infinite population that are generally not of interest, 
	                but are only included to build a correct covariance structure in the model.
	                They are not controlled for by the experimenter and its levels have different interpretations each time the experiment would be repeated
	                (e.g. \"mouse 1\" from the first experiment is not the same \"mouse 1\" in the repeated experiment).
	                Effects of biological or technical replicates are always random effects.
	                As peptide sequence effects are often very strong, they might overwhelm the remaining effects in the experiment.
	                Therefore, we suggest specifying the peptide sequence effect as a separate random effect to allow the remaining fixed effects of interested to be penalized independently of the peptide effect.
                  When in doubt which effects should be included as random, please contact us!"))
  )
  ),
  
	#checkboxInput("borrowFixed", "Borrow information across fixed effects", value = FALSE, width = NULL),
	#checkboxInput("borrowRandom", "Borrow information across random effects", value = FALSE, width = NULL),
  div(class="MSqRob_input_container",
      list(
        tags$label("Save/load options:", `for`="save", class="MSqRob_label"),
        tags$button(id="button_save", tags$sup("[?]"), class="MSqRob_tooltip"),
  radioButtons("save", label=NULL,
                   c("Save the models" = 1,  #Conditioneel stukje toevoegen bij downloadHandler dat hij ook model.RData downloadt
                     "Load existing models" = 2,
                     "Don't save the models" = 3
                      )),
  hidden(helpText(id="tooltip_save","
                  Select one of three options:
                  (1) Save the models: the data and the fitted models will be saved as a \".rDatas\" file which allows the user to redo the analysis without having to redo the preprocessing and model fitting.
                  (2) Load existing models: allows the user to upload a \".rDatas\" file from a previous analysis. When this option is selected, all input except the type of analysis and the contrast options will be disabled because the model has already been fitted to the data.
                  (3) Don't save the models: No \".rDatas\" file will be saved (but a results.xlsx file will be produced nonetheless).
                  "))
      )
  ),
  
	#load saved model
  div(class="MSqRob_input_container",
      list(
        conditionalPanel(
        	condition = "input.save == 2",
        	tags$label("Model file", `for`="load_model", class="MSqRob_label"),
        	tags$button(id="button_load_model", tags$sup("[?]"), class="MSqRob_tooltip"),
        	fileInput(inputId="load_model", label=NULL, multiple = FALSE, accept = NULL, width = NULL),
        	hidden(helpText(id="tooltip_load_model","Specify the location of your saved model file."))
      	)
        )
  ),
  
	#Type of analysis
     	#selectInput("analysis_type", "Select the type of analysis", c("standard", "stagewise", "ANOVA")),
  div(class="MSqRob_input_container",
      list(
        tags$label("Analysis type", `for`="analysis_type", class="MSqRob_label"),
        tags$button(id="button_analysis_type", tags$sup("[?]"), class="MSqRob_tooltip"),
  selectInput("analysis_type", label=NULL, c("standard", "ANOVA")),
  hidden(helpText(id="tooltip_analysis_type","
                  In a standard analysis, each contrast will be estimated separately (default).
                  In an ANOVA analysis, one assesses whether there is on average an effect over the specified contrasts.
                  It is often more of an exploratory tool to identify proteins that are significant in many contrasts at once.
                  Note that the average expression for an ANOVA can be close to er even equal to zero for very significant effects 
                  (i.e. when some contrasts are highly positive and others highly negative).
                  "))
      )
  ),

  #Number of contrasts
  div(class="MSqRob_input_container",
      list(
        tags$label("Number of contrasts", `for`="nContr", class="MSqRob_label"),
        tags$button(id="button_nContr", tags$sup("[?]"), class="MSqRob_tooltip"),
      	numericInput("nContr", label=NULL, value=1, min = 1, max = NA, step = 1, width = NULL),
        hidden(helpText(id="tooltip_nContr","
                        The number of research hypotheses you want to test.
                        Each research hypothesis corresponds to a so-called statistical contrast.",
                        br(),
                        br(),
                        span("How to specify a contrast: ", class="bold"),
                        br(),
                        br(),
                        "If you followed the default preprocessing, data is log2 transformed. Therefore, a simple log2 fold change between two conditions 
                        can be tested by specifying \"1\" for the first condition and \"-1\" for the second condition
                        (a difference of log2 transformed values is equal to a log2 of their quotient).
                        An average over e.g. 3 conditions can be tested by specifying \"1/3\" for each of these conditions.
                        As a general rule of thumb, the sum of each contrast should amount to \"0\" or \"1\". If not, you are probably doing something wrong.",
                        strong("Looking at the right contrast is crucial to get the right output!"),
                        "Therefore, if you need help in specifying your contrasts, please do not hesitate to contact us!
                        "))
      )
  ),
  
        #Specification of contrasts
      	htmlOutput("selectLevels"),
	#Run button
	actionButton(inputId="go", label="Start the analysis!", class="MSqRob_button_space") #,
	#Download button
      	# htmlOutput("download_button")
    ),

    #Main panel with results and plots
	mainPanel(
	  
	  h3("Results", class="MSqRob_topheader"),
	  
	verbatimTextOutput("nText"),
	htmlOutput("plot_contrast"),
	verbatimTextOutput("contrastL"),

	#Volcano plot
       	fluidRow(
        column(width = 6,  #6 out of 12 => half the screen!

        	div(class="MSqRob_h4_container",
        	    list(
        	      h4("Volcano plot"),
        	      tags$button(id="button_h4_volcano_plot",tags$sup("[?]"), class="MSqRob_tooltip")
        	    )
        	),
        	hidden(helpText(id="tooltip_h4_volcano_plot","
        	                The volcano plot plots the estimated fold change or average expression as a function of log10 of the p-value.
        	                Points shown in red are significant at the specified FDR level, points shown in black are not.
        	                Select and deselect points by clicking on them either in the volcano plot or in the results table.
        	                Brush and double-click on the selected area to zoom in. Double click outside the selected area to zoom out.
        	                Adjust the significance level to visualize features with an FDR level below alpha.")),
        	
              	plotOutput("plot1", height = 300,
                           click = "plot1_click",
                           dblclick = "plot1_dblclick",
                           brush = brushOpts(
                             id = "plot1_brush",
                             resetOnNew = TRUE
                           )
                ),
		actionButton("add_area_selection", "Add selected area to selection", class="MSqRob_button_space"),
    actionButton("remove_area_selection", "Remove selected area from selection", class="MSqRob_button_space"),
		actionButton("remove_all_selection", "Remove everything from selection", class="MSqRob_button_space"),
		
		#Input significance level (alpha)
		div(class="MSqRob_input_container",
		    list(
		      tags$label("Significance level", `for`="alpha", class="MSqRob_label"),
		      tags$button(id="button_alpha", tags$sup("[?]"), class="MSqRob_tooltip"),
		      numericInput("alpha", label=NULL, value=.05, min = 0, max = 1, step = 0.01, width = NULL),
		      hidden(helpText(id="tooltip_alpha","Select the significance level (alpha) at which the type I error needs to be performed.
		                      Tests are traditionally performed at the 5% false discovery rate (FDR) level, but more stringent control (e.g. 1% FDR or even less) is sometimes adopted in experiments where false positives are highly unwanted (e.g. clinical settings).
		                      The lower this level, the more stringent the cut-off and thus the less proteins that will be declared significant, but the higher the statistical certainty of the differential abundance of these proteins.
		                      An FDR of 5% means that on average an expected 5% of the proteins that are called significant will be in fact false positives."))
		    )
		)
	),
        #Detail plot
	column(width = 6, #6 out of 12 => half the screen!

	       
	       div(class="MSqRob_h4_container",
	           list(
	             h4("Detail plot"),
	             tags$button(id="button_h4_detail_plot",tags$sup("[?]"), class="MSqRob_tooltip")
	           )
	       ),
	       hidden(helpText(id="tooltip_h4_detail_plot","
	                       If only one data point is selected in the results table, 
	                       this plot shows the individual log2 peptide intensities.
	                       This plot is extremely useful to visually assess the evidence for differential abundance in the data.")),
	       
	       
                plotOutput("plot2", height = 300,
                           click = "plot2_click",
                           dblclick = "plot2_dblclick",
                           brush = brushOpts(
                             id = "plot2_brush",
                             resetOnNew = TRUE
                           )
                ),
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


