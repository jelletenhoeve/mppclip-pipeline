library(shiny)

# Define UI for dataset viewer application



ranking_table_choices <- names(run_config[['ranking_tables']])
names(ranking_table_choices) <- sapply(run_config[['ranking_tables']], function(x) paste(names(x), collapse=', '))


shinyUI(pageWithSidebar(

	# Application title
	headerPanel("MicroRNA target explorer for breast cancer data"),

	# Sidebar with controls to select a dataset and specify the number
	# of observations to view
	sidebarPanel(
		


		selectInput("sample_group", "Choose a study / sample group:", choices = sample_group_choices),
		tags$hr(),
		tags$h4("Data information"),
		tags$h5("Sizes"),
		tags$span("- samples: "),
		textOutput("n_samples", inline=TRUE),
		tags$br(),
		tags$span("- genes: "),
		textOutput("n_genes", inline=TRUE),
		tags$br(),
		tags$span("- mirnas: "),
		textOutput("n_mirnas", inline=TRUE),
		tags$br(),
		
		tags$h5("Nomenclatures"),
		tags$span("- mirna: "),
		tags$span(run_config[['mir_nomenclature']]),
		tags$br(),
		tags$span("- gene: HUGO symbols"),
		
		tags$h5("Thresholds"),
		tags$span("- mirna: "),
		tags$span(run_config[['mir_threshold']]),
		tags$br(),
		tags$span("- gene: "),
		tags$span(run_config[['gene_threshold']]),
		
		tags$h5("Clinical characteristics"),
		tags$span(paste(names(run_config[['pheno_columns']]), collapse=', ') ),
		
		tags$h5("Clinical characteristics used in null model"),
		tags$span(paste(names(run_config[['null_pheno_columns']]), collapse=', ') ),

		tags$hr(),

		tags$h4("Ranking settings"),
		helpText("These settings only affect the MicroRNA ranking tables (marked with *)"),
		selectInput("target_set", "Target set", choices = global_config[['target_sets']]),
		selectInput("ranking_table", "Features set", choices = ranking_table_choices),

		tags$hr(),

		tags$h4("Model setting"),
		tags$h5("GLM with lasso or elasticnet regularization (glmnet)"),
		tags$span("- alpha: "),
		tags$span(run_config[['model_training_config']][['alpha']]),
		tags$br(),
		tags$span("- family: "),
		tags$span(run_config[['model_training_config']][['family']]),
		tags$br(),
		tags$span("- measure type: "),
		tags$span(run_config[['model_training_config']][['type.measure']]),
		tags$br(),
		tags$span("- lambda (penalty parameter): "),
		tags$span(run_config[['model_training_config']][['s']]),
		tags$br(),
		tags$span("- cutoff used for prediction: "),
		tags$span(run_config[['model_training_config']][['cutoff']]),

		tags$h5("Features used for training"),
		tags$span(paste(run_config[['model_training_config']][['features']], collapse=', ') ),
		tags$br(),

		tags$h5("Training set"),
		tags$span(run_config[['model_training_config']][['filter_variable']]),
		tags$h5("Response variable"),
		tags$span(run_config[['model_training_config']][['response_variable']]),


		tags$hr(),

		tags$h4("Download"),
		downloadButton('download_code_book', 'Code book'),
		downloadButton('download_figures_and_tables', 'Figures and tables'),

		tags$hr(),

		tags$h4("Manuscript"),
		tags$strong("Identification of distinct miRNA target regulation between breast cancer molecular subtypes using AGO2-PAR-CLIP and patient datasets"),
		tags$br(),
		tags$br(),
		tags$em("Thalia A. Farazi*, Jelle J. ten Hoeve*, Miguel Brown, Aleksandra Mihailovic, Hugo M. Horlings, Marc J. van de Vijver, Thomas Tuschl# and Lodewyk F.A. Wessels#"),
		tags$hr()
	),

#	mainPanel(
#		tabsetPanel(
#	  		tabPanel("Model prediction", dataTableOutput("model_prediction")),
#  			tabPanel("Interaction features", dataTableOutput("interaction_features")),
#  			tabPanel("Pairs overlap", imageOutput("pairs_overlap"))
#		)
#	)
	mainPanel(
		tabsetPanel(
	  		tabPanel("Threshold analysis",
	  			tabsetPanel(
		  			tabPanel("MicroRNA abundance thresholds", imageOutput("mirna_abundance_thresholds")),
		  			tabPanel("MicroRNA variance thresholds", imageOutput("mirna_variance_thresholds")),
		  			tabPanel("Gene abundance thresholds", imageOutput("gene_abundance_thresholds")),
		  			tabPanel("Gene variance thresholds", imageOutput("gene_variance_thresholds"))
	  			)
	  		),
	  		tabPanel("Model",
	  			tabsetPanel(
		  			tabPanel("Prediction and reponse", dataTableOutput("prediction")),
		  			tabPanel("Venn diagram", imageOutput("pairs_overlap")),
		  			tabPanel("Coefficients", dataTableOutput("coefficients")),
		  			tabPanel("Performance", imageOutput("glmnet_prediction_measures"))
	  			)
	  		),
	  		tabPanel("Genesets and pathways",
	  			tabsetPanel(
		  			tabPanel("Associations", dataTableOutput("mir_based_globaltests")),
		  			tabPanel("Pathway clustering", imageOutput("pathway_clustering")),
		  			tabPanel("Pathway characteristics", dataTableOutput("pathway_characteristics"))
	  			)
	  		),
	  		tabPanel("Survival and other clinical characteristics",
	  			tabsetPanel(
		  			tabPanel("Associations", dataTableOutput("phenotype_globaltests")),
			  		tabPanel("MicroRNA Ranking", dataTableOutput("mirna_ranking"))
	  			)
	  		),	  		
	  		tabPanel("Overall MicroRNA ranking",
	  			tabsetPanel(
		  			tabPanel("MicroRNA features", dataTableOutput("mirna_features")),
		  			tabPanel("PCA analysis", imageOutput("pca_biplot")),
		  			tabPanel("MicroRNA Ranking", dataTableOutput("mirna_ranking2"))
	  			)
	  		)
	  		,
	  		tabPanel("Targets",
	  			tabsetPanel(
		  			tabPanel("MP-PCLIP", dataTableOutput('mppclip_targets')),
		  			tabPanel("TS interaction features", dataTableOutput('ts_interaction_features')),
		  			tabPanel("Correlation analysis", imageOutput('correlation_analysis'))
	  			)
	  		)
		)
	)
))