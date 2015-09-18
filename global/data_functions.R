## helpers

getDataSet <- function(sample_group=NULL, interactions_sets=c(), showWarnings = TRUE, for_shiny=FALSE) {
	attach(run_config)


	if (is.null(sample_group)) {
		# get the overall data set
		if (for_shiny) {
			dataset_dir <- paste(
				'../../w_data/datasets/',
				dataset_name,
				sep=''
			)
		} else {
			dataset_dir <- paste(
				'../w_data/datasets/',
				dataset_name,
				sep=''
			)
		}

		dataset_file <- paste(
			dataset_name, '-',
			mir_nomenclature, '.rda', 
			sep=''
		)

	} else {
		# get the sample group specific dataset
		if (for_shiny) {
			dataset_dir <- paste(
				'../../w_data/datasets/',
				dataset_name, '.',
				mir_nomenclature, '-gt', mir_threshold, '.',
				'gene-gt', gene_threshold,
				sep=''
			)
		} else {
			dataset_dir <- paste(
				'../w_data/datasets/',
				dataset_name, '.',
				mir_nomenclature, '-gt', mir_threshold, '.',
				'gene-gt', gene_threshold,
				sep=''
			)
		}

		dataset_file <- paste(
			dataset_name, '.',
			mir_nomenclature, '-gt', mir_threshold, '.',
			'gene-gt', gene_threshold, '.',
			'samples-', sample_group, '.rda',
			sep=''
		)
	}

	load(paste(dataset_dir, '/', dataset_file, sep=''))

	for (interactions_set in interactions_sets) {

		if (interactions_set == 'model_prediction') {

		    prediction.df <- readTableOutput(sample_group = sample_group, name='prediction', pipeline_step='model-prediction', for_shiny=for_shiny)			
			preds <- prediction.df[prediction.df$prediction, c('mirna', 'gene')]
			preds[[1]] <- as.character(preds[[1]])
			preds[[2]] <- as.character(preds[[2]])

			interactions(mEset, 'model_prediction', showWarnings = showWarnings) <- preds

		} else if (interactions_set == 'model_prediction_parclip_union') {

		    prediction.df <- readTableOutput(sample_group = sample_group, name='prediction', pipeline_step='model-prediction', for_shiny=for_shiny)			
			preds <- prediction.df[prediction.df$prediction, c('mirna', 'gene')]
			preds[[1]] <- as.character(preds[[1]])
			preds[[2]] <- as.character(preds[[2]])

			loadExtData('parclip', mir_nomenclature=mir_nomenclature)

			interactions(mEset, 'model_prediction_parclip_union', showWarnings = showWarnings) <- unique(rbind(
				preds, 
				get('parclip')
			))

		} else {
			loadExtData(interactions_set, mir_nomenclature=mir_nomenclature)
			interactions(mEset, interactions_set, showWarnings = showWarnings) <- get(interactions_set)
		}
	}

	detach(run_config)

	mEset
}


existsObjectOutput <- function(sample_group=NULL, name, pipeline_step, for_shiny = FALSE) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=for_shiny)
	rda_file <- paste(output_dir, '/', name, '.rda', sep='')

	file.exists(rda_file)
}


readObjectOutput <- function(sample_group=NULL, name, pipeline_step, for_shiny = FALSE) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=for_shiny)
	rda_file <- paste(output_dir, '/', name, '.rda', sep='')

	load(rda_file)

	if (for_shiny & is.data.frame(x)) {
		# read the mapping table
		mapping <- read.table(file='column_names.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
		# filter
		mapping <- mapping[mapping$pipeline_step == pipeline_step & startsWith(rep(name, nrow(mapping)), mapping$name), ]

		if (nrow(mapping) == 0) {
			stop('No mapping found for ', pipeline_step, ', ', name)
		}


		idx <- match(colnames(x), mapping$column)
		if (any(is.na(idx))) {
			stop('pretty column names ', paste(colnames(x)[is.na(idx)], collapse=', '), ' not found in "column_names.txt" for ', pipeline_step, '-', name)
		}

		colnames(x) <- mapping$pretty_name[idx]
		x <- x[, order(mapping$order[idx])]


		datatable(
			x, filter = 'bottom', options = list(pageLength = 25)
		)
	} else {
		x
	}
	
}

readTableOutput <- function(sample_group=NULL, name, pipeline_step, for_shiny = FALSE) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=for_shiny)
	table_file <- paste(output_dir, '/', name, '.csv', sep='')

	df <- read.table(file=table_file, header=TRUE, sep=',', row.names=1)

	if (for_shiny) {

		# read the mapping table
		mapping <- read.table(file='column_names.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
		# filter
		mapping <- mapping[mapping$pipeline_step == pipeline_step & startsWith(rep(name, nrow(mapping)), mapping$name), ]

		if (nrow(mapping) == 0) {
			stop('No mapping found for ', pipeline_step, ', ', name)
		}


		idx <- match(colnames(df), mapping$column)
		if (any(is.na(idx))) {
			stop('Not all columns were found in "column_names.txt" for ', pipeline_step, '-', name)
		}

		colnames(df) <- mapping$pretty_name[idx]
		df <- df[, order(mapping$order[idx])]

		datatable(
			df, filter = 'bottom', options = list(pageLength = 25)
		)
	} else {
		df		
	}
}


writeTableOutput <- function(x, sample_group=NULL, name, pipeline_step) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=FALSE, create_if_not_exists=TRUE)
	table_file <- paste(output_dir, '/', name, '.csv', sep='')

	write.csv(x, file=table_file, row.names=TRUE)
}

writeObjectOutput <- function(x, sample_group=NULL, name, pipeline_step) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=FALSE, create_if_not_exists=TRUE)
	rda_file <- paste(output_dir, '/', name, '.rda', sep='')

	save(x, file=rda_file)
}

startPdfOutput <- function(sample_group=NULL, name, pipeline_step, ...) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=FALSE, create_if_not_exists=TRUE)
	pdf_file <- paste(output_dir, '/', name, '.pdf', sep='')

	pdf(file=pdf_file, ...)
}

startPngOutput <- function(sample_group=NULL, name, pipeline_step, ...) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=FALSE, create_if_not_exists=TRUE)
	png_file <- paste(output_dir, '/', name, '.png', sep='')

	png(file=png_file, ...)
}

getImageSrc <- function(sample_group=NULL, name, pipeline_step, for_shiny = FALSE) {
	output_dir <- getOutputDir(sample_group=sample_group, pipeline_step=pipeline_step, for_shiny=for_shiny)
	png_file <- paste(output_dir, '/', name, '.png', sep='')
	
	list(src = png_file, contentType = 'image/png')
}


getOutputDir <- function(sample_group=NULL, pipeline_step, for_shiny=FALSE, create_if_not_exists=FALSE) {

	pipeline_step <- match.arg(pipeline_step, choices=global_config[['pipeline_steps']])

	nomthres <- paste(run_config[['mir_nomenclature']], '-gt', run_config[['mir_threshold']], '.', 'gene-gt', run_config[['gene_threshold']], sep ='')
	
	if (for_shiny) {
		if (is.null(sample_group)) {
			output_dir <- paste('data/', run_config[['dataset_name']], '/_all/', pipeline_step, '', sep='')
		} else {
			output_dir <- paste('data/', run_config[['dataset_name']], '/', nomthres, '/', sample_group, '/', pipeline_step, '', sep='')
		}
	} else {
		if (is.null(sample_group)) {
			output_dir <- paste(run_config[['output_dir']], '/data/', run_config[['dataset_name']], '/_all/', pipeline_step, '', sep='')
		} else {
			output_dir <- paste(run_config[['output_dir']], '/data/', run_config[['dataset_name']], '/', nomthres, '/', sample_group, '/', pipeline_step, '', sep='')
		}
	}

	if (create_if_not_exists) {
		dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
	}

	output_dir
}


loadGenesets <- function() {
	### load all genesets
	genesets_dir <- "data/genesets"
	genesets <- sapply(run_config[['genesets']], function(geneset) {
		read.delim(
			file=paste(genesets_dir, '/', geneset, '.txt', sep=''),
			header=FALSE,
			stringsAsFactors = FALSE
		)$V1
	})
	names(genesets) <- run_config[['genesets']]
	genesets
}

loadPathways <- function() {
	if (run_config[['pathways']] == 'kegg_pathways') {
		library("org.Hs.eg.db")
		library(KEGG.db)

		eg2symbol <- as.list(org.Hs.egSYMBOL)

		kegg_pathways        <- as.list(KEGGPATHID2EXTID)                                       # get all KEGG pathways
		kegg_pathways        <- kegg_pathways[grep('hsa', names(kegg_pathways))]                # filter for 'hsa'
		kegg_pathways        <- lapply(kegg_pathways, function(x) unlist(eg2symbol[x], use.names = FALSE))         # map to GENE symbol
		names(kegg_pathways) <- paste(names(kegg_pathways), as.list(KEGGPATHID2NAME)[substring(names(kegg_pathways), 4)], sep=' - ')
		rm(eg2symbol)
		kegg_pathways
	} else {
		stop('No loading functionality for ', run_config[['pathways']])
	}
}

loadIntronicMirs <- function() {
	# todo
}


loadExtData <- function(ext_dataset_name, mir_nomenclature=NULL, as.list.if.possible=FALSE) {


	#
	# TargetScan data
	#

	if (ext_dataset_name == 'targetscan' & mir_nomenclature == 'targetscan_family') {

		load('data/targetscan/sc_hsa.rda')
		sc_hsa$Transcript.ID <- NULL
		sc_hsa$Species.ID <- NULL

		colnames(sc_hsa)[1] <- 'gene'
		colnames(sc_hsa)[2] <- 'mirna'

		assign("targetscan", sc_hsa, envir = .GlobalEnv)
	}
	
	else if (ext_dataset_name == 'targetscan' & mir_nomenclature == 'mirbase') {
		
		load('data/targetscan/sc_hsa.rda')
		# swith columns miRNA.family and Representative.miRNA (this is the mirbase nomenclature)
		sc_hsa <- sc_hsa[, c(1,2,13,4,5,6,7,8,9,10,11,12,3,14,15)]

		sc_hsa$Transcript.ID <- NULL
		sc_hsa$Species.ID <- NULL


		colnames(sc_hsa)[1] <- 'gene'
		colnames(sc_hsa)[2] <- 'mirna'

		assign("targetscan", sc_hsa, envir = .GlobalEnv)
	}
	
	else if (ext_dataset_name == 'targetscan_predicted_conserved' & mir_nomenclature == 'targetscan_family') {
		targetscan_predicted_conserved <- read.delim(file="data/targetscan/targetscan_predicted_conserved.txt", stringsAsFactors=FALSE)
		assign("targetscan_predicted_conserved", targetscan_predicted_conserved, envir = .GlobalEnv)
	}


	else if (ext_dataset_name == 'targetscan_predicted_conserved' & mir_nomenclature == 'mirbase') {
		load('data/targetscan/mirbase_targetscan_predicted_conserved.rda')
		
		colnames(mirbase_targetscan_predicted_conserved)[1] <- 'mirna'
		colnames(mirbase_targetscan_predicted_conserved)[2] <- 'gene'

		assign("targetscan_predicted_conserved", mirbase_targetscan_predicted_conserved, envir = .GlobalEnv)
	}



	#
	# MCF7 data
	#
	else if (ext_dataset_name == 'mcf7') {
		data(mcf7)
		colnames(mcf7)[1] <- 'gene'

		assign("mcf7", mcf7, envir = .GlobalEnv)
	}

	


	#
	# PARCLIP data
	#

	else if (ext_dataset_name == 'parclip' & mir_nomenclature == 'targetscan_family') {
		data(parclip)
		colnames(parclip) <- c('mirna', 'gene')
		if (as.list.if.possible) {
			parclip.list <- split(parclip$gene, parclip$mirna)
			assign("parclip", parclip.list, envir = .GlobalEnv)
		} else {
			assign("parclip", parclip, envir = .GlobalEnv)
		}
	}
	
	else if (ext_dataset_name == 'parclip' & mir_nomenclature == 'mirbase') {
		data(parclip_mirbase_all)
		colnames(parclip_mirbase_all) <- c('mirna', 'gene')
		if (as.list.if.possible) {
			parclip_mirbase_all.list <- split(parclip_mirbase_all$gene, parclip_mirbase_all$mirna)
			assign("parclip", parclip_mirbase_all.list, envir = .GlobalEnv)
		} else {
			assign("parclip", parclip_mirbase_all, envir = .GlobalEnv)
		}
	}

	#
	# Throw an error, because requested data is not available..
	#

	else {
		stop('No external data available for "ext_dataset_name"=', ext_dataset_name, ' and "mir_nomenclature"=', mir_nomenclature) 
	}
}
