library(Biobase)
library(gdata)


source('../../lib/AllClasses.R')
source('../../lib/AllGenerics.R')
source('../../lib/MatchedExpressionSet.R')
source('../../lib/helpers.R')

source('../../w_pipeline/global/data_functions.R')

load('config.rda')

sample_group_choices <- run_config[['sample_groups']]
names(sample_group_choices) <- paste(run_config[['dataset_name']], '-', run_config[['sample_groups']])



#getSubSeries <- function(serie) {
#	if (serie == "") {
#		"[select serie above]"
#	} else {
#		dataset_files <- dir('../w_data/datasets/')
#		
#		series <- sapply(strsplit(dataset_files, "\\."), function(x) x[[1]])
#
#		sub_series <- sapply(strsplit(dataset_files, ".samples-"), function(x) x[[2]])
#		sub_series <- sapply(sub_series, function(x) substr(x, 1, nchar(x)-4))
#
#		names(sub_series) <- sub_series
#
#		c("", sub_series[series == serie])
#	}
#}

#getNomenclatures <- fucntion(serie, sub_serie) {
#	if (serie == "") {
#		"[select sub serie above]"
#	} else {
#		dataset_files <- dir('../w_data/datasets/')
#		
#		series <- sapply(strsplit(dataset_files, "\\."), function(x) x[[1]])
#
#		sub_series <- sapply(strsplit(dataset_files, ".samples-"), function(x) x[[2]])
#		sub_series <- sapply(sub_series, function(x) substr(x, 1, nchar(x)-4))
#		names(sub_series) <- sub_series
#
#
#		sub_series <- sapply(strsplit(dataset_files, ".mirbase-"), function(x) x[[2]])
#
#		c("", sub_series[series == serie])
#	}
#}

