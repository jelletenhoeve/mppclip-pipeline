#run.R


# load configs
#source('run_configs/metabric_mirbase_92_60.R')
#source('run_configs/metabric_tsfam_92_60.R')
source('run_configs/tcga_tsfam_1e-04_0.R')


source('global/global_config.R')



# load libraries
library(Biobase)
library(globaltest)
library(reshape)
library(ROCR)
library(glmnet)
library(parallel)
library(survival)
library(gdata)

source('../lib/AllClasses.R')
source('../lib/AllGenerics.R')
source('../lib/MatchedExpressionSet.R')
source('../lib/helpers.R')
source('../w_pipeline/global/data_functions.R')





# 
# Run each step of the pipeline
# 
ordered_step_scripts <- c(
	step0 = 'do_preparation.R',

	step1 = 'do_threshold-analysis.R',

	step2 = 'do_interaction-features.R',

	step3 = 'do_model-prediction.R',
	step4 = 'do_model-prediction_venn-diagram.R',
	step5 = 'do_model-prediction_coefficients-and-performance.R',

	step6 = 'do_genesets-and-pathways-associations.R',
	step7 = 'do_genesets-and-pathways-associations_associations-table-and-heatmap.R',
	
	step8 = 'do_phenotype-associations.R',

	step9 = 'do_mirna-ranking.R',

	step10 = 'do_tables-and-figures_pathway-characteristics-table.R', 
	step11 = 'do_tables-and-figures_mppclip-table.R',
	step12 = 'do_tables-and-figures_ts-interaction-features-table.R',
	step13 = 'do_tables-and-figures_correlation_analysis.R',
	
	stepFinal = 'do_finalization.R'
)

ordered_step_scripts <- ordered_step_scripts[c('step0')]

for (step in names(ordered_step_scripts)) {

	if (step != 'step0') {
		zz <- file(paste(run_config[['output_dir']], "/log/pipeline-", step, ".txt", sep=''), open = "wt")
		sink(zz)
		sink(zz, type = "message")
	}
	script <- ordered_step_scripts[step]

	message('Start: ', step, ' = ', script)
	source(paste('steps/', script, sep=''))
	message('End: ', step, ' = ', script)

	if (step != 'step0') {
		sink(type = "message")
		sink()
	}
}




