
###
### This does the rank tests
###




mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	## load the sample data
	mEset <- getDataSet(sample_group, c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'), showWarnings=FALSE)

    interaction_features <- readObjectOutput(sample_group = sample_group, name='interaction_features', pipeline_step='interaction-features', for_shiny=FALSE)    
	prediction.df <- readObjectOutput(sample_group, name='prediction', pipeline_step='model-prediction')


	if (!all(paste(prediction.df[, c('gene', 'mirna')], sep=':') == paste(interaction_features[, c('gene', 'mirna')], sep=':')))
		stop('ordering of "interaction_features" and "prediction.df" is different!')


	full_table <- cbind(prediction.df, interaction_features[, c('Correlation', 'MCF7_mRNA_RPKM')])

	
	message(sample_group, ': writing tables...')
	writeTableOutput(full_table, sample_group, name='ts_interaction_features', pipeline_step='tables-and-figures')



	return(paste(sample_group, 'done!'))
})
