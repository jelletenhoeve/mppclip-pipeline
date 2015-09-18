###
### This script trains the model and makes a prediction 
###



genesets <- loadGenesets()
pathways <- loadPathways()


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {

	## load the data
	mEset <- getDataSet(sample_group, c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'))

	# do global tests for both genesets and pathways
	mirGlobalTests <- globalTest(mEset, genericSubsets=c(genesets, pathways))

	message(sample_group, ': writing files')
	writeObjectOutput(mirGlobalTests, sample_group, name='mir-based_globaltests', pipeline_step='genesets-and-pathways-associations')

	return(paste(sample_group, 'done!'))
})

