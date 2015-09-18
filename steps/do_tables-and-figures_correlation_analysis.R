
###
### This does the rank tests
###


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	## load the sample data
	mEset <- getDataSet(sample_group, c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'), showWarnings=FALSE)
	

	startPngOutput(sample_group=sample_group, name='correlation_analysis', pipeline_step='tables-and-figures', width=640, height=640)
	plotEcdfs(mEset, interaction.names=c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'), alternative='less', xlim=c(-1, 1))
	dev.off()

})
