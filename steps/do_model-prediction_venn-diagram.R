###
### This script generates a Venn Diagram allowing to inspect the overlap of
### predicted, data, parclip and target scan data
###

model_training_config <- run_config[['model_training_config']]



loadExtData('targetscan_predicted_conserved', mir_nomenclature=run_config[['mir_nomenclature']])
loadExtData('parclip', mir_nomenclature=run_config[['mir_nomenclature']])
loadExtData('targetscan', mir_nomenclature=run_config[['mir_nomenclature']])

sc  <- unique(targetscan[, c('mirna', 'gene')])

tsi_pairs <- paste(targetscan_predicted_conserved$mirna, targetscan_predicted_conserved$gene, sep=':')
pc_pairs  <- paste(parclip[[1]], parclip[[2]], sep=':')
sc_pairs  <- paste(sc$mirna, sc$gene, sep=':')

rm(sc)




#
# Start the loop
#

mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	set.seed(run_config[['seed']])


	## load the data
	mEset <- getDataSet(sample_group)

	message('Loading prediction from table file')
    prediction.df <- readObjectOutput(sample_group = sample_group, name='prediction', pipeline_step='model-prediction', for_shiny=FALSE)    




	###
	### Do Venn diagrams
	###

	message('Writing Venn diagrams...')


	# Pairs Venn
	mirs  <- featureNames(eset(mEset, id = 'mirna'))
	genes <- featureNames(eset(mEset, id = 'mrna'))

	data_pairs <- paste(
		rep(mirs, each=length(genes)),
		rep(genes, length(mirs)),
		sep=':'
	)



	
	predicted.idx <- prediction.df$prediction
	predicted_pairs <- paste(prediction.df[[1]][predicted.idx], prediction.df[[2]][predicted.idx], sep=':')
	


	pairs_list <- list(
		parclip       = pc_pairs,
		ts_pred_cons  = tsi_pairs,
		data          = data_pairs,
		ts_sum_counts = sc_pairs,
		prediction    = predicted_pairs
	)


#	idx <- predicted_pairs %in% pairs_list[['parclip']] &
#		predicted_pairs %in% pairs_list[['ts_pred_cons']] &
#		predicted_pairs %in% pairs_list[['prediction']]
#
#	write.table(predicted_pairs[idx], file=interactions_in_all_file, quote=FALSE, row.names=FALSE, col.names=FALSE)

	names(pairs_list)[3] <- name(mEset)
	names(pairs_list)[5] <- paste('prediction', model_training_config[['cutoff']], sep='-')
	sizes <- sapply(pairs_list, length)

	names(pairs_list) <- paste(names(pairs_list), ' (', sizes, ')', sep='')


	startPdfOutput(sample_group, name='pairs_overlap', pipeline_step='model-prediction', width=10, height=10)
	venn(pairs_list[c(1, 2, 4, 5, 3)])
	dev.off()

	startPngOutput(sample_group, name='pairs_overlap', pipeline_step='model-prediction', width=960, height=960)
	old_par <- par(mar=c(1,1,1,1), cex=.5)
	venn(pairs_list[c(1, 2, 4, 5, 3)])
	par(old_par)
	dev.off()




	return(paste(sample_group, 'done!'))
})


