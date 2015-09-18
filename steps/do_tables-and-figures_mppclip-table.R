
###
### This does the rank tests
###


loadExtData('mcf7', mir_nomenclature=run_config[['mir_nomenclature']])


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	## load the sample data
	mEset <- getDataSet(sample_group, c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'), showWarnings=FALSE)


	## select mirs and genes for which we have data 
	mirna_exprs <- exprs(eset(mEset, id = 'mirna'))
	mrna_exprs  <- exprs(eset(mEset, id = 'mrna'))

	#scaled_mirna_exprs <- t(scale(t(mirna_exprs)))
	#scaled_mrna_exprs <- t(scale(t(abs(mrna_exprs))))
	
	mirs <- rownames(mirna_exprs)
	genes <- rownames(mrna_exprs)

	mppclip_features <- interactions(mEset, 'model_prediction_parclip_union')



	message(sample_group, ': adding parclip and model predicted columns...')
	pclip.list <- split(interactions(mEset, 'parclip')[['gene']], interactions(mEset, 'parclip')[['mirna']])
	mp.list <- split(interactions(mEset, 'model_prediction')[['gene']], interactions(mEset, 'model_prediction')[['mirna']])

	mppclip_features$'Parclipped'       <- mapply(as.vector(mppclip_features[['mirna']]), as.vector(mppclip_features[['gene']]), FUN=function(mir, gene) gene %in% pclip.list[[mir]] )
	mppclip_features$'Model predicted'  <- mapply(as.vector(mppclip_features[['mirna']]), as.vector(mppclip_features[['gene']]), FUN=function(mir, gene) gene %in% mp.list[[mir]] )



	message(sample_group, ': adding means and variances...')

	mirna_mean  <- rowMeans(mirna_exprs, na.rm=TRUE)
	mirna_var   <- apply(mirna_exprs, 1, var, na.rm=TRUE)

	mrna_mean  <- rowMeans(mrna_exprs, na.rm=TRUE)
	mrna_var   <- apply(mrna_exprs, 1, var, na.rm=TRUE)


	mppclip_features$mirna.mean <- mirna_mean[as.vector(mppclip_features[['mirna']])]
	mppclip_features$mirna.var  <- mirna_var[as.vector(mppclip_features[['mirna']])]
	mppclip_features$mrna.mean  <- mrna_mean[as.vector(mppclip_features[['gene']])]
	mppclip_features$mrna.var   <- mrna_var[as.vector(mppclip_features[['gene']])]




	message(sample_group, ': adding mcf7 data...')
	idx <- match(mppclip_features$gene, mcf7$gene)
	mppclip_features$MCF7_mRNA_RPKM <- mcf7$RPKM[idx]


	message(sample_group, ': adding correlation...')
	cm.melt <- melt(correlationMatrix(mEset))
	cor.vec <- cm.melt$value
	names(cor.vec) <- paste(cm.melt$mirna, cm.melt$mrna, sep='|')
	mppclip_features$'Correlation' <- cor.vec[paste(mppclip_features[['mirna']], mppclip_features[['gene']], sep='|')]



	message(sample_group, ': adding interaction terms...')
	mppclip_features$'Interaction'  <- mapply(as.vector(mppclip_features[['mirna']]), as.vector(mppclip_features[['gene']]), FUN=function(mir, gene) mean(mrna_exprs[gene, ] * mirna_exprs[mir, ], na.rm=TRUE))
	mppclip_features$'Interaction2' <- mapply(as.vector(mppclip_features[['mirna']]), as.vector(mppclip_features[['gene']]), FUN=function(mir, gene) {
		mean(abs(mrna_exprs[gene, ]) * mirna_exprs[mir, ], na.rm=TRUE)
	})

	
	message(sample_group, ': writing tables...')
	writeTableOutput(mppclip_features, sample_group, name='mppclip_targets', pipeline_step='tables-and-figures')



	return(paste(sample_group, 'done!'))
})
