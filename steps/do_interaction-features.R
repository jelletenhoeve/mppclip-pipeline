###
### This script generates are target feature table (formely know)
### as the summary count table.
###



#
# Load the external datasets in which the features can be found.
#

for (ext_dataset_name in run_config[['ext_data_features']]) {
	loadExtData(ext_dataset_name, mir_nomenclature=run_config[['mir_nomenclature']], as.list.if.possible=TRUE)
}


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {
	
	## load the data
	mEset <- getDataSet(sample_group)

	## select mirs and genes for which we have data 
	mirna_exprs <- exprs(eset(mEset, id = 'mirna'))
	mrna_exprs  <- exprs(eset(mEset, id = 'mrna'))

	#scaled_mirna_exprs <- t(scale(t(mirna_exprs)))
	#scaled_mrna_exprs <- t(scale(t(abs(mrna_exprs))))
	
	mirs <- rownames(mirna_exprs)
	genes <- rownames(mrna_exprs)

	interaction_features <- get(run_config[['base_interaction_features_table']])
	interaction_features <- interaction_features[interaction_features$mirna %in% mirs & interaction_features$gene %in% genes, ]


	interaction_features <- unique(interaction_features) # keep unique rows only

	# interaction_features can be transcript based,
	# to make it gene based, we select the first transcript as representative.
	# TODO: implement a method with a rationale instead picking a random interaction
	unique.idx <- sapply(split(1:nrow(interaction_features), paste(interaction_features[[1]], interaction_features[[2]], sep=':')), function(x) x[1])
	interaction_features <- interaction_features[unique.idx, ]




	#
	# DATA FEATURES
	# - mirna.mean
	# - mirna.var
	# - mrna.mean
	# - mrna.var
	# - Correlation
	# - Interaction
	# - Interaction2
	#

	message(sample_group, ': adding means and variances...')

	mirna_mean  <- rowMeans(mirna_exprs, na.rm=TRUE)
	mirna_var   <- apply(mirna_exprs, 1, var, na.rm=TRUE)
	mrna_mean  <- rowMeans(mrna_exprs, na.rm=TRUE)
	mrna_var   <- apply(mrna_exprs, 1, var, na.rm=TRUE)

	interaction_features$mirna.mean <- mirna_mean[as.vector(interaction_features$mirna)]
	interaction_features$mirna.var  <- mirna_var[as.vector(interaction_features$mirna)]
	interaction_features$mrna.mean  <- mrna_mean[as.vector(interaction_features$gene)]
	interaction_features$mrna.var   <- mrna_var[as.vector(interaction_features$gene)]


	message(sample_group, ': adding correlation...')

	cm.melt <- melt(correlationMatrix(mEset))
	cor.vec <- cm.melt$value
	names(cor.vec) <- paste(cm.melt$mirna, cm.melt$mrna, sep='|')
	interaction_features$'Correlation' <- cor.vec[paste(interaction_features$mirna, interaction_features$gene, sep='|')]


	message(sample_group, ': adding interaction terms...')
	interaction_features$'Interaction'  <- mapply(as.vector(interaction_features$mirna), as.vector(interaction_features$gene), FUN=function(mir, gene) mean(mrna_exprs[gene, ] * mirna_exprs[mir, ], na.rm=TRUE))
	interaction_features$'Interaction2' <- mapply(as.vector(interaction_features$mirna), as.vector(interaction_features$gene), FUN=function(mir, gene) {
		mean(abs(mrna_exprs[gene, ]) * mirna_exprs[mir, ], na.rm=TRUE)
	})




	#
	# MCF7 features
	#

	if ('mcf7' %in% run_config[['ext_data_features']]) {
		message(sample_group, ': adding mcf7 data...')
		idx <- match(interaction_features$gene, mcf7$gene)
		interaction_features$MCF7_mRNA_RPKM <- mcf7$RPKM[idx]
	}

	#
	# PARCLIP features
	#

	if ('parclip' %in% run_config[['ext_data_features']]) {
		message(sample_group, ': adding parclip columns...')
		interaction_features$'Parclipped'  <- mapply(as.vector(interaction_features$mirna), as.vector(interaction_features$gene), FUN=function(mir, gene) gene %in% parclip[[mir]] )
		interaction_features$'Parclipped.by.any.miRNA' <- mapply(as.vector(interaction_features$mirna), as.vector(interaction_features$gene), FUN=function(mir, gene) gene %in% unlist(parclip))
	}


	message(sample_group, ': writing tables...')
	writeObjectOutput(interaction_features, sample_group=sample_group, name='interaction_features', pipeline_step='interaction-features')
	writeTableOutput(interaction_features, sample_group=sample_group, name='interaction_features', pipeline_step='interaction-features')


	return(paste(sample_group, 'done!'))
})

