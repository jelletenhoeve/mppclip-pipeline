
###
### This does the rank tests
###


pathways <- loadPathways()


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	## load the sample data
	mEset <- getDataSet(sample_group, c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'), showWarnings=FALSE)


	## load the mir global test data
	mirGlobalTests <- readObjectOutput(sample_group, name='mir-based_globaltests', pipeline_step='genesets-and-pathways-associations')


	gt_mir_pvals <- t(sapply(mirGlobalTests, function(gto) if(is.null(gto)) NA else result(gto)$'p-value') )
	subsets <- subsets(mirGlobalTests[[1]])
	colnames(gt_mir_pvals) <- names(subsets)
	mirs <- rownames(gt_mir_pvals)


	split_interaction_counts <- interactionCounts(mEset, pathways, split=TRUE, negative.correlations.only=TRUE)
	pathway_pvals_mat <- gt_mir_pvals[, grep('hsa', colnames(gt_mir_pvals))]
	pathway_pvals <- melt(t(pathway_pvals_mat))
	colnames(pathway_pvals) <- c('pathway', 'miRNA', 'Global Test p-value')
	pathway_pvals$'Global Test p-value adjusted FDR' <- p.adjust(pathway_pvals$'Global Test p-value', method='fdr')


	message('Counts...')
	counts <- sapply(split_interaction_counts, function(ct) {
		mapply(as.vector(pathway_pvals$pathway), as.vector(pathway_pvals$miRNA), FUN=function(pw, mir) {
			if(mir %in% rownames(ct) & pw %in% colnames(ct)) {
				ct[mir, pw]
			} else {
				0
			}
		})
	})



	message('Genes...')
	split_interaction_genes <- interactionCounts(mEset, pathways, split=TRUE, return.genes=TRUE, negative.correlations.only=TRUE)
	genes <- sapply(split_interaction_genes, function(gt) {
		mapply(as.vector(pathway_pvals$pathway), as.vector(pathway_pvals$miRNA), FUN=function(pw, mir) {
			if(mir %in% rownames(gt) & pw %in% colnames(gt)) {
				gt[mir, pw]
			} else {
				''
			}
		})
	})
	pathway_size <- sapply(pathways, length)[as.vector(pathway_pvals$pathway)]
	pathway_size_in_data <- sapply(pathways, function(set) {length(intersect(featureNames(eset(mEset, 'mrna')), set))})[as.vector(pathway_pvals$pathway)]

	message('Combining into pathway table...')
	rownames(counts) <- NULL
	rownames(genes) <- NULL
	pw_table <- cbind(pathway_pvals, pathway_size, pathway_size_in_data, counts, genes)


	message(sample_group, " writing pathway count table ...")
	writeTableOutput(pw_table, sample_group, name='pathway_characteristics', pipeline_step='tables-and-figures')

	return(paste(sample_group, 'done!'))
})
