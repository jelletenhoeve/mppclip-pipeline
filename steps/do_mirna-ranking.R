
###
### This does the rank tests
###


pathways <- loadPathways()
genesets <- loadGenesets()


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	## load the sample data
	mEset <- getDataSet(sample_group, c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'), showWarnings=FALSE)


	## load the mir global test data
	mirGlobalTests <- readObjectOutput(sample_group, name='mir-based_globaltests', pipeline_step='genesets-and-pathways-associations')



	gt_mir_pvals <- t(sapply(mirGlobalTests, function(gto) if(is.null(gto)) NA else result(gto)$'p-value') )
	subsets <- subsets(mirGlobalTests[[1]])
	colnames(gt_mir_pvals) <- names(subsets)
	mirs <- rownames(gt_mir_pvals)

	# KEGG mininum pathway minimum p-value
	pathway.min.pvalue <- apply(gt_mir_pvals[, grep('hsa', colnames(gt_mir_pvals))], 1, min)


	interaction_counts <- interactionCounts(mEset, pathways, split=TRUE, return.genes=FALSE, negative.correlations.only=FALSE)

	for (id in interactionNames(mEset)) {
		message('Compiled mirna features for ', id, '...')

		#
		# Micorna features table
		#

		# KEGG targetted mininum pathway minimum p-value
		targetted.pathway.min.pvalue <- sapply(mirs, function(mir) {
			if ( mir %in% rownames(interaction_counts[[id]])) {
				pathways_targetted <- colnames(interaction_counts[[id]])[interaction_counts[[id]][mir, ] > 0]
				if (length(pathways_targetted) == 0) {
					NA
				} else {
					idx <- sapply(pathways_targetted, function(x) which(startsWith(colnames(gt_mir_pvals), x)))
					min(gt_mir_pvals[mir, idx], na.rm=TRUE)
				}
			} else {
				NA
			}
		})


		# Wilcoxon Test pvalue
		wilcoxList <- wilcoxTest(mEset, id=id, split=TRUE, alternative='less')
		#tt.pvals <- sapply(ttList, '[[', 'p.value')
		wilcox.pvals <- sapply(wilcoxList, function(x) {
			if (is.null(x)) {
				NA
			} else {
				x[['p.value']]
			}
		})

		
		wilcox.pvals <- wilcox.pvals[mirs]
		names(wilcox.pvals) <- mirs

		# Construct the table
		mirna_features <- data.frame(matrix(ncol = 0, nrow = length(mirs)))
		rownames(mirna_features) <- mirs

		mirna_features[['Wilcoxon Test p-value']] <- wilcox.pvals
		if (length(which(id == colnames(gt_mir_pvals))) != 1) stop('Internal error..')
		mirna_features[['Global Test miR prediction p-value']] <- gt_mir_pvals[, which(id == colnames(gt_mir_pvals))]

		for (geneset_name in names(genesets)) {
			if ( length(which(geneset_name == colnames(gt_mir_pvals))) != 1 ) stop('Internal error..')
			mirna_features[[paste('Global Test', geneset_name, 'p-value')]] <- gt_mir_pvals[, which(geneset_name == colnames(gt_mir_pvals))]
		}

#		mirna_features[['Global Test minimum p-value of pathways']] <- pathway.min.pvalue
		mirna_features[['Global Test minimum p-value of targetted pathways']] <- targetted.pathway.min.pvalue


		# add the mean, variance, min and maximum, expression
		mirna_exprs <- exprs(eset(mEset, id = 'mirna'))[mirs, ]

		mirna_features$'Expression mean'     <- rowMeans(mirna_exprs, na.rm=TRUE)
		mirna_features$'Expression variance' <- apply(mirna_exprs, 1, var, na.rm=TRUE)
		mirna_features$'Expression min'      <- apply(mirna_exprs, 1, min, na.rm=TRUE)
		mirna_features$'Expression max'      <- apply(mirna_exprs, 1, max, na.rm=TRUE)


		# add pheno global tests if they exist
		if (existsObjectOutput(sample_group, name='phenotype_globaltests', pipeline_step='phenotype-associations')) {

			# Add clinical Global Test results
			message("Loading phenotype global tests...")
			phenoGlobalTests <- readObjectOutput(sample_group, name='phenotype_globaltests', pipeline_step='phenotype-associations')

			pheno.pvals <- sapply(names(phenoGlobalTests), function(phenotype) {
				idx <- grep(id, names(subsets(phenoGlobalTests[[phenotype]])))
				pvals <- result(phenoGlobalTests[[phenotype]])[idx, 'p-value']
				names(pvals) <- substring(names(subsets(phenoGlobalTests[[phenotype]]))[idx], nchar(id)+2)
				pvals[mirs]
			})

			rownames(pheno.pvals) <- mirs
			pheno.pvals <- data.frame(pheno.pvals)

			mirna_features <- cbind(mirna_features, pheno.pvals)

		} else {
			message('No phenotype global tests found for ', sample_group)
		}

		message(sample_group, " writing mirna features table ...")
		writeTableOutput(mirna_features, sample_group, name=paste('mirna_features', id, sep='-'), pipeline_step='mirna-ranking')

		# trick to make column names have . instead of [spaces]
		mirna_features <- readTableOutput(sample_group, name=paste('mirna_features', id, sep='-'), pipeline_step='mirna-ranking')


		message(sample_group, " writing pca-biplot of the mirna features ...")

		startPngOutput(sample_group, name=paste('pca_biplot' , id, sep='-'), pipeline_step='mirna-ranking', width=1240, height=1240)
		pcdat <- sapply(mirna_features, function(x) rank(x) / nrow(mirna_features))		
		rownames(pcdat) <- rownames(mirna_features)
		loadings <- prcomp(pcdat)$rotation
		biplot(prcomp(pcdat, scale=TRUE))
		dev.off()


		# make the ranking tables
		for(rt_name in names(run_config[['ranking_tables']])) {


			if (!all(run_config[['ranking_tables']][[rt_name]] %in% colnames(mirna_features))) {
				warning('Not all columns for ', rt_name, ' are present in "mirna_features')
				next
			}

			cols <- unlist(run_config[['ranking_tables']][[rt_name]][run_config[['ranking_tables']][[rt_name]] %in% colnames(mirna_features)])
			rank_table <- mirna_features[, cols, drop=FALSE]
			names(rank_table) <- names(cols)

			rank_table <- cbind(rank_table, rank.test(rank_table, TRUE), mirna_features[, c('Expression.mean', 'Expression.variance', 'Expression.min', 'Expression.max')])
			rank_table <- rank_table[order(rank_table$'rank-test-p-value'), ]

			message(sample_group, " writing mirna ranking table ...")
			writeTableOutput(rank_table, sample_group, name=paste('rank_table', rt_name, id, sep='-'), pipeline_step='mirna-ranking')
		}


	}


	return(paste(sample_group, 'done!'))
})
