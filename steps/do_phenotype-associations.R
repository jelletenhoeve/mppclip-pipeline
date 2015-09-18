###
### This script trains the model and makes a prediction 
###




mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	## load the data
	mEset <- getDataSet(sample_group, c('targetscan_predicted_conserved', 'parclip', 'model_prediction', 'model_prediction_parclip_union'))

	# do global tests for phenotypes
	message('Starting pheno global tests')

	phenoGlobalTests <- globalTestPheno(mEset, run_config[['pheno_columns']], nullPhenoColumns=run_config[['null_pheno_columns']], include.mirna=TRUE)
	phenoGlobalTests[sapply(phenoGlobalTests, is.null)] <- NULL


	if (length(phenoGlobalTests) == 0) {
		message('No clinical global tests for ', sgID)
	} else {


		gts.mat <- do.call('rbind', lapply(phenoGlobalTests, function(res) res@result))
		gts.df <- data.frame(
			phenotype          = do.call('c', lapply(names(phenoGlobalTests), function(pheno) rep(pheno, nrow(phenoGlobalTests[[pheno]]@result)))),
			subset             = rownames(gts.mat),
			'gt-p.value'       = gts.mat[, 'p-value'],
			'gt-statistic'     = gts.mat[, 'Statistic'], 
			'gt-expected'      = gts.mat[, 'Expected'],
			'gt-stdv'          = gts.mat[, 'Std.dev'],
			'gt-n'             = gts.mat[, '#Cov'],
			'p.value.adjusted' = p.adjust(gts.mat[, 'p-value'], method='fdr'),
			stringsAsFactors   = FALSE
		)

		# Add multiple test adjusted pvalue's per target set	

#	 	adj.pvals <- lapply(interactionNames(mEset), function(id) {
#			idx <- startsWith(gts.df$subset, paste(id, '.', sep=''))
#			p.value.adjusted.targetset <- rep(NA, nrow(gts.df))
#			
#			p.value.adjusted.targetset[idx] <- p.adjust(gts.df$gt.p.value[idx], method='fdr')
#			p.value.adjusted.targetset
#	 	})

#	 	names(adj.pvals) <- paste('p.value.adjusted', interactionNames(mEset), sep='.')

#		gts.df <- cbind(gts.df, data.frame(adj.pvals))	 	

		adj.pvals.per.id <- rep(NA, nrow(gts.df))
		for (id in interactionNames(mEset)) {
			idx <- startsWith(gts.df$subset, paste(id, '.', sep=''))
			adj.pvals.per.id[idx] <- p.adjust(gts.df$gt.p.value[idx], method='fdr')
	 	}
		gts.df[['adjusted.pvalue.per.id']] <- adj.pvals.per.id

	
		message(sample_group, ': writing files')
	
		writeObjectOutput(phenoGlobalTests, sample_group, name='phenotype_globaltests', pipeline_step='phenotype-associations')
		
		writeTableOutput(gts.df, sample_group, name='phenotype_globaltests', pipeline_step='phenotype-associations')


		# make a ranking table (for the phenotype characteristics only)
		for (id in interactionNames(mEset)) {

			pheno.pvals <- sapply(names(phenoGlobalTests), function(phenotype) {
				idx <- startsWith(names(subsets(phenoGlobalTests[[phenotype]])), paste(id, '.', sep=''))
				pvals <- result(phenoGlobalTests[[phenotype]])[idx, 'p-value']
				names(pvals) <- substring(names(subsets(phenoGlobalTests[[phenotype]]))[idx], nchar(id)+2)
				pvals
			})

			pheno.pvals <- data.frame(pheno.pvals)

			mirna_ranking <- cbind(pheno.pvals, rank.test(pheno.pvals, TRUE))
			
			# Sort by p-value
			mirna_ranking <- mirna_ranking[order(mirna_ranking$'rank-test-p-value'), ]

			message(sample_group, " writing rank table ...")
			writeTableOutput(mirna_ranking, sample_group, name=paste('mirna_ranking', id, sep='-'), pipeline_step='phenotype-associations')
		}
	}


	return(paste(sample_group, 'done!'))
})

