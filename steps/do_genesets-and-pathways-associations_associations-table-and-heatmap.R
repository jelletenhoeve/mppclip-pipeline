##
## This script outputs:
## - a table with the globaltest outcomes of the tartgets (model-prediction, paclip, targetscan), genesets, and pathways
## - a heatmap with the globatest for pathways
##



genesets <- loadGenesets()
pathways <- loadPathways()


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	mirGlobalTests <- readObjectOutput(sample_group, name='mir-based_globaltests', pipeline_step='genesets-and-pathways-associations')

	gts.mat <- do.call('rbind', lapply(mirGlobalTests, function(res) res@result))
	gts.df <- data.frame(
		mir                = do.call('c', lapply(names(mirGlobalTests), function(mir) rep(mir, nrow(mirGlobalTests[[mir]]@result)))),
		subset             = rownames(gts.mat),
		'gt-p.value'       = gts.mat[, 'p-value'],
		'gt-statistic'     = gts.mat[, 'Statistic'], 
		'gt-expected'      = gts.mat[, 'Expected'],
		'gt-stdv'          = gts.mat[, 'Std.dev'],
		'gt-n'             = gts.mat[, '#Cov'],
		'p.value.adjusted' = p.adjust(gts.mat[, 'p-value'], method='fdr'),
		stringsAsFactors   = FALSE
	)


	# write outputs
	# full table
	message(sample_group, ': writing files...')
	writeTableOutput(gts.df, sample_group, name='mir-based_globaltests', pipeline_step='genesets-and-pathways-associations')


	# the heatmap
	hm.m <- sapply(mirGlobalTests, function(res) result(res)$'p-value')
	rownames(hm.m) <- rownames(result(mirGlobalTests[[1]]))

	hm.m <- hm.m[grep('hsa', rownames(hm.m)), ]                   # only select pathways
	hm.m <- hm.m[apply(hm.m, 1, function(x) !all(is.na(x))), ]    # remove pathway for which no p-value were calculated


	startPdfOutput(sample_group, name='pathway_clustering', pipeline_step='genesets-and-pathways-associations', width=15, height=15)
	heatmap.2(log10(hm.m), margins=c(15,20), cexRow=.6, cexCol=.6, lhei=c(1, 10), lwid=c(1, 5), main=sample_group, trace='none', symbreaks=FALSE)
	dev.off()

	startPngOutput(sample_group, name='pathway_clustering', pipeline_step='genesets-and-pathways-associations', width=1280, height=1280)
	heatmap.2(log10(hm.m), margins=c(15,20), cexRow=.6, cexCol=.6, lhei=c(1, 10), lwid=c(1, 5), main=sample_group, trace='none', symbreaks=FALSE)
	dev.off()


	return(paste(sample_group, 'done!'))
})

