###
### Pipeline execution settings
###


run_config <- list(

	# processing parameters
	n_cores = 5,
	seed = 1234,
	output_dir = '../w_shiny/tcga_tsfam_1e-04_0',
	


	# data set parameters
	dataset_name     = 'TCGA',
	mir_nomenclature = 'targetscan_family',
	mir_threshold    = 1e-04,
	gene_threshold   = 0,


	sample_groups = c(
		'basal','her2', 'luminala', 'luminalb', 'all'
	),
	# removed "LumA" because it hangs at genesets and pathways associations

	# interaction feature parameters
	base_interaction_features_table = 'targetscan',  # should be any of 'ext_data_features'
	ext_data_features = list(
		'targetscan',
		'mcf7',
		'parclip'
	),
	
	# model prediction parameters
	model_training_config = list(
		filter_variable = 'Parclipped.by.any.miRNA',
		features = c(
			"Total.num.conserved.sites",
			"Number.of.conserved.8mer.sites",
			"Number.of.conserved.7mer.m8.sites",
			"Number.of.conserved.7mer.1a.sites",
			"Total.num.nonconserved.sites",  
			"Number.of.nonconserved.8mer.sites",
			"Number.of.nonconserved.7mer.m8.sites",
			"Number.of.nonconserved.7mer.1a.sites",
			"Total.context.score",
			"Aggregate.PCT",
			"mirna.mean",
			"mirna.var",
			"mrna.mean",
			"mrna.var",
			"Interaction",
			"Interaction2"
		),
		response_variable = 'Parclipped',
		alpha = 0.5,
		s = 'lambda.min',
		cutoff = 0.5,
		family = 'binomial',
		type.measure = 'deviance'
	),


	# genesets and pathway assocation features
	genesets = c('CGC', '70genes'),
	pathways = 'kegg_pathways',


	# phenotype associations features
	pheno_columns     = list(
		overall_survival = c('overall_survival', 'event_death')
	),

	null_pheno_columns  = list(
#		survival    = c('size_mm', 'grade', 'lymph_nodes_positive'),
#		metastasis  = c('size_mm', 'grade', 'lymph_nodes_positive')
	),

	ranking_tables = list(
		table1 = list(
			target_activity     = 'Global.Test.miR.prediction.p.value',
			pathway_activity    = 'Global.Test.minimum.p.value.of.targetted.pathways',
			cancer_activity     = 'Global.Test.CGC.p.value',
			prognostic_activity = 'overall_survival'
		),
		table2 = list(
			target_activity     = 'Global.Test.miR.prediction.p.value',
			pathway_activity    = 'Global.Test.minimum.p.value.of.targetted.pathways',
			cancer_activity     = 'Global.Test.CGC.p.value'
		),
		table3 = list(
			target_activity     = 'Global.Test.miR.prediction.p.value',
			pathway_activity    = 'Global.Test.minimum.p.value.of.targetted.pathways'
		)
	)

)
