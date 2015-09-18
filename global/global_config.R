


global_config <- list(


	pipeline_steps = c(
		'threshold-analysis',
		'interaction-features',
		'model-prediction',
		'genesets-and-pathways-associations',
		'phenotype-associations',
		'mirna-ranking',
		'tables-and-figures'
	),

	target_sets = c(
		'targetscan_predicted_conserved',
		'parclip',
		'model_prediction',
		'model_prediction_parclip_union'
	),

	feature_sets = list(
		targetscan_features = c(
			"Total.num.conserved.sites",
			"Number.of.conserved.8mer.sites",
			"Number.of.conserved.7mer.m8.sites",
			"Number.of.conserved.7mer.1a.sites",
			"Total.num.nonconserved.sites",  
			"Number.of.nonconserved.8mer.sites",
			"Number.of.nonconserved.7mer.m8.sites",
			"Number.of.nonconserved.7mer.1a.sites",
			"Total.context.score",
			"Aggregate.PCT"
		),
		miranda_features = c(
		)
	)

)