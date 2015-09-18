###
### This script trains the model and makes a prediction 
###

model_training_config <- run_config[['model_training_config']]


mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {


	set.seed(run_config[['seed']])



	message('Loading interaction_features from .rda file')
    interaction_features <- readObjectOutput(sample_group = sample_group, name='interaction_features', pipeline_step='interaction-features', for_shiny=FALSE)    


	#
	# Create a train set by setting a filter (typicall Parclipped by any)
	# And do the prediction for all data we have interatcion features for.
	#
	train.idx <- interaction_features[[model_training_config[['filter_variable']]]]
	train_set <- as.matrix(interaction_features[train.idx, model_training_config[['features']]])

	# remove row containing NA's
	n.train <- nrow(train_set)
	train.idx[which(train.idx)[apply(train_set, 1, function(x) any(is.na(x)))]] <- FALSE
	train_set <- as.matrix(interaction_features[train.idx, model_training_config[['features']]])

	message(n.train - nrow(train_set), ' train interactions were left out due to NA values.')

	# make the response vector
	train_y   <- interaction_features[train.idx, model_training_config[['response_variable']]]

	message('Training glmnet model...')
	fit.glmnet <- cv.glmnet(train_set, train_y,
		family=model_training_config[['family']],
		alpha=model_training_config[['alpha']],
		type.measure=model_training_config[['type.measure']]
	)

	# additinally also train a GLM model which can be used for generating the coefficients table later
	message('Training glm model...')
	fit.glm <- glm(train_y ~ ., family=binomial, data=data.frame(train_set))


	message('Do prediction for all data...')
	newx <- as.matrix(interaction_features[, model_training_config[['features']]])
	
	response <- predict(fit.glmnet, newx, type='response', s=model_training_config[['s']])[,1]
	
	prediction <- response >= model_training_config[['cutoff']]
	prediction[is.na(prediction)] <- FALSE

	# make a comprehensive data frame as output table
	prediction.df <- cbind(
		interaction_features[, c('mirna', 'gene')],
		newx,
		interaction_features[, c(model_training_config[['response_variable']], model_training_config[['filter_variable']])],
		response, prediction
	)


	message(sample_group, ': writing files...')

	writeTableOutput(prediction.df, sample_group, name='prediction', pipeline_step='model-prediction')
	writeObjectOutput(prediction.df, sample_group, name='prediction', pipeline_step='model-prediction')


	writeObjectOutput(fit.glmnet, sample_group, name='glmnet-fit', pipeline_step='model-prediction')
	writeObjectOutput(fit.glm, sample_group, name='glm-fit', pipeline_step='model-prediction')

	writeObjectOutput(train_set, sample_group, name='train-set', pipeline_step='model-prediction')
	writeObjectOutput(train_y, sample_group, name='train-y', pipeline_step='model-prediction')


	return(paste(sample_group, 'done!'))
})


