###
### This script generates coefficients of the trained models
###

model_training_config <- run_config[['model_training_config']]



#
# Start the loop
#

mc.result <- mclapply(run_config[['sample_groups']], mc.cores=run_config[['n_cores']], function(sample_group) {



	message('Loading fitted models and train_set from .rda file')
    fit.glmnet <- readObjectOutput(sample_group = sample_group, name='glmnet-fit', pipeline_step='model-prediction', for_shiny=FALSE)    
    fit.glm <- readObjectOutput(sample_group = sample_group, name='glm-fit', pipeline_step='model-prediction', for_shiny=FALSE)    
    train_set <- readObjectOutput(sample_group = sample_group, name='train-set', pipeline_step='model-prediction', for_shiny=FALSE)    
    train_y <- readObjectOutput(sample_group = sample_group, name='train-y', pipeline_step='model-prediction', for_shiny=FALSE)    



	## Performance plots
#	message('Writing ', glmnet_train_performance_plots_file, '...')
#	pdf(file=glmnet_train_performance_plots_file, width=10.5, height=7)
	startPngOutput(sample_group, name='glmnet_prediction_measures', pipeline_step='model-prediction', width=960, height=640)
	par(mfrow=c(2,3))
	glmnet_train_response          <- predict(fit.glmnet, train_set, s=model_training_config[['s']], type='response')[,1]
	glmnet_train_prediction_object <- plotPredictionMeasures(glmnet_train_response, train_y, cutoff=model_training_config[['cutoff']])
	dev.off()

	startPngOutput(sample_group, name='glm_prediction_measures', pipeline_step='model-prediction', width=960, height=640)
	par(mfrow=c(2,3))
	glm_train_response             <- predict(fit.glm, data.frame(train_set), type='response')
	glm_train_prediction_object    <- plotPredictionMeasures(glm_train_response, train_y, cutoff=model_training_config[['cutoff']])
	dev.off()


	## Univariate analyses
	message('Doing univariate analyses...')

	features <- model_training_config[['features']]

	glmnet.multivar.coef <- coef(fit.glmnet, s = model_training_config[['s']])[features,]

	glm.multivar.coef <- summary.glm(fit.glm)$coefficients[, 'Estimate'][features]
	glm.multivar.pval <- summary.glm(fit.glm)$coefficients[, 'Pr(>|z|)'][features]
	
	univar.auc <- sapply(features, function(f) {
		prediction_object <- prediction(train_set[, f], train_y)
		performance(prediction_object, "auc")@y.values[[1]]
	})

	univar.correlation <- sapply(features, function(f) {
		cor(train_set[, f], train_y, method='pearson', use='pairwise.complete.obs')
	})

	univar.correlation.pval <- sapply(features, function(f) {
		cor.test(train_set[, f], as.numeric(train_y), method='pearson', use='pairwise.complete.obs')$p.value
	})

	mean.parclip <- sapply(features, function(f) {
		mean(train_set[train_y, f], na.rm=TRUE)
	})
	mean.non.parclip <- sapply(features, function(f) {
		mean(train_set[!train_y, f], na.rm=TRUE)
	})

	coef.df <- data.frame(
		variable                    = features,
		mean.parclip                = mean.parclip,
		mean.non.parclip            = mean.non.parclip,
		univar.correlation          = univar.correlation,
		univar.correlation.pval     = univar.correlation.pval,
		univar.correlation.pval.adj = p.adjust(univar.correlation.pval, method='fdr'),
		univar.auc                  = univar.auc,
		glmnet.multivar.coef        = glmnet.multivar.coef,
		glmnet.included             = abs(glmnet.multivar.coef) > 0,
		glm.multivar.coef           = glm.multivar.coef,
		glm.multivar.pval           = glm.multivar.pval
	)

	glmnet_train.auc <- performance(glmnet_train_prediction_object, "auc")@y.values[[1]]
	glm_train.auc <- performance(glm_train_prediction_object, "auc")@y.values[[1]]

	train.auc.df <- data.frame(
		variable                    = 'AUC ON TRAIN SET',
		mean.parclip                = '',
		mean.non.parclip            = '',
		univar.correlation          = '',
		univar.correlation.pval     = '',
		univar.correlation.pval.adj = '',
		univar.auc                  = '',
		glmnet.multivar.coef        = glmnet_train.auc,
		glmnet.included             = '',
		glm.multivar.coef           = glm_train.auc,
		glm.multivar.pval           = ''
	)

	coef.df <- rbind(coef.df, train.auc.df)

	message('Writing coefficients file...')
	
	message(sample_group, ': writing files...')
	
	rownames(coef.df) <- NULL
	
	writeTableOutput(coef.df, sample_group, name='coefficients', pipeline_step='model-prediction')



	return(paste(sample_group, 'done!'))
})


