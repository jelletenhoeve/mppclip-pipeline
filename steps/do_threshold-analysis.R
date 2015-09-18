

mEset <- getDataSet(interactions_sets = c('targetscan_predicted_conserved'), showWarnings=FALSE)


message('mirna_abundance_thresholds...')
startPngOutput(name='mirna_abundance_thresholds', pipeline_step='threshold-analysis', width=960, height=960)
thresholdVariationPlot(mEset, id='targetscan_predicted_conserved',
  fixed.eset               = 'mrna',
  fixed.assayDataElement   = 'exprs',
  fixed.threshold          = 0,
  varying.assayDataElement = 'exprs',
  n.varying.thresholds     = 20,
  alternative='less',
  aggregation.function = mean,
  y = c('shift', 'p-value', 'n.varying')
)
dev.off()


message('mirna_variance_thresholds...')
startPngOutput(name='mirna_variance_thresholds', pipeline_step='threshold-analysis', width=960, height=960)
thresholdVariationPlot(mEset, id='targetscan_predicted_conserved',
  fixed.eset               = 'mrna',
  fixed.assayDataElement   = 'exprs',
  fixed.threshold          = 0,
  varying.assayDataElement = 'exprs',
  n.varying.thresholds     = 20,
  alternative='less',
  aggregation.function = var,
  y = c('shift', 'p-value', 'n.varying')
)
dev.off()


message('gene_abundance_thresholds...')
startPngOutput(name='gene_abundance_thresholds', pipeline_step='threshold-analysis', width=960, height=960)
thresholdVariationPlot(mEset, id='targetscan_predicted_conserved',
  fixed.eset               = 'mirna',
  fixed.assayDataElement   = 'exprs',
  fixed.threshold          = 0,
  varying.assayDataElement = 'exprs',
  n.varying.thresholds     = 20,
  alternative='less',
  aggregation.function = mean,
  y = c('shift', 'p-value', 'n.varying')
)
dev.off()


message('gene_variance_thresholds...')
startPngOutput(name='gene_variance_thresholds', pipeline_step='threshold-analysis', width=960, height=960)
thresholdVariationPlot(mEset, id='targetscan_predicted_conserved',
  fixed.eset               = 'mirna',
  fixed.assayDataElement   = 'exprs',
  fixed.threshold          = 0,
  varying.assayDataElement = 'exprs',
  n.varying.thresholds     = 20,
  alternative='less',
  aggregation.function = var,
  y = c('shift', 'p-value', 'n.varying')
)
dev.off()

