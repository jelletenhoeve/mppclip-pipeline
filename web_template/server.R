library(shiny)

stopifnot(packageVersion('shiny') > '0.7')

#library(Biobase)


#source('../lib/AllClasses.R')
#source('../lib/AllGenerics.R')
#source('../lib/MatchedExpressionSet.R')
#source('../lib/helpers.R')


# Server Data loading functions 




# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output, session) {


  output$mirna_abundance_thresholds <- renderImage({
    getImageSrc(name='mirna_abundance_thresholds', pipeline_step='threshold-analysis', for_shiny=TRUE)    
  }, deleteFile = FALSE)
  output$mirna_variance_thresholds <- renderImage({
    getImageSrc(name='mirna_variance_thresholds', pipeline_step='threshold-analysis', for_shiny=TRUE)    
  }, deleteFile = FALSE)
  output$gene_abundance_thresholds <- renderImage({
    getImageSrc(name='gene_abundance_thresholds', pipeline_step='threshold-analysis', for_shiny=TRUE)    
  }, deleteFile = FALSE)
  output$gene_variance_thresholds <- renderImage({
    getImageSrc(name='gene_variance_thresholds', pipeline_step='threshold-analysis', for_shiny=TRUE)    
  }, deleteFile = FALSE)



  output$ts_interaction_features <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name='ts_interaction_features', pipeline_step='tables-and-figures', for_shiny=TRUE)
  })

  output$mppclip_targets <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name='mppclip_targets', pipeline_step='tables-and-figures', for_shiny=TRUE)    
  })
  output$correlation_analysis <- renderImage({
    getImageSrc(sample_group = input$sample_group, name='correlation_analysis', pipeline_step='tables-and-figures', for_shiny=TRUE)    
  }, deleteFile = FALSE)


  output$mirna_ranking2 <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name=paste('rank_table', input$ranking_table, input$target_set, sep='-'), pipeline_step='mirna-ranking', for_shiny=TRUE)    
  })

  output$mirna_features <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name=paste('mirna_features', input$target_set, sep='-'), pipeline_step='mirna-ranking', for_shiny=TRUE)    
  })

  output$pca_biplot <- renderImage({
    getImageSrc(sample_group = input$sample_group, name=paste('pca_biplot', input$target_set, sep='-'), pipeline_step='mirna-ranking', for_shiny=TRUE)    
  }, deleteFile = FALSE)


  output$pathway_characteristics <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name='pathway_characteristics', pipeline_step='tables-and-figures', for_shiny=TRUE)
  })

  output$pathway_clustering <- renderImage({
    getImageSrc(sample_group = input$sample_group, name='pathway_clustering', pipeline_step='genesets-and-pathways-associations', for_shiny=TRUE)
  }, deleteFile = FALSE)


  output$mirna_ranking <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name=paste('mirna_ranking', input$target_set, sep='-'), pipeline_step='phenotype-associations', for_shiny=TRUE)
  })


  output$phenotype_globaltests <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name='phenotype_globaltests', pipeline_step='phenotype-associations', for_shiny=TRUE)    
  })

  output$mir_based_globaltests <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name='mir-based_globaltests', pipeline_step='genesets-and-pathways-associations', for_shiny=TRUE)    
  })

  output$prediction <- renderDataTable({
    #readTableOutput(sample_group = input$sample_group, name='prediction', pipeline_step='model-prediction', for_shiny=TRUE)    
    readObjectOutput(sample_group = input$sample_group, name='prediction', pipeline_step='model-prediction', for_shiny=TRUE)    
  })

  output$pairs_overlap <- renderImage({
      getImageSrc(sample_group = input$sample_group, name='pairs_overlap', pipeline_step='model-prediction', for_shiny=TRUE)
  }, deleteFile = FALSE)

  output$coefficients <- renderDataTable({
    readTableOutput(sample_group = input$sample_group, name='coefficients', pipeline_step='model-prediction', for_shiny=TRUE)    
  })

  output$glmnet_prediction_measures <- renderImage({
      getImageSrc(sample_group = input$sample_group, name='glmnet_prediction_measures', pipeline_step='model-prediction', for_shiny=TRUE)
  }, deleteFile = FALSE)


  # sizes
  output$n_samples <- renderText({
    nrow(pData(mEset()))
  })
  output$n_mirnas <- renderText({
    nrow(eset(mEset(), 'mirna'))
  })
  output$n_genes <- renderText({
    nrow(eset(mEset(), 'mrna'))
  })

  mEset <- reactive({
    getDataSet(sample_group = input$sample_group, for_shiny=TRUE)
  })
 

  #
  # Downloads
  #
  output$download_figures_and_tables <- downloadHandler(
    filename = function() { paste('mppclip_figures_and_tables.zip', sep='') },
    content = function(file) {
      file.copy('download/figures_and_tables.zip', file)
    },
    contentType = "application/zip"
  )
})

