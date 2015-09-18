#
# Preparation
# - make the output directory
# - write the pipeline configs to a file
# - make a copy of the shinyR application files, copy web_template
#

dir.create(run_config[['output_dir']], showWarnings=FALSE, recursive=TRUE)
dir.create(paste(run_config[['output_dir']], '/data', sep=''), showWarnings=FALSE, recursive=TRUE)
dir.create(paste(run_config[['output_dir']], '/download', sep=''), showWarnings=FALSE, recursive=TRUE)
dir.create(paste(run_config[['output_dir']], '/log', sep=''), showWarnings=FALSE, recursive=TRUE)
save(global_config, run_config, file=paste(run_config[['output_dir']], '/config.rda', sep=''))
file.copy(from='web_template/global.R', to=run_config[['output_dir']], overwrite=TRUE)
file.copy(from='web_template/server.R', to=run_config[['output_dir']], overwrite=TRUE)
file.copy(from='web_template/ui.R', to=run_config[['output_dir']], overwrite=TRUE)
file.copy(from='web_template/column_names.txt', to=run_config[['output_dir']], overwrite=TRUE)
