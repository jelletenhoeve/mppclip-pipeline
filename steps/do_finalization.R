#
# Finalize
# - make a zip of the data for downloading
#
wd <- getwd()
setwd(run_config[['output_dir']])
csv_files <- list.files('data', pattern=".csv", full.names=TRUE, recursive=TRUE)
png_files <- list.files('data', pattern=".png", full.names=TRUE, recursive=TRUE)

zip('download/figures_and_tables.zip', c(png_files, csv_files))

setwd(wd)
