###
### Parclip data
###

#parclip_targets <- read.csv(file='parclip/2010-Jul-12_MCF7B_mir_Targets.csv', stringsAsFactors=FALSE)
parclip_targets <- rbind(read.csv(file='parclip/2010-Jul-12_MCF7B_mir_Targets.csv', stringsAsFactors=FALSE), read.csv(file='parclip/new_targets.csv', stringsAsFactors=FALSE))
parclip_targets <- split(parclip_targets$GeneID, parclip_targets$miRNA)
parclip_targets <- lapply(parclip_targets, unique)

mapping_file <- 'targetscan/miR_Family_Info.txt'
mapping <- read.delim(file=mapping_file, stringsAsFactors=FALSE)
fam2mirbase <- split(mapping$MiRBase.ID, mapping$miR.family)

parclip_targets <- lapply(fam2mirbase, function(x) {
	y <- unique(unlist(parclip_targets[x]))
	y[y != "" & y != "-"]
})
parclip_targets[sapply(parclip_targets, is.null) ] <- NULL

parclip <- do.call('rbind', lapply(names(parclip_targets), function(x) data.frame(miR.Family=x, Gene.Symbol=parclip_targets[[x]], stringsAsFactors=FALSE)))

rm(parclip_targets, mapping_file, mapping, fam2mirbase)