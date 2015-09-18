parclip <- rbind(read.csv(file='data/parclip/2010-Jul-12_MCF7B_mir_Targets.csv', stringsAsFactors=FALSE), read.csv(file='data/parclip/new_targets.csv', stringsAsFactors=FALSE))
parclip <- parclip[parclip$GeneID != "" & parclip$GeneID != "-", ]
parclip <- parclip[, c('miRNA', 'GeneID', 'CCR')]


locs <- read.csv(file='data/parclip/PAR_CLIP_genomic_location.csv', stringsAsFactors=FALSE, na.strings=c("NA", "NULL"))[,-1]


parclip <- cbind(parclip, locs[match(parclip$CCR, locs$cluster.name), ])


write.csv(parclip, file='data/parclip/parclip_all_info.csv')

parclip_mirbase_3utr <- unique(parclip[parclip$X3.UTR == 1, c('miRNA', 'GeneID')])
parclip_mirbase_cds  <- unique(parclip[parclip$CDS == 1,    c('miRNA', 'GeneID')])
parclip_mirbase_5utr <- unique(parclip[parclip$X5.UTR == 1, c('miRNA', 'GeneID')])
parclip_mirbase_all  <- unique(parclip[, c('miRNA', 'GeneID')])


save(parclip_mirbase_3utr, file='data/parclip_mirbase_3utr.rda')
save(parclip_mirbase_cds, file='data/parclip_mirbase_cds.rda')
save(parclip_mirbase_5utr, file='data/parclip_mirbase_5utr.rda')
save(parclip_mirbase_all, file='data/parclip_mirbase_all.rda')


