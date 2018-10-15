# setwd('C:\\Users\\jgwall\\Desktop\\000_ToUdate_PGRP\\Harris_SorghumQtl')
genos=read.delim('1a_genos_filtered.txt')


# Convert to hapmap
hmp = data.frame('rs#' = paste("s", 1:nrow(genos), sep=""))
hmp$alleles = paste(genos$REF, genos$ALT, sep="/")
hmp$chrom=genos$Chr
hmp$pos = genos$Pos
hmp$strand=hmp$'assembly#'=hmp$center=hmp$protLDIS=hmp$assayLDIS=hmp$panelLDIS=hmp$QCcode=NA

# Convert genotype calls
taxa = names(genos)[5:ncol(genos)]
calls = lapply(taxa, function(taxon){
  mycalls = genos[,taxon]
  mycalls[is.na(mycalls)] = -999 # To make subscripting below work; mycalls thrown out afterward so is okay
  alleles = rep(NA, length(mycalls))
  alleles[mycalls==0] = paste(genos$REF, genos$REF, sep="")[mycalls==0]
  alleles[mycalls==0.5] = paste(genos$REF, genos$ALT, sep="")[mycalls==0.5]
  alleles[mycalls==1] = paste(genos$ALT, genos$ALT, sep="")[mycalls==1]
  alleles[mycalls==-999] = "NN"
  return(alleles)
})
calls=as.data.frame(calls)
names(calls) = taxa

# Add in
hmp = cbind(hmp, calls)
write.table(hmp, file='1b_genos_filtered.hmp.txt', sep='\t', quote=F, row.names=F, col.names=T)
