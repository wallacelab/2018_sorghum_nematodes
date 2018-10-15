# setwd('C:\\Users\\jgwall\\Desktop\\000_ToUdate_PGRP\\Harris_SorghumQtl')
genos=read.csv('0_sorghum_genos.csv', comment.char="#", header=T)

# Filter parameters
max_het_sites = 0.55

# Change to more standard encoding
metadata = genos[,1:4]
calls = genos[,-c(1:4)]
calls[calls==0] = NA
calls[calls==3] = 1.5
calls = calls-1 # To put in a 0-0.5-1 range for convenience

# Filter by heterozygosity
hets_sites = rowSums(calls == 0.5, na.rm=T) / ncol(calls)
tokeep = hets_sites <= max_het_sites


# TODO: Can add more filtering if later decide to; this is all needed for now

# Do filtering
calls = calls[tokeep,]
metadata = metadata[tokeep,]

# Write filtered data
newgenos = cbind(metadata, calls)
write.table(newgenos, file='1a_genos_filtered.txt', sep='\t', quote=F, row.names=F, col.names=T)
