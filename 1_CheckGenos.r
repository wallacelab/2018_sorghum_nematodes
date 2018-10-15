# setwd('C:\\Users\\jgwall\\Desktop\\000_ToUdate_PGRP\\Harris_SorghumQtl')
genos=read.csv('0_sorghum_genos.csv', comment.char="#", header=T)

# According to comment, 0=missing, 1=ref, 2=alt, 3=het. REcode these to be additive
metadata = genos[,1:4]
calls = genos[,-c(1:4)]
calls[calls==0] = NA
calls[calls==3] = 1.5
calls = calls-1 # To put in a 0-0.5-1 range for convenience

# MDS of genotypes
mydists = dist(t(calls))
mds = cmdscale(mydists, eig=T)

# # # # # LD using simple distance matrix - Never got working quite right; just did in TASSEL
# # # # ld_calls = calls
# # # # ld_calls[ld_calls==0.5] = NA # Set hets to missing
# # # # ld = cor(t(ld_calls), use='pairwise')

# Missingness
missing_taxa= colSums(is.na(calls)) / nrow(calls)
missing_sites= rowSums(is.na(calls)) / ncol(calls)

# Heterozygosity
hets_taxa = colSums(calls == 0.5, na.rm=T) / nrow(calls)
hets_sites = rowSums(calls == 0.5, na.rm=T) / ncol(calls)


# Plot
nrow=5
ncol=2
png('1_genos_plots.png', width=3*ncol, height=3*nrow, units='in', res=150)
  par(mfrow=c(nrow,ncol))
  
  # MDS plot
  percents = round(mds$eig / sum(mds$eig) * 100, digits=1) # % variance explained
  parents=subset(mds$points, grepl(rownames(mds$points), pattern="^ENTRY"))
  plot(mds$points[,1], mds$points[,2], col="darkblue", pch=20, cex=1.5, main='MDS plot',
       xlab=paste("PC1 -", percents[1],"%"), ylab=paste("PC2 -", percents[2],"%"))
  points(parents[,1], parents[,2], col="red", pch=20, cex=2)
  plot(NA, NA, axes=F, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")  # Blank plot for filler
  
# # #   # LD heatmap - Never got working quite right.
# # #   # heatmap(as.matrix(ld), Rowv=NA, Colv=NA, labRow=NA, labCol=NA)
# # #   # plot(NA, NA, axes=F, ylim=c(0,1), xlim=c(0,1), ylab="", xlab="")# Blank plot
  
  # Missingness - Taxa
  hist(missing_taxa, col="darkred", main="Missingness - Taxa", breaks=40, ylab="Fraction missing")
  plot(sort(missing_taxa), col="darkred", main="Missingness - Taxa", xlab="entry", ylab="Fraction missing")
    
  # HEterozygosity - Taxa
  hist(hets_taxa, col="darkgreen", main="Heterozygosity - Taxa", breaks=40, ylab="Fraction het")
  plot(sort(hets_taxa), col="darkgreen", main="Heterozygosity - Taxa", xlab="entry", ylab="Fraction het")
  
  # Missingness - Sites
  hist(missing_sites, col="red", main="Missingness - Sites", breaks=40, ylab="Fraction missing")
  plot(sort(missing_sites), col="red", main="Missingness - Sites", xlab="sites", ylab="Fraction missing")
  
  # HEterozygosity - Sites
  hist(hets_sites, col="green", main="Heterozygosity - Sites", breaks=40, ylab="Fraction het")
  plot(sort(hets_sites), col="green", main="Heterozygosity - Sites", xlab="sites", ylab="Fraction het")
  
  
    
dev.off()
