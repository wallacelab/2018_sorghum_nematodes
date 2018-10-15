#! /usr/bin/Rscript

library(argparse)
library(qtl)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input ABH genotype file (converted via TASSEL & in csv format)")
parser$add_argument("-p", "--phenos", help="Phenotype file in TASSEL format")
parser$add_argument("-b", "--badmarkers", nargs="*",  help="List of known bad markers to drop")
parser$add_argument("-o", "--outprefix")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/Harris_SorghumQtl/")
# args=parser$parse_args(c("-i","1o_genos_imputed.abh.csv", "-o","99_tmp", "-p", "1h_good_phenos.tassel.txt"))

# Format phenotype file
phenofile=paste(args$outprefix, ".tmp_phenos.csv", sep="")
phenos=read.delim(args$phenos, skip=2)
names(phenos)[names(phenos)=='Taxon'] = "id"    # For rQTL to identify
write.csv(phenos, file=phenofile, row.names=F)

# Load genotypes
cat("Loading genotype data from",args$infile,"\n")
genos=read.cross(file=args$infile, format="csvs", phefile=phenofile)
cat("\tLoaded",totmar(genos),"genotypes\n")

# Drop supplied bad markers
if(length(args$badmarkers)>0){
    cat("Dropping",length(args$badmarkers),"bad markers supplied by user\n")
    genos = drop.markers(genos, args$badmarkers)
}

# Drop possible problematic markers due to weird segregation ratios
summary <- geno.table(genos)
todrop <- rownames(summary[summary$P.value < 0.05/totmar(genos),])  # Drop markers with Bonferroni-corrected p-value below 0.05
genos <- drop.markers(genos, todrop)
cat("Dropped", length(todrop),"markers with segregation distortion (Bonferroni-corrected p-values < 0.05)\n")

# Calculate genetic map
cat("Calculating pairwise recombination frequencies\n")
genos = est.rf(genos)
cat("Calculating genetic map (this will take a while\n")
newmap=est.map(genos, map.function="kosambi", verbose=T, error.prob=0.01, maxit=100, n.cluster=7)
png(paste(args$outprefix, ".maps.png",sep=""), width=1000, height=2000)
    par(mfrow=c(2,1))
    plotRF(genos)
    plotMap(newmap, show.marker.names=T)
dev.off()

# Write out new map
genos = replace.map(genos, newmap)
write.cross(genos, filestem=args$outprefix)


# # Do QTL mapping
# scan = scanone(genos, n.cluster=7, pheno.col=2:ncol(genos$pheno))
# nplots = ncol(genos$pheno)-1
# png(paste(args$outprefix, ".scans.png",sep=""), width=5, height=3 * nplots, units="in", res=150)
#     par(mfrow=c(nplots,1))
#     for(i in 1:nplots){
#         plot(scan, lodcolumn=i)
#     }
# dev.off()
# 
# # LOD interval
# lodint(scan, chr="05", lodcolumn=3, drop=3)

# Results of the above
#       chr pos      BRIX eggperroot_categories eggperroot_log noeggs_categories noeggs_log   Rootwt
# s1220  05 155 0.6146643              31.51105       32.27009          31.51105   32.88746 1.138557
# s1221  05 160 0.6755093              36.65542       41.04529          36.65542   42.06997 1.455162
# s1222  05 165 0.6176709              34.14317       35.91302          34.14317   37.93194 1.488995

# LOD interval
