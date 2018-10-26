#! /usr/bin/Rscript

library(argparse)
library(qtl)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of R/QTL CIM scans in RDS format (= saved R object from step 2b))")
parser$add_argument("-o", "--outprefix")
parser$add_argument("-c", "--chrom", help="Which chromosome to get a LOD score on")
parser$add_argument("-t", "--traits", nargs="*", help="Which traits to compute LOD intervals for")
parser$add_argument("-l", "--lod-cutoff", type="double", default=3, help="LOD-support interval cutoff")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/Harris_SorghumQtl/")
# args=parser$parse_args(c("-i","2b_rqtl_analysis.cim_results.rds", "-o","99_tmp", "-c", "05", "-t", "eggperroot_log", "noeggs_log"))

# Load genotypes
cat("Loading R/QTL CIM results from",args$infile,"\n")
cims = readRDS(args$infile)

# Cut down to just CIM results I want
targets = subset(cims, names(cims) %in% args$traits)
cat("\tIdentified",length(targets),"target traits to test\n")

# Compute LOD intervals
lods = lapply(targets, lodint, chr=args$chrom, drop=args$lod_cutoff)

# write out
lods.out = lapply(names(lods), function(trait){
    data.frame(trait, lods[[trait]])
})
lods.out=do.call(rbind, lods.out)
write.table(lods.out, file=paste(args$outprefix, ".lod_intervals.txt", sep=""), quote=F, row.names=F, col.names=T, sep='\t')
