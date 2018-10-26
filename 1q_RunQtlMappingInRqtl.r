#! /usr/bin/Rscript

library(argparse)
library(qtl)
library(parallel)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input R/QTL file (output from step 1p)")
parser$add_argument("-o", "--outprefix")
parser$add_argument("-p", "--perms", type="integer", default=10, help="Number of random permutations to run")
parser$add_argument("-z", "--zoom", help="Chromosome to zoom in on (currently only accepts 1)")
parser$add_argument("-c", "--num-cores", type="integer", default=1, help="Number of computer cores to use")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/Harris_SorghumQtl/")
# args=parser$parse_args(c("-i","1p_rqtl_map.csv", "-o","99_tmp", "-z", "05", "-c", "7"))

# Load genotypes
cat("Loading R/QTL map file from",args$infile,"\n")
map=read.cross(file=args$infile, format="csv", genotypes=c("AA","AB","BB"))
cat("\tLoaded",totmar(map),"genotypes and",nphe(map),"phenotypes\n")

# Do QTL mapping via composite interval mapping (CIM)
cat("Performing composite interval mapping in parallel on",args$num_cores,"cores\n")
map = calc.genoprob(map)
cims = mclapply(phenames(map)[-1], function(mytrait){
    cim(map, pheno.col=mytrait)
}, mc.cores=args$num_cores)
names(cims) = names(map$pheno)[-1]

# Run permutations
cat("Running",args$perms,"permutations in parallel on",args$num_cores,"cores\n")
perms = mclapply(phenames(map)[-1], function(mytrait){
    cim(map, pheno.col=mytrait, n.perm=args$perms) 
}, mc.cores=args$num_cores)

# Make plots
nplots = length(cims)
png(paste(args$outprefix, ".scans.png",sep=""), width=5 * nplots, height=6, units="in", res=150)
    par(mfcol=c(2, nplots))
    for(i in 1:nplots){
        
        # Get cutoffs
        permlods = sort(as.numeric(perms[[i]]))
        cutoff_95 = permlods[length(permlods) * 0.95]
        cutoff_99 = permlods[length(permlods) * 0.99]
        
        # Make plots
        mycim = cims[[i]]
        ymax = max(c(cutoff_99, mycim$lod))
        plot(mycim, ylim=c(0, ymax), main=names(cims)[i])
        abline(h=cutoff_95, col="blue", lwd=1, lty='dotted')
        abline(h=cutoff_99, col="red", lwd=1, lty='dashed')
        
        plot(mycim, ylim=c(0, ymax), chr=args$zoom, main=names(cims)[i])    # Hardwired hack to zoom in on chromo
        abline(h=cutoff_95, col="blue", lwd=1, lty='dotted')
        abline(h=cutoff_99, col="red", lwd=1, lty='dashed')
    }
dev.off()

# 
# # LOD interval
# lodint(cims, chr="05", lodcolumn=3, drop=3)

# Results of the above
#       chr pos      BRIX eggperroot_categories eggperroot_log noeggs_categories noeggs_log   Rootwt
# s1220  05 155 0.6146643              31.51105       32.27009          31.51105   32.88746 1.138557
# s1221  05 160 0.6755093              36.65542       41.04529          36.65542   42.06997 1.455162
# s1222  05 165 0.6176709              34.14317       35.91302          34.14317   37.93194 1.488995

# LOD interval
