#! /usr/bin/Rscript

library(argparse)
library(qtl)
library(parallel)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input R/QTL file (output from step 1p)")
parser$add_argument("-o", "--outprefix")
parser$add_argument("-p", "--perms", type="integer", default=10, help="Number of random permutations to run")
# parser$add_argument("-z", "--zoom", help="Chromosome to zoom in on (currently only accepts 1)")
parser$add_argument("-c", "--num-cores", type="integer", default=1, help="Number of computer cores to use")
parser$add_argument("-s", "--seed", type="integer", default=1, help="Random seed to use")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/Harris_SorghumQtl/")
# args=parser$parse_args(c("-i","1p_rqtl_map.csv", "-o","99_tmp", "-c", "7"))

# Load genotypes
set.seed(args$seed)
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
dim = ceiling(sqrt(nplots))
png(paste(args$outprefix, ".scans.png",sep=""), width=5 * dim, height=3*dim, units="in", res=150)
    par(mfrow=c(dim, dim))
    for(i in 1:nplots){
        
        # Get cutoffs
        permlods = sort(as.numeric(perms[[i]]))
        cutoff_95 = permlods[length(permlods) * 0.95]
        
        # Make plots
        mycim = cims[[i]]
        ymax = max(c(cutoff_95, mycim$lod))
        plot(mycim, ylim=c(0, ymax), main=names(cims)[i])
        abline(h=cutoff_95, col="darkred", lwd=1, lty='dotted')
        
    }
dev.off()



# Make plots
nplots = length(cims)
dim = ceiling(sqrt(nplots))
svg(paste(args$outprefix, ".scans.svg",sep=""), width=5 * dim, height=3*dim)
    par(mfrow=c(dim, dim))
    for(i in 1:nplots){
        
        # Get cutoffs
        permlods = sort(as.numeric(perms[[i]]))
        cutoff_95 = permlods[length(permlods) * 0.95]
        
        # Make plots
        mycim = cims[[i]]
        ymax = max(c(cutoff_95, mycim$lod))
        plot(mycim, ylim=c(0, ymax), main=names(cims)[i])
        abline(h=cutoff_95, col="darkred", lwd=1, lty='dotted')
        
    }
dev.off()

# Write data and perms to get an idea of cutoffs
saveRDS(cims, file=paste(args$outprefix, ".cim_results.rds",sep=""))

cims.out = lapply(names(cims), function(trait){
    data.frame(trait=trait, cims[[trait]])
})
cims.out = do.call(rbind, cims.out)
write.table(cims.out, file=paste(args$outprefix, ".cim_results.txt",sep=""), sep='\t', quote=F, row.names=F, col.names=T)

perms.out = do.call(cbind, perms)
colnames(perms.out) = names(cims)
write.table(perms.out, file=paste(args$outprefix, ".perm_results.txt",sep=""), sep='\t', quote=F, row.names=F, col.names=T)
