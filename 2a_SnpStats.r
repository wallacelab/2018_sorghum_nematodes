#! /usr/bin/Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input genotype file (hapmap format)")
parser$add_argument("-p", "--phenofile", help="Input phenotype file (TASSEL format)")
parser$add_argument("-s", "--snp", help="SNP to get stats of")
parser$add_argument("-l", "--reverse-log", default=FALSE, action="store_true", help="Whether to also reverse the log-transformation of phenotypes")
parser$add_argument("-o", "--outprefix")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/Harris_SorghumQtl/")
# args=parser$parse_args(c("-i","1d_genos.imputed.hmp.txt", "-p", "1h_good_phenos.tassel.txt", "-o","99_tmp", "-s", "s1221", "-l"))

# Load data
cat("Loading data to calculate SNP stats\n")
genos = read.delim(args$infile)
phenos = read.delim(args$phenofile, skip=2)
phenonames = names(phenos)[-1]

# Pull out target SNP
mysnp = subset(genos, genos$rs. == args$snp)
if (nrow(mysnp) != 1){
    warning("Trying to pull out target SNP did not identify a single row: ",nrow(mysnp)," found instead")
}

# Cut out metadata and simplify
mycalls = mysnp[12:ncol(mysnp)]
mycalls = droplevels(unlist(mycalls))

# Match up to phenotypes
phenos$genotype = mycalls[match(phenos$Taxon, names(mycalls))]

# Graphical output
cat("\tOutputting graphical summary\n")
nplots = length(phenonames)
png(paste(args$outprefix, ".boxplots.png", sep=""), width=5, height=nplots * 4, units="in", res=300)
    par(mfrow=c(nplots, 1), cex=1.2)
    for(trait in phenonames){
        plot(phenos[[trait]] ~ phenos$genotype, xlab=paste("Genotype at", args$snp), ylab=trait)
    }
dev.off()

# Text output
cat("\tOutputting text summary\n")
output = lapply(phenonames, function(trait){
    mydata = tapply(phenos[[trait]], phenos$genotype, summary)  # Get a summary of trait based on genotype levels
    combined = do.call(rbind, mydata)
    combined = data.frame(trait=trait, genotype=rownames(combined), combined)   # Reformat
    if("NA.s" %in% names(combined)){combined$NA.s = NULL}    # Remove column of NAs that only appaers for some traits
    return(combined)
})
output = do.call(rbind, output)

# Undo log transformation if requested
if(args$reverse_log){
    cat("\tIncluding text output with log transformation reversed\n")
    islog = subset(output, grepl("_log$", output$trait))
    islog$trait = paste(islog$trait, "_unlogged", sep="")
    mystats = islog[, 3:ncol(islog)]
    unlogged = data.frame(trait=islog$trait, genotype=islog$genotype, exp(mystats))
    output=rbind(output, unlogged)
}


# Write out
write.table(output, file=paste(args$outprefix, ".summary.txt", sep=""), quote=F, row.names=F, col.names=T, sep='\t')


# Write out phenotype table with genotypes in case want to do manual checks
write.table(phenos, file=paste(args$outprefix, ".raw.txt", sep=""), quote=F, row.names=F, col.names=T, sep='\t')
