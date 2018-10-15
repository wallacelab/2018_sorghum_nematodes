#! /usr/bin/Rscript
 
library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input results from TASSEL GLM")
parser$add_argument("-g", "--genos", help="Genotype file in numeric format")
parser$add_argument("-t", "--traits", default=NULL, nargs="*", help="Only calculate covariates for these traits")
parser$add_argument("-o", "--outprefix", help="Output file")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/Harris_SorghumQtl/")
# args=parser$parse_args(c("-i","1i_glm.sitefile.txt", "-g", "1d_genos.imputed.numeric.txt", "-o","99_tmp", "-t", "BRIX", "eggperroot_log"))


# Load data
cat("Identifying best markers for covariates\n")
glm=read.delim(args$infile)
genos=read.delim(args$genos, skip=1)
names(genos)[1]="Taxon"

# Filter for specified traits
if(!is.null(args$traits)){
    cat("Filtering for traits",args$traits,"\n")
    glm = subset(glm, glm$Trait %in% args$traits)
    glm=droplevels(glm) # To clear unneded factor levels
}

# Find best marker for each phenotype
traits = split(glm, glm$Trait)
best=lapply(traits, function(mytrait){
    best_p = subset(mytrait, mytrait$p == min(mytrait$p))
    return(best_p$Marker)
})

# Make sure only have one each
n_markers=sapply(best, length)
if(any(n_markers >1)){
    cat("\tWARNING. Trait(s) have >1 best marker:", names(best)[n_markers>1],"\n")
    cat("\t\tOnly the first will be used\n")
    best = lapply(best, function(x){x[1]})
}

# Write out covariate files
cat("\tWriting covariates to file\n")
tmp = lapply(names(best), function(mytrait){
    # Subset genotypes
    mymarker = as.character(best[[mytrait]]) # Make sure is a character, not a factor
    myscores = genos[, names(genos) %in% c("Taxon", mymarker)]
    names(myscores)[names(myscores)==mymarker] = paste(mymarker, "_additive", sep="")
    myscores[,3] = ifelse(myscores[,2]==0.5, yes=1, no=0)
    names(myscores)[3] = paste(mymarker, "_dominant", sep="")
        
    # Write out
    outfile=paste(args$outprefix, mytrait, "txt", sep=".")
    write("<Phenotype>", file=outfile)
    write("taxa\tcovariate\tcovariate", file=outfile, append=T)
    suppressWarnings(write.table(myscores, file=outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T))
    
})
