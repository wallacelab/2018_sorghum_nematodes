#! /usr/bin/Rscript
 
library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input numeric TASSEL genotype file (Reference probability)")
parser$add_argument("-r", "--min-rsq", type='double', default=1, help="Minimum r2 value to be considered redundnat")
parser$add_argument("-o", "--outfile", help="List of redundant markers to be removed")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/Harris_SorghumQtl/")
# args=parser$parse_args(c("-i","1d_genos.imputed.numeric.txt", "-o","99_tmp.txt"))


# Load data
genos=read.delim(args$infile, skip=1, row.names=1)

# Calculate correlations
cat("Calculating all pairwise correlations\n")
cors = cor(genos, use='pairwise')

# Identify redundant markers. (Will always keep the first one found)
cat("Identifying redundant markers\n")
to_remove=list()
for(i in 1:nrow(cors)){
    mycors = cors[i,i:ncol(cors)]   # Correlations from here to end
    redundant = mycors[is.finite(mycors) & mycors >= args$min_rsq]
    if(length(redundant)>1){
        to_remove[[i]] = redundant[-1]  # Skip the first item since that's the marker we're checking against
    }
}
to_remove = unique(names(unlist(to_remove)))    #collapse down to unique names

# Write out
write(to_remove, args$outfile)
