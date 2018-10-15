# setwd('C:\\Users\\jgwall\\Desktop\\000_ToUdate_PGRP\\Harris_SorghumQtl')
phenos=read.csv('0_sorghum_phenos.csv', row.names=1)

log_trans = names(phenos) # Log-transform all phenotypes
logistic_trans = c('noeggs', 'eggperroot') #Turn into 0/1 phenotypes
category_trans = c('noeggs', 'eggperroot') #turn into 0/1/2 phenotypes
cat_cutoffs = list(noeggs = 20000, eggperroot = 1000)# Cutoffs between categories 1 & 2

# Make new data frame of phenotypes
#newphenos=data.frame(row.names=rownames(phenos))
newphenos=phenos# Decided to keep raw values in there just as a check for heritability

# Log transform
for(myname in log_trans){
  newname=paste(myname, '_log', sep="")
  newphenos[,newname] = log(phenos[,myname])
  # Set values from log-transforming 0s to missing
  infinite = !is.finite(newphenos[,newname])
  newphenos[infinite, newname]= NA
}

# Logistic transform
for(myname in logistic_trans){
  newname=paste(myname, '_binary', sep="")
  newphenos[,newname] = ifelse(phenos[,myname]>0, yes=1, no=0)
}

# Category transform
for(myname in category_trans){
  newname=paste(myname, '_categories', sep="")
  # Separte into categories; 0=0, middle=1, high=2
  upper_category=ifelse(phenos[,myname]>cat_cutoffs[[myname]], yes=2, no=1)
  newphenos[,newname] = ifelse(phenos[,myname]==0, yes=0, no=upper_category)
}

# Sort & write out
newphenos = newphenos[,order(names(newphenos))]
write.table(newphenos, file='1a_phenos_transformed.txt', sep='\t', quote=F, row.names=T, col.names=T)

# Write out in TASSEL format
outfile="1a_phenos_transformed.tassel.txt"
write("<Phenotype>", outfile)
write(paste(c('taxa', rep("data", ncol(newphenos))), collapse='\t'), outfile, append=T)
write(paste(c('Taxon', names(newphenos)), collapse='\t'), outfile, append=T)
write.table(newphenos, file=outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

# Histograms of phenotypes, split by original phenotype on rows
orig = sub(names(newphenos), pattern="_.+", rep="")
ncol=max(table(orig))
nrow=length(unique(orig))
mycolor=1

png('1_phenos_hists.transformed.png', width=3*ncol, height=3*nrow, units='in', res=150)
  par(mfrow=c(nrow, ncol))
  for(o in unique(orig)){
    mycolor=mycolor+1
    mydata = subset(newphenos, select = orig==o)
    null=lapply(names(mydata), function(p){
      hist(mydata[,p], main=p, breaks=40, col=mycolor)
    })
    # Add blank spots to right of plots if needed so next pheno starts on next line
    blank_count = ncol - sum(orig==o)
    if(blank_count>0){
      for(i in 1:blank_count){plot(NA, NA, axes=F, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")}
    }
  }
dev.off()
