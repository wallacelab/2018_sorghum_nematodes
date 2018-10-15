# setwd('C:\\Users\\jgwall\\Desktop\\000_ToUdate_PGRP\\Harris_SorghumQtl')
phenos=read.csv('0_sorghum_phenos.csv', row.names=1)

# Histograms of phenotypes
ncol=3
nrow=ncol(phenos)
png('1_phenos_hists.png', width=3*ncol, height=3*nrow, units='in', res=150)
  par(mfrow=c(nrow, ncol))
  null=lapply(names(phenos), function(p){
    mydata=phenos[,p]
    hist(mydata, main=p, breaks=40, col='darkblue')
    hist(log(mydata, base=10), main=paste("log",p), breaks=40, col='darkgreen')
    hist(log(mydata+1, base=10), main=paste("log",p,"+1"), breaks=40, col='darkred')
  })
dev.off()
