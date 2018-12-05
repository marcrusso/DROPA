#! /usr/bin/Rscript
myArgs <- commandArgs(trailingOnly = TRUE)
x11()
fileUnder=paste('./', myArgs[1],'/',myArgs[3],'tableFoldChange', sep="")
fileOver=paste('./', myArgs[1],'/',myArgs[2],'tableFoldChange', sep="")

###### GET THE AVERAGE NUMBER OF INTERGENIC PEAKS FROM SHUFFLE AND ORIGINAL ####
if(file.exists(fileUnder)){
UnderTable=read.table(fileUnder,header = T)
} else {UnderTable = matrix(0,nrow=3,ncol = 3)}
if(file.exists(fileUnder)){
OverTable=read.table(fileOver,header = T)
} else {OverTable = matrix(0,nrow=3,ncol = 3)}

totalperRegionShuffle=c(UnderTable[2,1],UnderTable[2,2]+OverTable[2,2],UnderTable[2,3]+OverTable[2,3],UnderTable[2,4]+OverTable[2,4],UnderTable[2,5]+OverTable[2,5],UnderTable[2,6]+OverTable[2,6],UnderTable[2,7]+OverTable[2,7])

totalperRegionNormal=c(UnderTable[3,1],UnderTable[3,2]+OverTable[3,2],UnderTable[3,3]+OverTable[3,3],UnderTable[3,4]+OverTable[3,4],UnderTable[3,5]+OverTable[3,5],UnderTable[3,6]+OverTable[3,6],UnderTable[3,7]+OverTable[3,7])

totalProportion=totalperRegionNormal/totalperRegionShuffle

totalProportion=replace(totalProportion,totalProportion==Inf,0)
totalProportion=replace(totalProportion,is.nan(totalProportion),0)

filename=paste('./', myArgs[1],'/Enrichment_over_expected.table', sep="")

matrixresults=matrix(c('Region',"Intergenic","Upstream","5UTR","Exon","Intron","3UTR",'Downstream','Results /Shuffle',totalProportion,'Shuffle',totalperRegionShuffle,'Results',totalperRegionNormal),nrow=4,ncol=8,byrow = T)

write(matrixresults,file = filename,ncolumns = 4,sep='\t')
namegraph=paste('./', myArgs[1],'/Enrichment_over_expected.pdf', sep="")
pdf(namegraph)
barplot(totalProportion,main=paste("Relative proportion of peaks in gene regions of Total Peaks"),xlab="Gene Region",ylab="Input Peaks/Shuffled Peaks",col = rainbow(length(totalProportion)),names.arg=c("Intergenic","Upstream","5UTR","Exon","Intron","3UTR",'Downstream'),cex.names=0.9)
invisible(dev.off())


###### SUM THE SHUFFLE PEAKS BETWEEN UNDER AND OVER
filename=paste('./', myArgs[1],'/Totalshuffletable', sep="")
fileUnder=paste('./', myArgs[1],'/',myArgs[3],'shuffletable', sep="")
fileOver=paste('./', myArgs[1],'/',myArgs[2],'shuffletable', sep="")
if(file.exists(fileUnder)){
  UnderTable=read.table(fileUnder,header = T)
  if(file.exists(fileOver)){
  OverTable=read.table(fileOver,header = T)
  shufflelist = list(UnderTable,OverTable)
  totalperRegionShuffle = Reduce('+',shufflelist)
  totalperRegionShuffle['Shuffle']=totalperRegionShuffle['Shuffle']/2
  totalperRegionShuffle['Intergenic']=totalperRegionShuffle['Intergenic']/2
  totalperRegionShuffle = as.data.frame(do.call(cbind, totalperRegionShuffle))
  write.table(totalperRegionShuffle,file = filename,sep='\t', row.names = FALSE)
  } else {
    UnderTable=read.table(fileUnder,header = T)
    write.table(UnderTable,file = filename,sep='\t', row.names = FALSE)
  }
} else {if(file.exists(fileOver)) {    
  OverTable=read.table(fileOver,header = T)
  write.table(OverTable,file = filename,sep='\t', row.names = FALSE)
  }}

