# Title     : TODO
# Objective : TODO
# Created by: bruno
# Created on: 07.07.18
#! /usr/bin/Rscript

# Plotting the graph that shows the frecuency of peak in each zone of the genes.
# ExInput: myArgs[1]=FGloss myArgs[2]=Over_FPKMFGloss
x11()
myArgs <- commandArgs(trailingOnly = TRUE)
#myArgs=c('Normal','Under_FPKMNormal')
library(UpSetR)
a=paste('./', myArgs[1],'/',myArgs[2],'.csv', sep="")
tableoffrequency=read.table(textConnection(gsub("\"", "\t", gsub(",", "\t", readLines(a)))),header = T)
attach(tableoffrequency)
colu=c()
Downstream=factor(Downstream)
if (nlevels(Downstream)==2) {
  colu=c(colu,'Downstream')
}
UTR3=factor(UTR3)
if (nlevels(UTR3)==2) {
  colu=cbind(colu,'UTR3')
}
Introns=factor(Introns)
if (nlevels(Introns)==2) {
  colu=c(colu,'Introns')
}
Exons=factor(Exons)
if (nlevels(Exons)==2) {
  colu=c(colu,'Exons')
}
UTR5=factor(UTR5)
if (nlevels(UTR5)==2) {
  colu=c(colu,'UTR5')
}
Upstream=factor(Upstream)
if (nlevels(Upstream)==2) {
  colu=c(colu,'Upstream')
}
numberofsets=length(colu)
movies <- read.csv(a, header = T, sep = "\t",",")
a=paste('./', myArgs[1],'/',myArgs[2],'_Upsetplot.pdf', sep="")
pdf(a)
upset(movies,numberofsets,sets =colu,keep.order = T,order.by = "freq")
invisible(dev.off())

##### Plotting the Pie Chart for the frecuency of Intergenic and Intragenic peaks. #####

NameInter2=paste('./', myArgs[1],'/',myArgs[1],'_intergenic.bed', sep="")
if (file.exists(NameInter2)) {
  results2=read.table(NameInter2)
  colnames(results2) [1:6] = c("ChromoInter","BeginInter","EndInter","NameInter",'Name1','Name2')
  suppressMessages(attach(results2))
  a=paste('./', myArgs[1],'/',myArgs[1],'PeaksInGenes', sep="")
  results=read.table(a)
  colnames(results) [1:8] = c("Chromo","Begin","End","NameGen",'Strand','NamePeak','BeginPeak','EndPeak')
  suppressMessages(attach(results))
  sizegenes=c(nlevels(NamePeak),nlevels(NameInter))
  piepercent<- round(100*sizegenes/sum(sizegenes), 1)
  genname=c(paste('Intragenic:',nlevels(NamePeak)),paste('Intergenic:',nlevels(NameInter)))
  # pie(sizegenes,labels = paste(piepercent,'%'), main='Percent of Intergenic and Intragenic peaks',col = rainbow(length(sizegenes)))
  # legend("topright",genname, cex = 0.8,fill = rainbow(length(sizegenes)))
  # a=paste('./', myArgs[1],'/',myArgs[1],'Intergenic.png', sep="")
  # savePlot(a)
  } else {print('No Intergenic Regions')}

##### Pie charts of the Overlapping and No Overlapping Genes. ####

namefile=paste('./',myArgs[1],'/',myArgs[2],'.csv',sep='')
dat=read.table(textConnection(gsub("\"", "\t", gsub(",", "\t", readLines(namefile)))),header = T)
suppressMessages(attach(dat))
frecParts=c(table(Upstream)[2],table(UTR5)[2],table(Exons)[2],table(Introns)[2],table(UTR3)[2],table(Downstream)[2])
frecParts[is.na(frecParts)] <- 0


##### Bar charts of Intergenic and parts of the Genes. #####

if (file.exists(NameInter2)) {
  namefile=paste('./',myArgs[1],'/',myArgs[2],'.csv',sep='')
  dat=read.table(textConnection(gsub("\"", "\t", gsub(",", "\t", readLines(namefile)))),header = T)
  suppressMessages(attach(dat))
  frecParts=c(table(Upstream)[2],table(UTR5)[2],table(Exons)[2],table(Introns)[2],table(UTR3)[2],table(Downstream)[2],nlevels(NameInter))

  frecParts[is.na(frecParts)] <- 0

  piepercent<- round(100*frecParts/sum(nlevels(NamePeak),nlevels(NameInter)), 1)
  NumbersInLegend=c(paste('Upstream:',table(Upstream)[2]),paste('UTR5:',table(UTR5)[2]),paste('Exons:',table(Exons)[2]),paste('Introns:',table(Introns)[2]),paste('UTR3:',table(UTR3)[2]),paste('Downstream:',table(Downstream)[2]),paste('Intergenic:',nlevels(NameInter)))
  nametitle='Frequency of peaks in parts of the genes'
  adressfile=paste('./', myArgs[1],'/',myArgs[2],'_Histogram.pdf', sep="")
  pdf(adressfile)
  bar <- barplot(piepercent,names.arg = c("Upstream", "5UTR", "Exon", "Intron", "3UTR", "Downstream", "Intergenic"),cex.names = 0.8, ylab = "Percent of Total Overlaps Peaks", main=nametitle,col = rainbow(length(frecParts)), axisnames = TRUE, ylim = c(0,80))
  text(bar, piepercent, labels = piepercent, pos = 3)
  legend('topright',NumbersInLegend, cex = 0.8,fill = rainbow(length(frecParts)))
  invisible(dev.off())
  }

namefile=paste('./',myArgs[1],'/',myArgs[2],'.csv',sep='')
dat=read.table(textConnection(gsub("\"", "\t", gsub(",", "\t", readLines(namefile)))),header = T)
suppressMessages(attach(dat))
frecParts=c(table(Upstream)[2],table(UTR5)[2],table(Exons)[2],table(Introns)[2],table(UTR3)[2],table(Downstream)[2],nlevels(NameInter))

frecParts[is.na(frecParts)] <- 0

piepercent<- round(100*frecParts/sum(frecParts), 1)
NumbersInLegend=c(paste('Upstream:',table(Upstream)[2]),paste('5UTR:',table(UTR5)[2]),paste('Exon:',table(Exons)[2]),paste('Intron:',table(Introns)[2]),paste('3UTR:',table(UTR3)[2]),paste('Downstream:',table(Downstream)[2]),paste('Intergenic:',nlevels(NameInter)))
nametitle=paste('Frequency of peaks in parts of the genes')
adressfile=paste('./', myArgs[1],'/',myArgs[2],'Intergenic.pdf', sep="")
pdf(adressfile)
pie(frecParts,labels = paste(piepercent,'%'), main=nametitle,col = rainbow(length(frecParts)))
legend("bottomleft",NumbersInLegend, cex = 0.8,fill = rainbow(length(frecParts)))
invisible(dev.off())