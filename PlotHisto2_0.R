#! /usr/bin/Rscript
myArgs <- commandArgs(trailingOnly = TRUE)
x11()
n_of_inputs=(length(myArgs))
files_names <- c()
fileintergenicShuffle = c()
numberIntergenicshuffle = c()
n=c()
n_total = 0
numberIntergenicshuffleFinal = 0
names_columns = c("ChromoInter","BeginInter","EndInter","NameInter",'Name1','Name2')
count = 1
# READ THE FILES WHICH EACH TWO ARGUMENTS ARE THE NAME OF THE FOLDER AND THE NAME OF THE FILE IN THE FOLDER
for (x in seq(3,n_of_inputs,2)) {
  files_names[count]=paste('./', myArgs[x],'/',myArgs[x+1],'.csv', sep="")
  fileintergenicShuffle[count]=paste('./', myArgs[x],'/',myArgs[x],'_intergenic.bed', sep="")
  if(file.exists(fileintergenicShuffle[count])){
    intergenics_shuffles=read.table(fileintergenicShuffle[count])
    colnames(intergenics_shuffles)[1:6]=names_columns
    suppressMessages(attach(intergenics_shuffles))
    numberIntergenicshuffle[count]=nlevels(NameInter)
    n[count]=1
  } else {numberIntergenicshuffle[count]=0
  n[count]=0}
  count = count + 1
}
count = count - 1
filenormal=paste('./', myArgs[1],'/',myArgs[2],'.csv', sep="")

###### GET THE AVERAGE NUMBER OF INTERGENIC PEAKS FROM SHUFFLE AND ORIGINAL ####
for (x in 1:count) {numberIntergenicshuffleFinal = numberIntergenicshuffleFinal + numberIntergenicshuffle[x]
n_total = n_total + n[x]
}

if(n_total==0){
  n_total = 1}
numberIntergenicshuffleFinal = numberIntergenicshuffleFinal/n_total

fileintergenicNormal=paste('./', myArgs[1],'/',myArgs[1],'_intergenic.bed', sep="")
intergenicNormal=read.table(fileintergenicNormal)
colnames(intergenicNormal)[1:6]=c("ChromoInter","BeginInter","EndInter","NameInter",'Name1','Name2')
suppressMessages(attach(intergenicNormal))
numberIntergenicNormal=nlevels(NameInter)

###### INTERGENIC PROPORTION BETWEEN ORIGINAL AND SHUFFLE ####
  
intergenicProporcion=numberIntergenicNormal/numberIntergenicshuffleFinal
 
###### GET THE AVERAGE NUMBER OF INTRAGENIC PEAKS IN SHUFFLE AND ORIGINAL ####


# if(!file.exists(fileshuffle) & !file.exists(fileshuffle2) & !file.exists(fileshuffle3)){
 # print('No intragenic peaks in shuffles')
#} else {
count2 = 0
for (x in files_names){
  if (file.exists(x)){
    count2 = count2 + 1
    tableoffrequencyshuffle=read.table(textConnection(gsub("\"", "\t", gsub(",", "\t", readLines(x)))),header = T)
    suppressMessages(attach(tableoffrequencyshuffle))
    if (count2==1){
    frecPartsshuffleTotal=c(table(Upstream)[2],table(UTR5)[2],table(Exons)[2],table(Introns)[2],table(UTR3)[2],table(Downstream)[2])
    frecPartsshuffleTotal[is.na(frecPartsshuffleTotal)] <- 0
    formafinal = c(count2, numberIntergenicshuffle[count2],frecPartsshuffleTotal)
    } else {frecPartsshuffle=c(table(Upstream)[2],table(UTR5)[2],table(Exons)[2],table(Introns)[2],table(UTR3)[2],table(Downstream)[2])
    frecPartsshuffle[is.na(frecPartsshuffle)] <- 0
    frecPartsshuffleTotal =frecPartsshuffleTotal + frecPartsshuffle 
    formafinal = c(formafinal, count2, numberIntergenicshuffle[count2],frecPartsshuffle)
    }
  } else {frecPartsshuffle=c(0,0,0,0,0,0)
  frecPartsshuffleTotal =frecPartsshuffleTotal + frecPartsshuffle
  }
  shuffle_table = matrix(c("Shuffle","Intergenic","Upstream","5UTR","Exon","Intron","3UTR",'Downstream',formafinal), nrow = 1+count2, ncol = 8, byrow = F)
  filaname=paste('./', myArgs[1],'/',myArgs[2],'shuffletable', sep="")
  write(shuffle_table,file = filaname,ncolumns = 8,sep='\t')
}

averagefreqshuffle=frecPartsshuffleTotal/count2
  
  if (file.exists(filenormal)){
    tableoffrequencyNormal=read.table(textConnection(gsub("\"", "\t", gsub(",", "\t", readLines(filenormal)))),header = T)
    suppressMessages(attach(tableoffrequencyNormal))
    frecPartsNormal=c(table(Upstream)[2],table(UTR5)[2],table(Exons)[2],table(Introns)[2],table(UTR3)[2],table(Downstream)[2])
    frecPartsNormal[is.na(frecPartsNormal)] <- 0
    resultfreq=frecPartsNormal/averagefreqshuffle
  } else {  
    print(paste('There is no intragenic peaks in the ','./', myArgs[1],'/',myArgs[2],'.csv'))
    frecPartsNormal=c(0,0,0,0,0,0)
    }

###### CALCULATE RESULTS AND PLOT THE BAR GRAPH OF RELATIVE PROPORTION ####
  
resultfreq=frecPartsNormal/averagefreqshuffle
resultfreq=replace(resultfreq,resultfreq==Inf,0)
resultfreq=replace(resultfreq,is.nan(resultfreq),0)
resultfreqinter=c(intergenicProporcion,resultfreq)

filaname=paste('./', myArgs[1],'/',myArgs[2],'tableFoldChange', sep="")
matrixresults=matrix(c("Intergenic","Upstream","5UTR","Exon","Intron","3UTR",'Downstream',resultfreqinter,numberIntergenicshuffleFinal,averagefreqshuffle,numberIntergenicNormal,frecPartsNormal),nrow=4,ncol=7,byrow = F)
write(matrixresults,file = filaname,ncolumns = 7,sep='\t')
