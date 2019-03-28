#! /usr/bin/Rscript


from pathlib import Path
import argparse
import os
import DROPAcore as core
import warnings
import shutil
warnings.filterwarnings("ignore")

############ INPUT OF VARIABLES ###################

parser = argparse.ArgumentParser(description='Input program')
parser.add_argument('namefile', type=str, help='Query BED file')
parser.add_argument('-ref', type=str, help='Gene Reference Folder')
parser.add_argument('-o',type=str,help= 'Output directory')
parser.add_argument('-ex',type=str, help='Expression file. Default is "None"')
parser.add_argument('-lim',type=float, help='Expression Threshold. Default is 0.5')
parser.add_argument('-dis',type=int, help='Distance from TSS and TTS, to define Upstream and Downstream regions. Default is 5000')
parser.add_argument('-shuffle',type=int , help='Numbers for a shuffle comparison, 0 for none. Default is 0')
parser.add_argument('-gsize',type=str , help='Genome size file, used for shuffle')
parser.set_defaults(rnaseq='None',li=0.5,dis=5000,shuffle=False)

data=parser.parse_args().namefile
data2=parser.parse_args().o
rnaSeq=parser.parse_args().ex
distance=parser.parse_args().dis
limit=parser.parse_args().lim
shuffle=parser.parse_args().shuffle
annotation = parser.parse_args().ref
chromosome_size_file = parser.parse_args().gsize


################### BODY OF THE PROGRAM ###################
Utr3 = str(str(annotation) + '/'+'3UTR.bed')
Utr5 = str(str(annotation) + '/'+'5UTR.bed')
CDE = str(str(annotation) + '/'+'CDE.bed')
annotation = str('./' + str(annotation) + '/'+'refgene.bed')
if limit == None:
    limit = 0.5
if shuffle == None:
    shuffle = 0

gainresults = core.PeakOverlap(annotation,data,tssdistance=distance,peakname=data2)

checkresults=core.CheckExpression(gainresults[0],rnaSeq,limit=limit,peakname=data2,TSSTTSdistance=distance)

try:
    FeatureAssignOverResults=core.FeatureAssign(checkresults[0],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2,lap='Expressed_')
    ExonOverResults= core.TableCreator(FeatureAssignOverResults,namepeak=data2,lap='Expressed')
    fileOver = "./" + data2 + "/Expressed_" + data2 +".csv"
    core.Upsetplotting('./' + data2 + '/' + data2 + '_Expressed_Annotation.table',"Expressed_"+data2, data2)
    core.Dropa_pie_AllwithIntergenic('./' + data2 + '/' + data2 + '_intergenic.bed',
                                      './' + data2 + '/' + data2 + "_Expressed_Annotation.table", "Expressed_" + data2, data2, False)
except:
    print('There is no peaks in expressed genes')
    fileOver = "None"

try:
    FeatureAssignUnderResults=core.FeatureAssign(checkresults[1],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2,lap='NotExpressed_')
    ExonUnderResults=core.TableCreator(FeatureAssignUnderResults,namepeak=data2,lap='NotExpressed')
    fileUnder = "./" + data2 + "/NotExpressed_" + data2 +".csv"
    core.Upsetplotting('./' + data2 + '/' + data2 + '_NotExpressed_Annotation.table', 'NotExpressed_'+data2, data2)
    core.Dropa_pie_AllwithIntergenic('./' + data2 + '/' + data2 + '_intergenic.bed',
                                      './' + data2 + '/' + data2 + "_NotExpressed_Annotation.table", "NotExpressed_" + data2, data2, False)
except:
    print('There is no peaks in unexpressed genes')
    fileUnder = "None"

#### Make a plot of the total percent of each zone been overlap by a peak for express and not express genes together
h,k=1,1
if fileOver != "None":
    if fileUnder != "None":
        filenames = [fileOver, fileUnder]
        filenames2 = ['./' + data2 + '/' + data2 + "_Expressed" + '_Annotation.table',
                     './' + data2 + '/' + data2 + "_NotExpressed" + '_Annotation.table']
        with open('./' + data2 + '/' + data2 + '_All_Annotation.table', 'w') as outfile:
            for fname in filenames2:
                with open(fname) as infile:
                    for line in infile:
                        if h==1:
                            outfile.write(line)
                        else:
                            h=1
                    h=0
        with open("./" + data2 + "/All_" + data2 + ".csv", 'w') as outfile:
            h=1
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        if h==1:
                            outfile.write(line)
                        else:
                            h=1
                    h=0

    else:
        filenames2 = ['./' + data2 + '/' + data2 + "_Expressed" + '_Annotation.table']
        with open('./' + data2 + '/' + data2 + '_All_Annotation.table', 'w') as outfile:
            for fname in filenames2:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
        filenames = fileOver
        with open("./" + data2 + "/All_" + data2 + ".csv", 'w') as outfile:
                with open(filenames) as infile:
                    for line in infile:
                        outfile.write(line)

elif fileUnder != "None":
    filenames2 = ['./' + data2 + '/' + data2 + "_NotExpressed" + '_Annotation.table']
    with open('./' + data2 + '/' + data2 + '_All_Annotation.table', 'w') as outfile:
        for fname in filenames2:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    filenames = fileUnder
    with open("./" + data2 + "/All_" + data2 + ".csv", 'w') as outfile:
            with open(filenames) as infile:
                for line in infile:
                    outfile.write(line)

core.Upsetplotting('./' + data2 + '/' + data2 + '_All_Annotation.table',"All"+data2, data2)

core.Dropa_pie_Intergenic('./' + data2 + '/' + data2 + '_intergenic.bed','./' + data2 + '/' + data2 + 'PeaksInGenes',data2, data2)

core.Dropa_pie_AllwithIntergenic('./' + data2 + '/' + data2 + '_intergenic.bed','./' + data2 + '/' + data2 + "_All_Annotation.table","All_"+data2, data2, True)

core.Dropa_histogram('./' + data2 + '/' + data2 + '_intergenic.bed','./' + data2 + '/' + data2 + "_All_Annotation.table", data, "All_"+data2, data2)

############### MAKE THE COMPARATION WITH SHUFFLE #######################
if shuffle != 0:
    for x in range(0, shuffle):
        shufflepeaks=core.randpeak(data,x+1,data2,chromosome_size_file)

        resultsgain=core.PeakOverlap(annotation,shufflepeaks,tssdistance=distance,peakname=data2+'shuffle'+str(x+1))

        if os.stat(resultsgain[0]).st_size == 0:
            print('Shuffle1 has no Intergenic Peaks')
        else:
            resultscheck=core.CheckExpression(resultsgain[0],rnaSeq,limit=limit,peakname=data2+'shuffle'+str(x+1),TSSTTSdistance=5000)

        try:
            resultsFeatureAssignOver=core.FeatureAssign(resultscheck[0],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2+'shuffle'+str(x+1),lap='Expressed')

            finalResultsOver=core.TableCreator(resultsFeatureAssignOver,namepeak=data2+'shuffle'+str(x+1),lap='Expressed')

        except:
            print('There is no peaks in express genes in Shuffle '+ str(x+1))

        try:
            resultsFeatureAssignUnder=core.FeatureAssign(resultscheck[1],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2+'shuffle'+str(x+1),lap='NotExpressed')

            finalResultsUnder=core.TableCreator(resultsFeatureAssignUnder,namepeak=data2+'shuffle'+str(x+1),lap='NotExpressed')

        except:
            print('There is no peaks in unexpress genes in Shuffle ' + str(x+1))
        h, k = 1, 1
        if fileOver != "None":
            if fileUnder != "None":
                filenames = [fileOver, fileUnder]
                filenames2 = ['./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + "_Expressed" + '_Annotation.table',
                              './' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + "_NotExpressed" + '_Annotation.table']
                with open('./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + '_All_Annotation.table', 'w') as outfile:
                    for fname in filenames2:
                        with open(fname) as infile:
                            for line in infile:
                                if h == 1:
                                    outfile.write(line)
                                else:
                                    h = 1
                            h = 0
                with open("./" + data2+'shuffle'+str(x+1) + "/All_" + data2+'shuffle'+str(x+1) + ".csv", 'w') as outfile:
                    h = 1
                    for fname in filenames:
                        with open(fname) as infile:
                            for line in infile:
                                if h == 1:
                                    outfile.write(line)
                                else:
                                    h = 1
                            h = 0

            else:
                filenames2 = ['./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + "_Expressed" + '_Annotation.table']
                with open('./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + '_All_Annotation.table', 'w') as outfile:
                    for fname in filenames2:
                        with open(fname) as infile:
                            for line in infile:
                                outfile.write(line)
                filenames = fileOver
                with open("./" + data2+'shuffle'+str(x+1) + "/All_" + data2+'shuffle'+str(x+1) + ".csv", 'w') as outfile:
                    with open(filenames) as infile:
                        for line in infile:
                            outfile.write(line)

        elif fileUnder != "None":
            filenames2 = ['./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + "_NotExpressed" + '_Annotation.table']
            with open('./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + '_All_Annotation.table', 'w') as outfile:
                for fname in filenames2:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            filenames = fileUnder
            with open("./" + data2+'shuffle'+str(x+1) + "/All_" + data2+'shuffle'+str(x+1) + ".csv", 'w') as outfile:
                with open(filenames) as infile:
                    for line in infile:
                        outfile.write(line)
    core.Dropa_Enrichment('./' + data2 + '/' + data2 + '_intergenic.bed', './' + data2 + '/' + data2 + "_All_Annotation.table", shuffle, data2, data2)

    for x in range(1, shuffle+1):
        shutil.rmtree(data2 + 'shuffle' + str(x), ignore_errors=True)
        os.remove(data2 + 'shuffle' + str(x) + ".bed")

my_file = Path("./" + data2 + "/" + data2 +"Genes")
if my_file.is_file():
    os.remove("./" + data2 + "/" + data2 +"Genes")

my_file = Path("./" + data2 + "/" + data2 +"PeaksinGenes")
if my_file.is_file():
    os.remove("./" + data2 + "/" + data2 +"PeaksinGenes")

my_file = Path("./" + data2 + "/" + data2 + "Intergenic.png")
if my_file.is_file():
    os.remove("./" + data2 + "/" + data2 + "Intergenic.png")

my_file = Path("./" + data2 + "/infoGenes" + data2)
if my_file.is_file():
    os.remove("./" + data2 + "/infoGenes" + data2)

my_file = Path("./" + data2 + "/infoRepeatsExpressed" + data2 + ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/infoRepeatsExpressed" + data2 + ".csv")

my_file = Path("./" + data2 + "/infoRepeatsNotExpressed" + data2 + ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/infoRepeatsNotExpressed" + data2 + ".csv")

my_file = Path("./" + data2 + "/Expressed_" + data2 + "exon")
if my_file.is_file():
    os.remove("./" + data2 + "/Expressed_" + data2 + "exon")

my_file = Path("./" + data2 + "/NotExpressed_" + data2 + "exon")
if my_file.is_file():
    os.remove("./" + data2 + "/NotExpressed_" + data2 + "exon")

my_file = Path("./" + data2 + "/SameOverlapping" + data2)
if my_file.is_file():
    os.remove("./" + data2 + "/SameOverlapping" + data2)

my_file = Path("./" + data2 + "/Expressed_" + data2 + "tableFoldChange")
if my_file.is_file():
    os.remove("./" + data2 + "/Expressed_" + data2 + "tableFoldChange")

my_file = Path("./" + data2 + "/NotExpressed_" + data2 + "tableFoldChange")
if my_file.is_file():
    os.remove("./" + data2 + "/NotExpressed_" + data2 + "tableFoldChange")

my_file = Path("./" + data2 + "/All_" + data2 + ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/All_" + data2 + ".csv")

my_file = Path("./" + data2 + "/NotExpressed_" + data2)
if my_file.is_file():
    os.remove("./" + data2 + "/NotExpressed_" + data2)

my_file = Path("./" + data2 + "/Expressed_" + data2)
if my_file.is_file():
    os.remove("./" + data2 + "/Expressed_" + data2)

my_file = Path("./" + data2 + "/NotExpressed_" + data2+ ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/NotExpressed_" + data2+ ".csv")

my_file = Path("./" + data2 + "/Expressed_" + data2+ ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/Expressed_" + data2+ ".csv")

my_file = Path("./" + data2 + "/Expressed_" + data2+ "Intergenic.pdf")
if my_file.is_file():
    os.remove("./" + data2 + "/Expressed_" + data2+ "Intergenic.pdf")

my_file = Path("./" + data2 + "/NotExpressed_" + data2+ "Intergenic.pdf")
if my_file.is_file():
    os.remove("./" + data2 + "/NotExpressed_" + data2+ "Intergenic.pdf")
