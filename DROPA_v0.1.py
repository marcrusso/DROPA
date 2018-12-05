#! /usr/bin/Rscript


from pathlib import Path
import argparse
import os
import DROPAcore as ex
import shutil

############ INPUT OF VARIABLES ###################

parser = argparse.ArgumentParser(description='Input program')
parser.add_argument('namefile', type=str, help='bed file with peaks')
parser.add_argument('-a', type=str, help='folder with gen annotation files')
parser.add_argument('-o',type=str,help= 'Direction of the output')
parser.add_argument('-rnaseq',type=str, help='Rna-seq file, default is isoforms.fpkm_tracking. None for no FPKM comparison')
parser.add_argument('-li',type=float, help='FPKM limit, default is 0.5')
parser.add_argument('-dis',type=int, help='distance from TSS and TTS, default is 5000 bases')
parser.add_argument('-shuffles',type=int , help='Numbers for a shuffle comparison, 0 for none. Default is 0')
parser.set_defaults(rnaseq='isoforms.fpkm_tracking',li=0.5,dis=5000,shuffle=False)

data=parser.parse_args().namefile
data2=parser.parse_args().o
rnaSeq=parser.parse_args().rnaseq
distance=parser.parse_args().dis
limit=parser.parse_args().li
shuffle=parser.parse_args().shuffles
annotation = parser.parse_args().a


################### BODY OF THE PROGRAM ###################
Utr3 = str(str(annotation) + '/'+'3UTR.bed')
Utr5 = str(str(annotation) + '/'+'5UTR.bed')
CDE = str(str(annotation) + '/'+'CDE.bed')
annotation = str('./' + str(annotation) + '/'+'anotationgen.bed')
if shuffle == None:
    shuffle = 0

gainresults=ex.PeakOverlap(annotation,data,tssdistance=distance,peakname=data2)

checkresults=ex.CheckExpression(gainresults[0],rnaSeq,limit=limit,peakname=data2,TSSTTSdistance=distance)

command = 'Rscript'

path2script = './Rplot.R'

try:
    FeatureAssignOverResults=ex.FeatureAssign(checkresults[0],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2,lap='Expressed_')
    ExonOverResults=ex.TableCreator(FeatureAssignOverResults,namepeak=data2,lap='Expressed')
    # Variable number of args in a list
    args = [data2,'Expressed_'+data2]
    # Build subprocess command
    cmd = [command, path2script] + args
    cmd = ' '.join(cmd)
    # check_output will run the command and store to result
    os.system(cmd)
    fileOver = "./" + data2 + "/Expressed_" + data2 +".csv"
except:
    print('There is no peaks in express genes')
    fileOver = "None"

try:
    FeatureAssignUnderResults=ex.FeatureAssign(checkresults[1],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2,lap='NotExpressed_')
    ExonUnderResults=ex.TableCreator(FeatureAssignUnderResults,namepeak=data2,lap='NotExpressed')
    args = [data2,'NotExpressed_'+data2]
    cmd = [command, path2script] + args
    cmd = ' '.join(cmd)
    os.system(cmd)
    fileUnder = "./" + data2 + "/NotExpressed_" + data2 +".csv"
except:
    print('There is no peaks in unexpress genes')
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
        command = 'Rscript'

        path2script = './All_data_plot.R'
        args = [data2, 'All_' + data2]
        cmd = [command, path2script] + args
        cmd = ' '.join(cmd)
        os.system(cmd)
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
        command = 'Rscript'

        path2script = './All_data_plot.R'
        args = [data2, 'All_' + data2]
        cmd = [command, path2script] + args
        cmd = ' '.join(cmd)
        os.system(cmd)
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
    command = 'Rscript'

    path2script = './All_data_plot.R'
    args = [data2, 'All_' + data2]
    cmd = [command, path2script] + args
    cmd = ' '.join(cmd)
    os.system(cmd)

############### MAKE THE COMPARATION WITH SHUFFLE #######################
args = [data2, 'Expressed_' + data2] # Reset args for R
args2 = [data2, 'NotExpressed_' + data2] # Args for the under expression
if shuffle != 0:
    for x in range(0, shuffle):
        args = args + [data2 + 'shuffle' + str(x + 1), 'Expressed_' + data2 + 'shuffle' + str(x + 1)]
        args2 = args2 + [data2 + 'shuffle' + str(x + 1), 'NotExpressed_' + data2 + 'shuffle' + str(x + 1)]
        shufflepeaks=ex.randpeak(data,x+1,data2)

        resultsgain=ex.PeakOverlap(annotation,shufflepeaks,tssdistance=distance,peakname=data2+'shuffle'+str(x+1))

        if os.stat(resultsgain[0]).st_size == 0:
            print('Shuffle1 has no Intergenic Peaks')
        else:
            resultscheck=ex.CheckExpression(resultsgain[0],rnaSeq,limit=limit,peakname=data2+'shuffle'+str(x+1),TSSTTSdistance=5000)

        try:
            resultsFeatureAssignOver=ex.FeatureAssign(resultscheck[0],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2+'shuffle'+str(x+1),lap='Expressed')

            finalResultsOver=ex.TableCreator(resultsFeatureAssignOver,namepeak=data2+'shuffle'+str(x+1),lap='Expressed')

        except:
            print('There is no peaks in express genes in Shuffle '+ str(x+1))

        try:
            resultsFeatureAssignUnder=ex.FeatureAssign(resultscheck[1],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2+'shuffle'+str(x+1),lap='NotExpressed')

            finalResultsUnder=ex.TableCreator(resultsFeatureAssignUnder,namepeak=data2+'shuffle'+str(x+1),lap='NotExpressed')

        except:
            print('There is no peaks in unexpress genes in Shuffle ' + str(x+1))


    command = 'Rscript'

    path2script = './PlotHisto2_0.R'
    cmd = [command, path2script] + args
    cmd = ' '.join(cmd)
    os.system(cmd)

    cmd = [command, path2script] + args2
    cmd = ' '.join(cmd)
    os.system(cmd)

    path2script = './PlotAllHisto2_0.R'

    args = [data2, 'Expressed_' + data2,'NotExpressed_' + data2]
    cmd = [command, path2script] + args
    cmd = ' '.join(cmd)
    os.system(cmd)

    for x in range(0, shuffle):
        shutil.rmtree(data2 + 'shuffle' + str(x + 1), ignore_errors=True)
        os.remove(data2 + 'shuffle' + str(x + 1) + ".bed")

my_file = Path("./" + data2 + "/" + data2 +"Genes")
if my_file.is_file():
    os.remove("./" + data2 + "/" + data2 +"Genes")

my_file = Path("./" + data2 + "/" + data2 +"PeaksInGenes")
if my_file.is_file():
    os.remove("./" + data2 + "/" + data2 +"PeaksInGenes")

my_file = Path("./" + data2 + "/" + data2 + "Intergenic.png")
if my_file.is_file():
    os.remove("./" + data2 + "/" + data2 + "Intergenic.png")

my_file = Path("./" + data2 + "/infoGenes" + data2)
if my_file.is_file():
    os.remove("./" + data2 + "/infoGenes" + data2)

my_file = Path("./" + data2 + "/infoRepeatsExpressed_" + data2 + ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/infoRepeatsExpressed_" + data2 + ".csv")

my_file = Path("./" + data2 + "/infoRepeatsNotExpressed_" + data2 + ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/infoRepeatsNotExpressed_" + data2 + ".csv")

my_file = Path("./" + data2 + "/infoRepeatsNotExpressed_" + data2 + ".csv")
if my_file.is_file():
    os.remove("./" + data2 + "/infoRepeatsExpressed_" + data2 + ".csv")

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



