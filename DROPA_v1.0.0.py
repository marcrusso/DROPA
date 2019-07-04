#! /usr/bin/Rscript


from pathlib import Path
import argparse
import os
import DROPAcore as core
import warnings
import shutil
warnings.filterwarnings("ignore")

################### INPUT OF VARIABLES ###################

parser = argparse.ArgumentParser(description='Input program')
parser.add_argument('namefile', type=str, help='Query BED file')
parser.add_argument('-ref', type=str, help='Gene Reference Folder')
parser.add_argument('-o',type=str,help= 'Output directory')
parser.add_argument('-ex',type=str, help='Expression file. Default is "None"')
parser.add_argument('-lim',type=float, help='Expression Threshold. Default is 0.5')
parser.add_argument('-dis', help='Distance from TSS and TTS, to define Upstream and Downstream regions. Default is 5000')
parser.add_argument('-shuffle',type=int , help='Numbers for a shuffle comparison, 0 for none. Default is 0')
parser.add_argument('-gsize',type=str , help='Genome size file, used for shuffle')
parser.set_defaults(ex=None,lim=0.5,dis=[5000,5000],shuffle=0)

data=parser.parse_args().namefile
data2=parser.parse_args().o
rnaSeq=parser.parse_args().ex
distance=parser.parse_args().dis
limit=parser.parse_args().lim
shuffle=parser.parse_args().shuffle
annotation = parser.parse_args().ref
chromosome_size_file = parser.parse_args().gsize


################### BODY OF THE PROGRAM ###################

if ',' in distance:
    distance = distance.split(',')
    distance = [int(distance[0]),int(distance[1])]
elif type(distance) == str:
    distance = [int(distance), int(distance)]

Utr3 = str(str(annotation) + '/'+'3UTR.bed')
Utr5 = str(str(annotation) + '/'+'5UTR.bed')
CDE = str(str(annotation) + '/'+'CDE.bed')
annotation = str('./' + str(annotation) + '/'+'refgene.bed')

gainresults = core.PeakOverlap(annotation,data,tssdistance=distance,peakname=data2)


checkresults = core.CheckExpression(gainresults[0],rnaSeq,limit=limit,peakname=data2,TSSTTSdistance=distance)

try:
    FeatureAssignOverResults=core.FeatureAssign(checkresults[0],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2,lap='Expressed_')
    ExonOverResults= core.TableCreator(FeatureAssignOverResults,namepeak=data2,lap='Expressed', peaksingenes= gainresults[0], warnings= True)
    core.Upsetplotting(ExonOverResults,"Expressed_"+data2, data2)
    core.Dropa_pie_AllwithIntergenic(gainresults[1],ExonOverResults, "Expressed_" + data2, data2, False)
except:
    print('There are no peaks overlapping expressed genes')
    ExonOverResults = "None"

try:
    FeatureAssignUnderResults=core.FeatureAssign(checkresults[1],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2,lap='NotExpressed_')
    ExonUnderResults=core.TableCreator(FeatureAssignUnderResults,namepeak=data2,lap='NotExpressed', peaksingenes= gainresults[0], warnings= True)
    core.Upsetplotting(ExonUnderResults, 'NotExpressed_'+data2, data2)
    core.Dropa_pie_AllwithIntergenic(gainresults[1],ExonUnderResults, "NotExpressed_" + data2, data2, False)
except:
    print('There are no peaks overlapping unexpressed genes')
    ExonUnderResults = "None"

#### Make a plot of the total percent of each zone been overlap by a peak for express and not express genes together
h,k=1,1
if ExonOverResults != "None":
    if ExonUnderResults != "None":
        filenames = [ExonOverResults, ExonUnderResults]
        with open('./' + data2 + '/' + data2 + '_All_Annotation.table', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        if h==1:
                            outfile.write(line)
                        else:
                            h=1
                    h=0
    else:
        filenames = ExonOverResults
        with open('./' + data2 + '/' + data2 + '_All_Annotation.table', 'w') as outfile:
            with open(filenames) as infile:
                for line in infile:
                    outfile.write(line)

elif ExonUnderResults != "None":
    filenames = ExonUnderResults
    with open('./' + data2 + '/' + data2 + '_All_Annotation.table', 'w') as outfile:
        with open(filenames) as infile:
            for line in infile:
                outfile.write(line)

core.Upsetplotting('./' + data2 + '/' + data2 + '_All_Annotation.table',"All"+data2, data2)

core.Dropa_pie_Intergenic('./' + data2 + '/' + data2 + '_intergenic.bed','./' + data2 + '/' + data2 + 'PeaksInGenes',data2, data2)

core.Dropa_pie_AllwithIntergenic('./' + data2 + '/' + data2 + '_intergenic.bed','./' + data2 + '/' + data2 + "_All_Annotation.table","All_"+data2, data2, True)

core.Dropa_histogram('./' + data2 + '/' + data2 + '_intergenic.bed','./' + data2 + '/' + data2 + "_All_Annotation.table", data, "All_"+data2, data2)

################### COMPARISON WITH SHUFFLE ###################

if shuffle != 0:
    for x in range(0, shuffle):
        finalResultsOver, finalResultsUnder = None, None
        shufflepeaks=core.randpeak(data,x+1,data2,chromosome_size_file)

        resultsgain=core.PeakOverlap(annotation,shufflepeaks,tssdistance=distance,peakname=data2+'shuffle'+str(x+1))

        if os.stat(resultsgain[0]).st_size == 0:
            print('Shuffle1 has no Intergenic Peaks')
        else:
            resultscheck=core.CheckExpression(resultsgain[0],rnaSeq,limit=limit,peakname=data2+'shuffle'+str(x+1),TSSTTSdistance=distance)

        try:
            resultsFeatureAssignOver=core.FeatureAssign(resultscheck[0],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2+'shuffle'+str(x+1),lap='Expressed')

            finalResultsOver=core.TableCreator(resultsFeatureAssignOver,namepeak=data2+'shuffle'+str(x+1),lap='Expressed', peaksingenes= resultsgain[0])

        except:
            print('There are no peaks in expressed genes in Shuffle '+ str(x+1))

        try:
            resultsFeatureAssignUnder=core.FeatureAssign(resultscheck[1],UTR3=Utr3,UTR5=Utr5,CDE=CDE,peakname=data2+'shuffle'+str(x+1),lap='NotExpressed')

            finalResultsUnder=core.TableCreator(resultsFeatureAssignUnder,namepeak=data2+'shuffle'+str(x+1),lap='NotExpressed', peaksingenes=resultsgain[0])

        except:
            print('There are no peaks in unexpressed genes in Shuffle ' + str(x+1))
        h, k = 1, 1
        if finalResultsOver != None:
            if finalResultsUnder != None:
                filenames = [finalResultsOver, finalResultsUnder]
                with open('./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + '_All_Annotation.table', 'w') as outfile:
                    for fname in filenames:
                        with open(fname) as infile:
                            for line in infile:
                                if h == 1:
                                    outfile.write(line)
                                else:
                                    h = 1
                            h = 0

            else:
                filenames = finalResultsOver
                with open('./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + '_All_Annotation.table', 'w') as outfile:
                    with open(filenames) as infile:
                        for line in infile:
                            outfile.write(line)

        elif finalResultsUnder != None:
            filenames = finalResultsUnder
            with open('./' + data2+'shuffle'+str(x+1) + '/' + data2+'shuffle'+str(x+1) + '_All_Annotation.table', 'w') as outfile:
                with open(filenames) as infile:
                    for line in infile:
                        outfile.write(line)

    core.Dropa_Enrichment('./' + data2 + '/' + data2 + '_intergenic.bed', './' + data2 + '/' + data2 + "_All_Annotation.table", shuffle, data2, data2)

################### Erase shuffle folders and excess files ###################

    for x in range(1, shuffle+1):
        shutil.rmtree(data2 + 'shuffle' + str(x), ignore_errors=True)
        os.remove(data2 + 'shuffle' + str(x) + ".bed")

my_file = Path("./" + data2 + "/" + data2 +"PeaksinGenes")
if my_file.is_file():
    os.remove(my_file)

my_file = Path("./" + data2 + "/Expressed_" + data2 + "exon")
if my_file.is_file():
    os.remove(my_file)

my_file = Path("./" + data2 + "/NotExpressed_" + data2 + "exon")
if my_file.is_file():
    os.remove(my_file)

my_file = Path("./" + data2 + "/All_" + data2 + ".csv")
if my_file.is_file():
    os.remove(my_file)

my_file = Path("./" + data2 + "/NotExpressed_" + data2)
if my_file.is_file():
    os.remove(my_file)

my_file = Path("./" + data2 + "/Expressed_" + data2)
if my_file.is_file():
    os.remove(my_file)

################### -~- ###################