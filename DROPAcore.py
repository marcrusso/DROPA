
# coding: utf8
import os
import numpy as np
from tqdm import tqdm
import pandas as pd
import time
from collections import Counter
import csv
### INTERVALTREE ONLY WORKS ON PYTHON 3.5 OR MINOR VERSION

from intervaltree import IntervalTree

def PeakOverlap(genesfile, peaksfile,tssdistance=0,peakname='null'):
    LuckPeak, LuckGen, LuckTree, LuckBegin , Genlist = {},{},{},{},{}

####### CREATE A INTERVALTREE VARIABLE

    tree = IntervalTree()
    n = 0
    m = 0
    intergenic = set()
    intergenic_output = {}
    for lines in open(peaksfile):
        fields = lines.split()
        namegain, chromogain, begingain, endgain = fields[3], fields[0], int(fields[1]), int(fields[2])
        space4, space5 = fields[4], fields[5]
        LuckPeak[namegain] = [chromogain, begingain, endgain, namegain, space4, space5]
        LuckBegin[begingain] = [namegain,begingain,endgain]
        intergenic = intergenic|set([namegain])

        if chromogain not in LuckTree:
            print('Chromosome '+chromogain+' of ' +peakname+'...')
            LuckTree[chromogain] = 0
            if n == 1:
                for lines2 in open(genesfile):
                    fields2 = lines2.split()
                    if fields2[0] != k:
                        continue
                    else:
                        nameid = fields2[3]
                        begingen = int(fields2[1]) - tssdistance
                        endgen = int(fields2[2]) + tssdistance
                        chromogen = fields2[0]
                        strand = fields2[5]
                        if tree.overlap(begingen, endgen) != set():
                            Genlist[nameid] = [chromogen] + [fields2[1]] + [fields2[2]] + [nameid] + [strand]
                            for x in tree.overlap(begingen, endgen):
                                LuckGen[m] = [chromogen] + [fields2[1]] + [fields2[2]] + [nameid] + [strand] + LuckBegin[x.begin]
                                intergenic = intergenic - set([LuckBegin[x.begin][0]])
                                m+=1
            else:
                tree[begingain:endgain] = (begingain, endgain)
            n = 1
            ### RESET THE TREE EACH TIME BEFORE START A NEW CHROMOSOME
            tree = IntervalTree()
            tree[begingain:endgain] = (begingain, endgain)
        ### get all the peaks of the chromosome to fill the tree until the next item of the field is another chromosome. Then start to compare all items of the tree with all the genes int he same chromosome
        else:
            k = chromogain
            tree[begingain:endgain] = (begingain,endgain)
    for lines2 in open(genesfile):
        fields2 = lines2.split()
        if fields2[0] != k:
            continue
        else:
            nameid = fields2[3]
            begingen = int(fields2[1]) - tssdistance
            endgen = int(fields2[2]) + tssdistance
            chromogen = fields2[0]
            strand = fields2[5]
            if tree.overlap(begingen, endgen) != set():
                Genlist[nameid] = [chromogen] + [fields2[1]] + [fields2[2]] + [nameid] + [strand]
                for x in tree.overlap(begingen, endgen):
                    LuckGen[m] = [chromogen] + [fields2[1]] + [fields2[2]] + [nameid] + [strand] + LuckBegin[x.begin]
                    intergenic = intergenic - set([LuckBegin[x.begin][0]])
                    m += 1


    for x in intergenic:
        intergenic_output[x] = LuckPeak[x]

    ### OUTPUT
    if not os.path.exists(peakname):
        os.makedirs(peakname)

    if len(intergenic) == 0:
        print('No Intergenic peaks')
    else:
        results_intragenic = pd.DataFrame(list(intergenic_output.values())).sort_values(by=[0])
        results_intragenic.to_csv('./' + peakname + '/' + peakname + '_intergenic.bed', index=None, sep='\t', header=False)

    results = pd.DataFrame(list(LuckGen.values()))
    results_genes = pd.DataFrame(list(Genlist.values())).sort_values(by=0)
    results.to_csv('./' + peakname + '/' + peakname + 'PeaksInGenes', index=None, sep='\t', header= False)
    results_genes.to_csv('./' + peakname + '/' + peakname + 'Genes', index=None, sep='\t', header=False)
    return ('./' + peakname + '/' + peakname + 'PeaksInGenes', './' + peakname + '/' + peakname + 'Genes',
            './' + peakname + '/' + peakname + '_intergenic.bed')


def CheckExpression(peakgenes, transcribe, limit=0.5,peakname='null',TSSTTSdistance=0):
    # LOAD THE TRANSCRIPT/GENES DATA
    tic = time.clock()
    transcript_info = {}
    count_lines = 0
    print('Checking Transcripts/genes overlapped by peaks and selecting the best match')
    if transcribe != 'None':
        if os.stat(transcribe).st_size == 0:
            return print('RNA-seq file is empty')
        for line in open(transcribe, "r"):
            fields = line.split()
            if len(fields) != 2:
                print(fields[0])
                return print('Please input a correct RNA-seq file')
            try:
                FPKM = float(fields[1])
                pass
            except:
                count_lines +=1
                if count_lines != 1:
                    return print('Please input a correct RNA-seq file')
                continue
            trascr_id = fields[0]

            transcript_info[trascr_id] = [FPKM]
    else:
        for line in open(peakgenes, "r"):
            fields = line.split()
            FPKM = float(0)
            trascr_id = fields[3]

            transcript_info[trascr_id] = [FPKM]

    # LOAD THE PEAK-IN-GENE DATA
    transcr_id_not_found = 0
    peaks_to_transcripts = {}
    total = 0
    for line in open(peakgenes, "r"):
        total += 1
        chrom, transcr_start, transcr_end, transcr_id, strand, peak_name, peak_start, peak_end = line.split()
        if peak_name not in peaks_to_transcripts:
            peaks_to_transcripts[peak_name] = []

        if transcr_id not in transcript_info:
            transcr_id_not_found += 1
            continue

        transcript_data = transcript_info[transcr_id]
        peaks_to_transcripts[peak_name].append(
            transcript_data + [chrom, transcr_start, transcr_end, transcr_id, strand, peak_start, peak_end])
    print('Files Loaded')

###### SELECT TRANSCRIPTS OR GENES BY EXPRESION. ALL MAX IN THE BEST_TRANSCRIPT, IF THERE IS NO MAX (FPKM = 0): SAVE IT IN NOMAX, MORE THAN 1 MAX: SAVE IN PEAK2MAX.
    best_transcript_under = {}
    best_transcript_over = {}
    nobest_transcript = {}
    over_over = 0
    over_only = 0
    under_over = 0
    under_only = 0
    for peak in peaks_to_transcripts:
        transcripts = peaks_to_transcripts[peak]
        fpkm = [x[0] for x in transcripts]
        maxfpkm = max(fpkm)
        indices = [i for i, x in enumerate(fpkm) if x == maxfpkm]
        ### IF MAX OF FPKM IS LESS THAN THE LIMIT, USE THE QUANTITTY OF OVERLAPPING BASES AS A SELECTIVE PROCEDURE
        if maxfpkm < limit and len(fpkm) > 1:
            under_over+=1
            overlapbases = []
            for x in transcripts:
                beginpeak = int(x[6])
                endpeak = int(x[7])
                begingen = int(x[2]) - TSSTTSdistance
                endgen = int(x[3]) + TSSTTSdistance
                if beginpeak in range(begingen, endgen) and endpeak not in range(begingen, endgen):
                    overlapbases.append(endgen - beginpeak)
                elif beginpeak in range(begingen, endgen) and endpeak in range(begingen, endgen):
                    overlapbases.append(endpeak - beginpeak)
                elif beginpeak not in range(begingen, endgen) and endpeak in range(begingen, endgen):
                    overlapbases.append(endpeak - begingen)
                else:
                    overlapbases.append(endgen - begingen)
            best_transcript_under[peak] = transcripts[overlapbases.index(max(overlapbases))]
            indices_max = [j for j, m in enumerate(overlapbases) if m == max(overlapbases)]
            if len(indices_max) > 1:
                for y in indices_max:
                    if peak not in nobest_transcript:
                        nobest_transcript[peak] = [transcripts[y][5]]
                    else:
                        nobest_transcript[peak] = nobest_transcript[peak] + [transcripts[y][5]]
        ### SEE IF THERE IS MORE THAN 1 MAX AND SELECT THAT IT HAS MORE OVERLAPPING BASES
        elif len(indices) > 1:
            overlapbases = []
            aux_trans = [0 for x in range(len(indices))]
            n=0
            for x in indices:
                aux_trans[n] = transcripts[x]
                beginpeak = int(aux_trans[n][6])
                endpeak = int(aux_trans[n][7])
                begingen = int(aux_trans[n][2]) - TSSTTSdistance
                endgen = int(aux_trans[n][3]) + TSSTTSdistance
                if beginpeak in range(begingen, endgen) and endpeak not in range(begingen, endgen):
                    overlapbases.append(endgen - beginpeak)
                elif beginpeak in range(begingen, endgen) and endpeak in range(begingen, endgen):
                    overlapbases.append(endpeak - beginpeak)
                elif beginpeak not in range(begingen, endgen) and endpeak in range(begingen, endgen):
                    overlapbases.append(endpeak - begingen)
                else:
                    overlapbases.append(endgen - begingen)
                n+=1
            if max(overlapbases) > limit:
                over_over += 1
                best_transcript_over[peak] = aux_trans[overlapbases.index(max(overlapbases))]
            else:
                best_transcript_under[peak] = aux_trans[overlapbases.index(max(overlapbases))]
                under_over += 1
        ### THIS IS PART I USE THE ALGORITHM
        else:
            max_trans = max(transcripts, key=lambda x: x[0])

            if max_trans[0] > limit:
                best_transcript_over[peak] = max_trans
                if len(fpkm) > 1:
                    over_over += 1
                else:
                    over_only += 1
            else:
                best_transcript_under[peak] = max_trans
                if len(fpkm) > 1:
                    under_over += 1
                else:
                    under_only += 1
    toc = time.clock()
    print('Done in ' + str(round(toc - tic, 2)) + ' sec. Printing the results of the checking')


    ### ORDER THE RESULTS FOR OUTPUT AND SAVE THEM IN A FILE
    if len(nobest_transcript) > 0:
        print('There are '+ str(len(nobest_transcript)) + ' peaks overlapping 2 or more genes/transcript with the same quantity of overlapping bases and no difference in expression level')
        w = csv.writer(open('./' + peakname + '/SameOverlapping' + peakname, "w"))
        for key, val in nobest_transcript.items():
            w.writerow([key, val])

    if len(best_transcript_under) > 0:
        best_transcript_under = pd.DataFrame(best_transcript_under).transpose()
        best_transcript_under['Peak_Name'] = best_transcript_under.index
        best_transcript_under.columns=['FPKM','Chromo', 'Gen_Start','Gen_End','Gen_TransID','Strand','Peak_Start','Peak_End','Peak_Name']
        best_transcript_under = best_transcript_under[['Chromo','Gen_Start','Gen_End','Gen_TransID','Strand','Peak_Name','Peak_Start','Peak_End','FPKM']]
        best_transcript_under.to_csv('./' + peakname + '/NotExpressed_' + peakname, index=None, sep='\t')

    else:
        print('No express genes with overlapping peaks')

    if len(best_transcript_over) > 0:
        best_transcript_over = pd.DataFrame(best_transcript_over).transpose()
        best_transcript_over['Peak_Name'] = best_transcript_over.index
        best_transcript_over.columns = ['FPKM', 'Chromo', 'Gen_Start', 'Gen_End', 'Gen_TransID', 'Strand',
                                         'Peak_Start', 'Peak_End', 'Peak_Name']
        best_transcript_over = best_transcript_over[
            ['Chromo', 'Gen_Start', 'Gen_End', 'Gen_TransID', 'Strand', 'Peak_Name', 'Peak_Start', 'Peak_End', 'FPKM']]
        best_transcript_over.to_csv('./' + peakname + '/Expressed_' + peakname, index=None, sep='\t')
    else:
        print('No Unexpress genes with overlapping peaks')

    print(' Number of peaks: ' + str(len(peaks_to_transcripts)))
    print(' There are ' + str(over_over) + ' peaks with overlapping genes and over the FPKM limit \n There are ' + str(under_over) + ' peaks with overlapping genes and under the FPKM limit'
              + '\n There are ' + str(
          over_only) + ' peaks with no overlapping genes and over the FPKM limit \n There are ' + str(
          under_only) + ' peaks with no overlapping genes and under the FPKM limit')
    results = open('./' + peakname + '/infoGenes' + peakname, 'w')
    results.write(' Number of peaks: ' + str(len(peaks_to_transcripts)) +'\n'+' There are '+str(over_over)+' peaks with overlapping genes and over the FPKM limit \n There are '+str(under_over)+ ' peaks with overlapping genes and under the FPKM limit'
        +'\n There are '+str(over_only)+' peaks with no overlapping genes and over the FPKM limit \n There are '+str(under_only)+' peaks with no overlapping genes and under the FPKM limit')

    print('total peaks overlapping transcripts/genes: ' + str(over_over+over_only+under_over+under_only))

    return ('./' + peakname + '/Expressed_' + peakname, './' + peakname + '/NotExpressed_' + peakname)


def TableCreator(peakexon,namepeak='null',lap='null'):
    exonchromo,exonbegin, exonend, exonname,exondirection,exontype,peakname,peakchromo,peakbegin, peakend,begingen, endgen=[[],[],[],[],[],[],[],[],[],[],[],[]]
    score=[]
    data=open(peakexon)
    a=1
    for lines in data.readlines():
        if a>1:
            exonchromo.append(lines.split('\t')[0])
            exonbegin.append(int(lines.split('\t')[9]))
            exonend.append(int(lines.split('\t')[10]))
            exonname.append(lines.split('\t')[3])
            peakend.append(int(lines.split('\t')[7]))
            peakbegin.append(int(lines.split('\t')[6]))
            exondirection.append(lines.split('\t')[4])
            exontype.append(lines.split('\t')[11].replace('\n',''))
            peakname.append(lines.split('\t')[5])
            score.append(lines.split('\t')[8])
            peakchromo.append(lines.split('\t')[0])
            begingen.append(int(lines.split('\t')[1]))
            endgen.append(int(lines.split('\t')[2]))
        a+=1
    data.close()

    begingen,endgen=[np.array(begingen),np.array(endgen)]
    exonbegin, exonend, peakbegin, peakend, exondirection = [np.array(exonbegin),np.array(exonend),np.array(peakbegin),np.array(peakend),np.array(exondirection)]
    countpeakexons, countpeakintrons, countpeak5UTR, countpeak3UTR,countpeakTSS, countpeakTTS =[[],[],[],[],[],[]]
    score=np.array(score)
    matrixlong=-1
    resnamegen, reschromo, resdirection,respeakstart,respeakends,respeakname =[[],[],[],[],[],[]]
    exonname,exontype = np.array(exonname),np.array(exontype)
    peakname, exonchromo = [np.array(peakname), np.array(exonchromo)]
    resscore=[]
    repeattimes=[]
    uniquepeak=np.unique(peakname)
    for peak in tqdm(range(len(uniquepeak)), desc='TableCreator'):
        indices = np.where(peakname==uniquepeak[peak])
        countrep = -1
        reschromo.append(exonchromo[indices][0])
        resnamegen.append(exonname[indices][0])
        resdirection.append(exondirection[indices][0])
        respeakstart.append(peakbegin[indices][0])
        respeakends.append(peakend[indices][0])
        respeakname.append(uniquepeak[peak])
        resscore.append(score[indices][0])
        matrixlong += 1
        countpeak5UTR.append(0)
        countpeak3UTR.append(0)
        countpeakexons.append(0)
        countpeakintrons.append(0)
        countpeakTSS.append(0)
        countpeakTTS.append(0)
        for indice in range(len(exonname[indices])):
            if ((peakbegin[indices][indice] <= begingen[indices][indice] and exondirection[indices][indice] == '+') or (
                    peakend[indices][indice] >= endgen[indices][indice] and exondirection[indices][indice] == '-'))and countpeakTSS[matrixlong] != 1:
                countpeakTSS[matrixlong] = 1
                countrep += 1
            if ((peakbegin[indices][indice] <= begingen[indices][indice] and exondirection[indices][indice] == '-') or (
                    peakend[indices][indice] >= endgen[indices][indice] and exondirection[indices][indice] == '+')) and countpeakTTS[matrixlong] != 1:
                countpeakTTS[matrixlong] = 1
                countrep += 1

            if (exonbegin[indices][indice] <= peakbegin[indices][indice] <= exonend[indices][indice] or exonbegin[indices][indice] <= peakend[indices][indice] <= exonend[indices][indice] or (
                    exonbegin[indices][indice] >= peakbegin[indices][indice] and peakend[indices][indice] >= exonend[indices][indice])):
                if exontype[indices][indice]=='5UTR' and countpeak5UTR[matrixlong] != 1:
                    countpeak5UTR[matrixlong] = 1
                    countrep += 1
                elif exontype[indices][indice]=='3UTR' and countpeak3UTR[matrixlong] != 1:
                    countpeak3UTR[matrixlong] = 1
                    countrep += 1
                elif countpeakexons[matrixlong] != 1:
                    countpeakexons[matrixlong] = 1
                    countrep += 1
            try:
                if ((exonend[indices][indice] <= peakbegin[indices][indice] <= exonbegin[indices][indice+1] or exonend[indices][indice] <= peakend[indices][indice] <= exonbegin[indices][indice+1] or (exonend[indices][indice] >= peakbegin[indices][indice] and peakend[indices][indice] >= exonbegin[indices][indice]))) and countpeakintrons[matrixlong] != 1:
                    countpeakintrons[matrixlong] = 1
                    countrep += 1
            except:
                pass
        repeattimes.append(countrep)
    repeatcounter = Counter(repeattimes)
    if not os.path.exists(namepeak):
        os.makedirs(namepeak)
    results=open('./'+namepeak+'/'+namepeak+"_"+lap+'_Annotation.table','w')
    results.write('Chromo'+'\t'+'Peak_Start'+'\t'+'Peak_End'+'\t'+'Peak_Name'+'\t'+'ID_transc'+'\t'+'FPKM'+'\t'+'Direction'+'\t'+'TSS'+'\t'+'5UTR'+'\t'+'Exons'+'\t'+'Introns'+'\t'+'3UTR'+'\t'+'TTS'+'\t'+'\n')
    for res in range(len(resnamegen)):
        results.write(reschromo[res]+'\t'+str(respeakstart[res])+'\t'+str(respeakends[res])+'\t'+respeakname[res]+'\t'+resnamegen[res]+'\t'+str(resscore[res])+'\t'+resdirection[res]+'\t'+str(countpeakTSS[res])+'\t'+str(countpeak5UTR[res])+'\t'+str(countpeakexons[res])+'\t'+str(countpeakintrons[res])+'\t'+str(countpeak3UTR[res])+'\t'+str(countpeakTTS[res])+'\t'+'\n')
    results.close()
    results=open('./'+namepeak+'/'+lap+"_"+namepeak+'.csv','w')
    results.write('\"'+'Peak_Name'+'\t'+'Upstream'+'\t'+'UTR5'+'\t'+'Exons'+'\t'+'Introns'+'\t'+'UTR3'+'\t'+'Downstream'+'\t'+'\"'+','+','+'\n')
    for res in range(len(resnamegen)):
        results.write('\"'+respeakname[res]+'\t'+str(countpeakTSS[res])+'\t'+str(countpeak5UTR[res])+'\t'+str(countpeakexons[res])+'\t'+str(countpeakintrons[res])+'\t'+str(countpeak3UTR[res])+'\t'+str(countpeakTTS[res])+'\t'+'\"'+','+','+'\n')
    results.close()
    w = csv.writer(open('./'+namepeak+'/infoRepeats'+lap+namepeak+'.csv', "w"))
    for key, val in repeatcounter.items():
        w.writerow([key, val])

    return('./'+namepeak+'/'+lap+namepeak+'.csv','./'+namepeak+'/'+namepeak+"_"+lap+'_Annotation.table')

# This program take the list of exons that are in the genes and put them in a new list where the each row is a exon with the name of the file
# Need a bed file and for the 2 inputs, been one the list of the genes with the exons size and position and the other the clean list of gain or lost genes.

def FeatureAssign(cleanpeaklist,UTR5='5UTR.bed',UTR3='3UTR.bed',CDE='CDE.bed',peakname='null',lap='null'):
    print('Creating the list of Exons for Analisis')
    tic=time.clock()
    peaktable=pd.read_csv(cleanpeaklist, sep="\t")
    #   Getting the 5UTR data and merging with the table of peaks in genes
    UTR5table = pd.read_csv(UTR5, header=None, sep="\t")
    UTR5table = UTR5table.rename(columns={0:'Chromo',1: 'BeginExon',2: 'EndExon',3: 'Gen_TransID'})
    del UTR5table[4],UTR5table[5],UTR5table['Chromo']
    mergedUTR5peak = pd.merge(peaktable, UTR5table, left_on='Gen_TransID', right_on='Gen_TransID', how='inner')
    mergedUTR5peak['Type'] = pd.Series('5UTR', index=mergedUTR5peak.index)

    #   Getting the 3UTR data and merging with the table of peaks in genes
    UTR3table = pd.read_csv(UTR3, header=None, sep="\t")
    UTR3table = UTR3table.rename(columns={0:'Chromo',1: 'BeginExon',2: 'EndExon',3: 'Gen_TransID'})
    del UTR3table[4],UTR3table[5],UTR3table['Chromo']
    mergedUTR3peak = pd.merge(peaktable, UTR3table, left_on='Gen_TransID', right_on='Gen_TransID', how='inner')
    mergedUTR3peak['Type'] = pd.Series('3UTR', index=mergedUTR3peak.index)

    #   Getting the 3UTR data and merging with the table of peaks in genes
    CDEtable = pd.read_csv(CDE, header=None, sep="\t")
    CDEtable = CDEtable.rename(columns={0:'Chromo',1: 'BeginExon',2: 'EndExon',3: 'Gen_TransID'})
    del CDEtable[4],CDEtable[5],CDEtable['Chromo']
    mergedCDEpeak = pd.merge(peaktable, CDEtable, left_on='Gen_TransID', right_on='Gen_TransID', how='inner')
    mergedCDEpeak['Type'] = pd.Series('CDE', index=mergedCDEpeak.index)

    # Concatenate the tables
    finaltable=pd.concat([mergedUTR5peak,mergedUTR3peak,mergedCDEpeak])
    if not os.path.exists(peakname):
        os.makedirs(peakname)
    finaltable.to_csv('./'+peakname+'/'+lap+peakname+'exon', index=None, sep='\t')
    toc=time.clock()
    print('Done in: ' + str(round(toc-tic,2))+ ' sec')

    return ('./'+peakname+'/'+lap+peakname+'exon')


def randpeak(peaks=None, numbershuffle='', Namefile='', genome_file =""):
    command = 'bedtools shuffle -chrom -i \''+peaks+'\' -g ' + genome_file + ' -noOverlapping > '+Namefile+'shuffle'+str(numbershuffle)+'.bed'
    os.system(command)
    return(Namefile+'shuffle'+str(numbershuffle)+'.bed')

############################################################################################################################