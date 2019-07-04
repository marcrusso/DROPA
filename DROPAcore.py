
# coding: utf8
import os
import numpy as np
from tqdm import tqdm
import pandas as pd
import time
from upsetplot import plot
from matplotlib import pyplot
from scipy import stats

from intervaltree import IntervalTree

def PeakOverlap(genesfile, peaksfile,tssdistance=[0,0],peakname='null'):
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
                        begingen = int(fields2[1]) - tssdistance[0]
                        endgen = int(fields2[2]) + tssdistance[1]
                        chromogen = fields2[0]
                        strand = fields2[5]
                        if tree.overlap(begingen, endgen) != set():
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
            begingen = int(fields2[1]) - tssdistance[0]
            endgen = int(fields2[2]) + tssdistance[1]
            chromogen = fields2[0]
            strand = fields2[5]
            if tree.overlap(begingen, endgen) != set():
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
    results.to_csv('./' + peakname + '/' + peakname + 'PeaksInGenes', index=None, sep='\t', header= False)
    return ('./' + peakname + '/' + peakname + 'PeaksInGenes',
            './' + peakname + '/' + peakname + '_intergenic.bed')

def CheckExpression(peakgenes, transcribe, limit=0.5,peakname='null',TSSTTSdistance=[0,0]):

    overlap_expres_trans, peak_overlap_trans = [], []
    # LOAD THE TRANSCRIPT/GENES DATA
    tic = time.clock()
    transcript_info = {}
    count_lines = 0
    print('Checking Transcripts/genes overlapped by peaks and selecting the best match')

    if transcribe != None:
        if os.stat(transcribe).st_size == 0:
            return print('Expression file is empty')
        for line in open(transcribe, "r"):
            fields = line.split()
            if len(fields) != 2:
                print(fields[0])
                return print('Please input a correct expression file')
            try:
                FPKM = float(fields[1])
                pass
            except:
                count_lines +=1
                if count_lines != 1:
                    return print('Please input a correct expression file')
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
        peaks_to_transcripts[peak_name].append(transcript_data + [chrom, transcr_start, transcr_end, transcr_id, strand, peak_start, peak_end])

    print('Files Loaded')
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
                begingen = int(x[2]) - TSSTTSdistance[0]
                endgen = int(x[3]) + TSSTTSdistance[1]
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
        ### IF THERE IS MORE THAN 1 MAX, SELECT THE ONE THAT HAS MORE OVERLAPPING BASES
        elif len(indices) > 1:
            ## SAVE TRANSCRIPTS THAT ARE EXPRESSED (OVER FPKM)
            overlapbases = []
            aux_trans = [0 for x in range(len(indices))]
            n=0
            peak_overlap_trans.append(peak)
            overlap_expres_trans = []
            for x in indices:
                overlap_expres_trans = overlap_expres_trans + [transcripts[x]]
                aux_trans[n] = transcripts[x]
                beginpeak = int(aux_trans[n][6])
                endpeak = int(aux_trans[n][7])
                begingen = int(aux_trans[n][2]) - TSSTTSdistance[0]
                endgen = int(aux_trans[n][3]) + TSSTTSdistance[1]
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
    print('Done in ' + str(round(toc - tic, 2)) + ' sec. Printing results')


    ### ORDER THE RESULTS FOR OUTPUT AND SAVE THEM IN A FILE

    if len(best_transcript_under) > 0:
        best_transcript_under = pd.DataFrame(best_transcript_under).transpose()
        best_transcript_under['Peak_Name'] = best_transcript_under.index
        best_transcript_under.columns=['FPKM','Chromo', 'Gen_Start','Gen_End','Gen_TransID','Strand','Peak_Start','Peak_End','Peak_Name']
        best_transcript_under = best_transcript_under[['Chromo','Gen_Start','Gen_End','Gen_TransID','Strand','Peak_Name','Peak_Start','Peak_End','FPKM']]
        best_transcript_under.to_csv('./' + peakname + '/NotExpressed_' + peakname, index=None, sep='\t')

    if len(best_transcript_over) > 0:
        best_transcript_over = pd.DataFrame(best_transcript_over).transpose()
        best_transcript_over['Peak_Name'] = best_transcript_over.index
        best_transcript_over.columns = ['FPKM', 'Chromo', 'Gen_Start', 'Gen_End', 'Gen_TransID', 'Strand',
                                         'Peak_Start', 'Peak_End', 'Peak_Name']
        best_transcript_over = best_transcript_over[
            ['Chromo', 'Gen_Start', 'Gen_End', 'Gen_TransID', 'Strand', 'Peak_Name', 'Peak_Start', 'Peak_End', 'FPKM']]
        best_transcript_over.to_csv('./' + peakname + '/Expressed_' + peakname, index=None, sep='\t')

    print(' Number of peaks: ' + str(len(peaks_to_transcripts)))
    print(' There are ' + str(over_over) + ' peaks with overlapping genes and over the expression threshold \n There are ' + str(under_over) + ' peaks with overlapping genes and under the expression threshold'
              + '\n There are ' + str(
          over_only) + ' peaks with no overlapping genes and over the expression threshold \n There are ' + str(
          under_only) + ' peaks with no overlapping genes and under the expression threshold')

    print('Total peaks overlapping transcripts/genes: ' + str(over_over+over_only+under_over+under_only))

    return ('./' + peakname + '/Expressed_' + peakname, './' + peakname + '/NotExpressed_' + peakname)


def TableCreator(peakexon,namepeak='null',lap='null',peaksingenes = 'null', warnings = False):
    exonchromo,exonbegin, exonend, trans_name,exondirection,exontype,peakname,peakchromo,peakbegin, peakend,begingen, endgen=[[],[],[],[],[],[],[],[],[],[],[],[]]
    score, upstream_warning_peak, upstream_warning_tran, upstream_warning_chrom, upstream_warning_peak_start, upstream_warning_peak_end, upstream_warning_FPKM, upstream_warning_direction = [],[],[],[],[],[],[],[]
    downstream_warning_peak, downstream_warning_tran, downstream_warning_chrom, downstream_warning_peak_start, downstream_warning_peak_end, downstream_warning_FPKM, downstream_warning_direction =  [],[],[],[],[],[],[]
    data=open(peakexon)
    a=1
    for lines in data.readlines():
        if a>1:
            exonchromo.append(lines.split('\t')[0])
            exonbegin.append(int(lines.split('\t')[9]))
            exonend.append(int(lines.split('\t')[10]))
            trans_name.append(lines.split('\t')[3])
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
    TTSwarning, TSSwarning = [],[]
    trans_name,exontype = np.array(trans_name),np.array(exontype)
    peakname, exonchromo = [np.array(peakname), np.array(exonchromo)]
    resscore=[]
    uniquepeak=np.unique(peakname)
    for peak in tqdm(range(len(uniquepeak)), desc='TableCreator'):
        indices = np.where(peakname==uniquepeak[peak])
        reschromo.append(exonchromo[indices][0])
        resnamegen.append(trans_name[indices][0])
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
        for indice in range(len(trans_name[indices])):
            if ((peakbegin[indices][indice] <= begingen[indices][indice] and exondirection[indices][indice] == '+') or (
                    peakend[indices][indice] >= endgen[indices][indice] and exondirection[indices][indice] == '-'))and countpeakTSS[matrixlong] != 1:
                countpeakTSS[matrixlong] = 1
            if ((peakbegin[indices][indice] <= begingen[indices][indice] and exondirection[indices][indice] == '-') or (
                    peakend[indices][indice] >= endgen[indices][indice] and exondirection[indices][indice] == '+')) and countpeakTTS[matrixlong] != 1:
                countpeakTTS[matrixlong] = 1

            if (exonbegin[indices][indice] <= peakbegin[indices][indice] <= exonend[indices][indice] or exonbegin[indices][indice] <= peakend[indices][indice] <= exonend[indices][indice] or (
                    exonbegin[indices][indice] >= peakbegin[indices][indice] and peakend[indices][indice] >= exonend[indices][indice])):
                if exontype[indices][indice]=='5UTR' and countpeak5UTR[matrixlong] != 1:
                    countpeak5UTR[matrixlong] = 1
                elif exontype[indices][indice]=='3UTR' and countpeak3UTR[matrixlong] != 1:
                    countpeak3UTR[matrixlong] = 1
                elif countpeakexons[matrixlong] != 1:
                    countpeakexons[matrixlong] = 1
            try:
                if ((exonend[indices][indice] <= peakbegin[indices][indice] <= exonbegin[indices][indice+1] or exonend[indices][indice] <= peakend[indices][indice] <= exonbegin[indices][indice+1] or (exonend[indices][indice] >= peakbegin[indices][indice] and peakend[indices][indice] >= exonbegin[indices][indice]))) and countpeakintrons[matrixlong] != 1:
                    countpeakintrons[matrixlong] = 1
            except:
                pass
            if countpeakintrons[matrixlong] == 0 and countpeak5UTR[matrixlong] == 0 and countpeak3UTR[matrixlong] == 0 and countpeakexons[matrixlong] == 0 and countpeakTSS[matrixlong] == 0 and countpeakTTS[matrixlong] == 0:
                countpeakintrons[matrixlong] = 1
        if countpeakexons[matrixlong] == 0 and countpeakintrons[matrixlong] == 0 and countpeak3UTR[matrixlong] == 0 and countpeak5UTR[matrixlong] == 0:
            if countpeakTSS[matrixlong] == 1:
                upstream_warning_peak.append(uniquepeak[peak])
                upstream_warning_tran.append(trans_name[indices][0])
                upstream_warning_chrom.append(exonchromo[indices][0])
                upstream_warning_direction.append(exondirection[indices][0])
                upstream_warning_FPKM.append(score[indices][0])
                upstream_warning_peak_end.append(peakbegin[indices][0])
                upstream_warning_peak_start.append(peakend[indices][0])
                TSSwarning.append(True)
                TTSwarning.append(False)

            elif countpeakTTS[matrixlong] == 1:
                downstream_warning_peak.append(uniquepeak[peak])
                downstream_warning_tran.append(trans_name[indices][0])
                downstream_warning_chrom.append(exonchromo[indices][0])
                downstream_warning_direction.append(exondirection[indices][0])
                downstream_warning_FPKM.append(score[indices][0])
                downstream_warning_peak_end.append(peakbegin[indices][0])
                downstream_warning_peak_start.append(peakend[indices][0])
                TTSwarning.append(True)
                TSSwarning.append(False)
            else:
                TTSwarning.append(False)
                TSSwarning.append(False)
        else:
            TTSwarning.append(False)
            TSSwarning.append(False)

    ##### Look peaks with antisense transcripts #####

    peaks_dir, peaks_inv_dir, peaks_inv_tran, peaks_tran = {},{},{},{}
    for lines in open(peaksingenes):
        fields = lines.split()
        if fields[5] in peaks_tran.keys():
            if fields[4] not in peaks_dir[fields[5]]:
                peaks_tran.update({fields[5]: peaks_tran[fields[5]]+[fields[3]]})
                peaks_dir.update({fields[5]: peaks_dir[fields[5]]+[fields[4]]})
                peaks_inv_tran[fields[5]] = peaks_tran[fields[5]]
                peaks_inv_dir[fields[5]] = peaks_dir[fields[5]]
            else:
                peaks_tran.update({fields[5]: peaks_tran[fields[5]]+[fields[3]]})
                peaks_dir.update({fields[5]: peaks_dir[fields[5]]+[fields[4]]})
        else:
            peaks_dir.update({fields[5]: [fields[4]]})
            peaks_tran.update({fields[5]: [fields[3]]})
    antisense_warning = []
    for peak in respeakname:
        if peak in peaks_inv_tran.keys():
            antisense_warning.append('True')
        else:
            antisense_warning.append('False')

    ###### -~- ######

    if not os.path.exists(namepeak):
        os.makedirs(namepeak)

    ##### Create the files with the Upstream and Downstream warnings ######
    if warnings == True:
        results=open('./'+namepeak+'/'+namepeak+'_Upstream_Warnings','w')
        results.write('Chromo'+'\t'+'Peak_Start'+'\t'+'Peak_End'+'\t'+'Peak_Name'+'\t'+'ID_transc'+'\t'+'FPKM'+'\t'+'Direction'+'\n')
        for x in range(len(upstream_warning_peak)):
            results.write(upstream_warning_chrom[x] + '\t' + str(upstream_warning_peak_start[x]) + '\t' + str(upstream_warning_peak_end[x]) +'\t' + str(upstream_warning_peak[x]) +'\t' + upstream_warning_tran[x] +'\t' + str(upstream_warning_FPKM[x]) +'\t' + upstream_warning_direction[x] + "\n")
        results.close()

        results=open('./'+namepeak+'/'+namepeak+'_Downstream_Warnings','w')
        results.write('Chromo'+'\t'+'Peak_Start'+'\t'+'Peak_End'+'\t'+'Peak_Name'+'\t'+'ID_transc'+'\t'+'FPKM'+'\t'+'Direction'+'\n')
        for x in range(len(downstream_warning_peak)):
            results.write(downstream_warning_chrom[x] + '\t' + str(downstream_warning_peak_start[x]) + '\t' + str(downstream_warning_peak_end[x]) +'\t' + downstream_warning_peak[x] +'\t' + downstream_warning_tran[x] +'\t' + str(downstream_warning_FPKM[x]) +'\t' + downstream_warning_direction[x] + "\n")
        results.close()

        results=open('./'+namepeak+'/'+namepeak+'_Antisense_Warnings','w')
        results.write('Peak'+'\t'+'transcripts' + '\t' + 'Strand' + '\n')
        for x in peaks_inv_tran.keys():
            results.write(x + '\t' + str(peaks_tran[x])[1:-1][0:] + "\t" + str(peaks_dir[x])[1:-1][0:] + "\n")
        results.close()

    ##### -~- #####

    results=open('./'+namepeak+'/'+namepeak+"_"+lap+'_Annotation.table','w')
    results.write('Chromo'+'\t'+'Peak_Start'+'\t'+'Peak_End'+'\t'+'Peak_Name'+'\t'+'ID_transc'+'\t'+'FPKM'+'\t'+'Direction'+'\t'+'TSS'+'\t'+'5UTR'+'\t'+'Exons'+'\t'+'Introns'+'\t'+'3UTR'+'\t'+'TTS'+'\t'+'Warning_Upstream'+'\t'+'Warning_downstream'+'\t'+'Warning_Antisense'+'\t'+'\n')
    for res in range(len(resnamegen)):
        results.write(reschromo[res]+'\t'+str(respeakstart[res])+'\t'+str(respeakends[res])+'\t'+respeakname[res]+'\t'+resnamegen[res]+'\t'+str(resscore[res])+'\t'+resdirection[res]+'\t'+str(countpeakTSS[res])+'\t'+str(countpeak5UTR[res])+'\t'+str(countpeakexons[res])+'\t'+str(countpeakintrons[res])+'\t'+str(countpeak3UTR[res])+'\t'+str(countpeakTTS[res])+'\t'+str(TSSwarning[res])+'\t'+str(TTSwarning[res])+'\t'+str(antisense_warning[res])+'\t'+'\n')
    results.close()

    return('./'+namepeak+'/'+namepeak+"_"+lap+'_Annotation.table')

# This program take the list of exons that are in the genes and put them in a new list where the each row is a exon with the name of the file
# Need a bed file and for the 2 inputs, been one the list of the genes with the exons size and position and the other the clean list of gain or lost genes.

def FeatureAssign(cleanpeaklist,UTR5='5UTR.bed',UTR3='3UTR.bed',CDE='CDE.bed',peakname='null',lap='null'):
    print('FeatureAssign')
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

def Upsetplotting(name_file,name_output,folder):
    UTR5,Exon,Intron,UTR3,Upstream,Downstream = [],[],[],[],[],[]
    data=open(name_file)
    next(data)
    for lines in data.readlines():
        Upstream.append(lines.split('\t')[7])
        UTR5.append(lines.split('\t')[8])
        Exon.append(lines.split('\t')[9])
        Intron.append(lines.split('\t')[10])
        UTR3.append(lines.split('\t')[11])
        Downstream.append(lines.split('\t')[12])

    Upstream = pd.Series([True if x=="1" else False for x in Upstream])
    UTR5 = pd.Series([True if x=="1" else False for x in UTR5])
    Exon = pd.Series([True if x=="1" else False for x in Exon])
    Intron = pd.Series([True if x=="1" else False for x in Intron])
    UTR3 = pd.Series([True if x=="1" else False for x in UTR3])
    Downstream = pd.Series([True if x=="1" else False for x in Downstream])

    concat = pd.concat([Upstream,UTR5,Exon,Intron,UTR3,Downstream],axis=1,keys=["Upstream","UTR5","Exon","Intron","UTR3","Downstream"])
    result = concat.groupby(["Upstream","UTR5","Exon","Intron","UTR3","Downstream"]).size()
    result = result.nlargest(12)


    plot(result, sort_by = "cardinality")
    pyplot.suptitle("Intersection size")
    pyplot.savefig(folder+"/"+name_output+"_Upsetplot.pdf")

def Dropa_pie_Intergenic(intergenic_file,peaks_in_genes,name_output,folder):
    intergenic_peaks = []
    intra_peaks = []
    data = open(intergenic_file)
    for lines in data.readlines():
        intergenic_peaks.append(lines.split('\t')[3])
    all_array=[len(intergenic_peaks)]
    data = open(peaks_in_genes)
    for lines in data.readlines():
        intra_peaks.append(lines.split('\t')[5])
    all_array.append(len(set(intra_peaks)))
    labels = "Intergenic", "Intragenic"
    fig1, ax1 = pyplot.subplots()
    fig1.suptitle('Proportion of Intergenic and Intragenic peaks')
    ax1.axis('equal')
    wedges, texts, autotexts = ax1.pie(all_array,autopct='%1.1f%%')
    for w in wedges:
        w.set_linewidth(1)
        w.set_edgecolor('black')
    ax1.legend(wedges, [labels[0]+"="+str(all_array[0]),labels[1]+"="+ str(all_array[1])],
              loc="upper right",bbox_to_anchor=(1.1, 1))
    pyplot.savefig(folder+"/"+name_output+"_Intergenic.pdf")

def Dropa_pie_AllwithIntergenic(intergenic_file,All_Anotation_table, name_output, folder, bool_for_intergenic):
    intergenic_peaks = []
    TSS_peaks = []
    UTR5_peaks = []
    Exon_peaks = []
    Intron_peaks = []
    UTR3_peaks = []
    TTS_peaks = []

    data = open(All_Anotation_table)
    for lines in data.readlines():
        TSS_peaks.append(lines.split('\t')[7])
        UTR5_peaks.append(lines.split('\t')[8])
        Exon_peaks.append(lines.split('\t')[9])
        Intron_peaks.append(lines.split('\t')[10])
        UTR3_peaks.append(lines.split('\t')[11])
        TTS_peaks.append(lines.split('\t')[12])

    TSS_ocurrence = TSS_peaks.count("1")
    UTR5_ocurrence = UTR5_peaks.count("1")
    Exon_ocurrence = Exon_peaks.count("1")
    Intron_ocurrence = Intron_peaks.count("1")
    UTR3_ocurrence = UTR3_peaks.count("1")
    TTS_ocurrence = TTS_peaks.count("1")

    if bool_for_intergenic:
        data = open(intergenic_file)
        for lines in data.readlines():
            intergenic_peaks.append(lines.split('\t')[3])
        all_array=len(intergenic_peaks)
        all_array = [all_array, TSS_ocurrence, UTR5_ocurrence, Exon_ocurrence, Intron_ocurrence, UTR3_ocurrence, TTS_ocurrence]
        labels = "Intergenic", "Upstream", "5`UTR", "Exon", "Intron","3`UTR", "Downstream"
        fig1, ax1 = pyplot.subplots()
        fig1.suptitle('Proportion of peaks over gene features')
        ax1.axis('equal')
        wedges, texts, autotexts = ax1.pie(all_array,autopct='%1.1f%%')
        for w in wedges:
            w.set_linewidth(1)
            w.set_edgecolor('black')
        ax1.legend(wedges, [labels[0]+"="+str(all_array[0]),labels[1]+"="+ str(all_array[1]),labels[2]+"="+ str(all_array[2]),labels[3]+"="+ str(all_array[3]),
                            labels[4]+"="+ str(all_array[4]),labels[5]+"="+ str(all_array[5]) , labels[6]+"="+ str(all_array[6])] , loc="upper right",bbox_to_anchor=(1.1, 1.11))
        pyplot.savefig(folder+"/"+name_output+"_Intergenic.pdf")

    ############################################################ PIE WITHOUT THE INTERGENIC ######################################################################################

    all_array2 = [TSS_ocurrence, UTR5_ocurrence, Exon_ocurrence, Intron_ocurrence, UTR3_ocurrence, TTS_ocurrence]
    labels = "Upstream", "5`UTR", "Exon", "Intron","3`UTR", "Downstream"
    fig2, ax2 = pyplot.subplots()
    fig2.suptitle('Proportion of peaks over gene features')
    ax2.axis('equal')
    wedges, texts, autotexts = ax2.pie(all_array2,autopct='%1.1f%%')
    for w in wedges:
        w.set_linewidth(1)
        w.set_edgecolor('black')
    ax2.legend(wedges, [labels[0]+"="+str(all_array2[0]),labels[1]+"="+ str(all_array2[1]),labels[2]+"="+ str(all_array2[2]),labels[3]+"="+ str(all_array2[3]),
                        labels[4]+"="+ str(all_array2[4]),labels[5]+"="+ str(all_array2[5])] , loc="upper right",bbox_to_anchor=(1.1, 1.11))
    pyplot.savefig(folder+"/"+name_output+".pdf")




def Dropa_histogram(intergenic_file,All_Anotation_table, Peak_file, name_output, folder):
    fig1, ax1 = pyplot.subplots()
    intergenic_peaks = []
    TSS_peaks = []
    UTR5_peaks = []
    Exon_peaks = []
    Intron_peaks = []
    UTR3_peaks = []
    TTS_peaks = []
    peak = []
    data = open(intergenic_file)
    for lines in data.readlines():
        intergenic_peaks.append(lines.split('\t')[3])
    Quantity_intergenic=len(intergenic_peaks)
    data = open(All_Anotation_table)
    for lines in data.readlines():
        TSS_peaks.append(lines.split('\t')[7])
        UTR5_peaks.append(lines.split('\t')[8])
        Exon_peaks.append(lines.split('\t')[9])
        Intron_peaks.append(lines.split('\t')[10])
        UTR3_peaks.append(lines.split('\t')[11])
        TTS_peaks.append(lines.split('\t')[12])
    data = open(Peak_file)
    for lines in data.readlines():
        peak.append(lines.split('\t')[3])

    quantity_peaks = len(peak)
    TSS_ocurrence = TSS_peaks.count("1")
    UTR5_ocurrence = UTR5_peaks.count("1")
    Exon_ocurrence = Exon_peaks.count("1")
    Intron_ocurrence = Intron_peaks.count("1")
    UTR3_ocurrence = UTR3_peaks.count("1")
    TTS_ocurrence = TTS_peaks.count("1")

    Intergenic_percentage = round((Quantity_intergenic/quantity_peaks)*100,2)
    TSS_percentage = round((TSS_peaks.count("1")/quantity_peaks)*100,2)
    UTR5_percentage = round((UTR5_peaks.count("1")/quantity_peaks)*100,2)
    Exon_percentage = round((Exon_peaks.count("1")/quantity_peaks)*100,2)
    Intron_percentage = round((Intron_peaks.count("1")/quantity_peaks)*100,2)
    UTR3_percentage = round((UTR3_peaks.count("1")/quantity_peaks)*100,2)
    TTS_percentage = round((TTS_peaks.count("1")/quantity_peaks)*100,2)
    x = np.arange(7)


    all_array = Intergenic_percentage, TSS_percentage, UTR5_percentage, Exon_percentage, Intron_percentage, UTR3_percentage, TTS_percentage
    labels = "Intergenic", "Upstream", "5`UTR", "Exon", "Intron","3`UTR", "Downstream"
    bar_list= pyplot.bar(x,all_array)
    bar_list[0].set_facecolor('b')
    bar_list[1].set_facecolor('g')
    bar_list[2].set_facecolor('r')
    bar_list[3].set_facecolor('c')
    bar_list[4].set_facecolor('m')
    bar_list[5].set_facecolor('y')
    bar_list[6].set_facecolor('b')
    fig1.suptitle('Frequency of peaks over gene features')
    pyplot.ylabel("Fraction (%)")
    pyplot.legend(bar_list, [labels[0]+"="+str(Quantity_intergenic),labels[1]+"="+ str(TSS_ocurrence),labels[2]+"="+ str(UTR5_ocurrence),labels[3]+"="+ str(Exon_ocurrence),labels[4]+"="+ str(Intron_ocurrence),labels[5]+"="+ str(UTR3_ocurrence) , labels[6]+"="+ str(TTS_ocurrence)] , loc="upper right",bbox_to_anchor=(1.1, 1.11))
    pyplot.xticks(x, labels)
    pyplot.savefig(folder+"/"+name_output+"Histogram.pdf")





def Dropa_Enrichment(intergenic_file,All_Anotation_table, number_shuffles, name_output, folder):
    intergenic_peaks = []
    TSS_peaks = []
    UTR5_peaks = []
    Exon_peaks = []
    Intron_peaks = []
    UTR3_peaks = []
    TTS_peaks = []
    data = open(intergenic_file)
    for lines in data.readlines():
        intergenic_peaks.append(lines.split('\t')[3])
    Quantity_intergenic=len(intergenic_peaks)
    data = open(All_Anotation_table)
    for lines in data.readlines():
        TSS_peaks.append(lines.split('\t')[7])
        UTR5_peaks.append(lines.split('\t')[8])
        Exon_peaks.append(lines.split('\t')[9])
        Intron_peaks.append(lines.split('\t')[10])
        UTR3_peaks.append(lines.split('\t')[11])
        TTS_peaks.append(lines.split('\t')[12])

    Original_TSS_ocurrence = TSS_peaks.count("1")
    Original_UTR5_ocurrence = UTR5_peaks.count("1")
    Original_Exon_ocurrence = Exon_peaks.count("1")
    Original_Intron_ocurrence = Intron_peaks.count("1")
    Original_UTR3_ocurrence = UTR3_peaks.count("1")
    Original_TTS_ocurrence = TTS_peaks.count("1")

    enrichment_intergenic_ocurrence, enrichment_TSS_ocurrence, enrichment_UTR5_ocurrence, enrichment_Exon_ocurrence, enrichment_Intron_ocurrence, enrichment_UTR3_ocurrence, enrichment_TTS_ocurrence = [],[],[],[],[],[],[]
    shuffle_TSS,shuffle_UTR5,shuffle_Exon,shuffle_Intron,shuffle_UTR3,shuffle_TTS,shuffle_intergenic = [], [], [], [], [], [], []
    for x in range(1,number_shuffles+1):
        TSS_peaks, UTR5_peaks, Exon_peaks, Intron_peaks, TTS_peaks, UTR3_peaks, intergenic_peaks = [], [], [], [], [], [], []
        data = open(name_output+"shuffle"+str(x)+"/"+name_output+"shuffle"+str(x)+"_All_Annotation.table")
        for lines in data.readlines():
            TSS_peaks.append(lines.split('\t')[7])
            UTR5_peaks.append(lines.split('\t')[8])
            Exon_peaks.append(lines.split('\t')[9])
            Intron_peaks.append(lines.split('\t')[10])
            UTR3_peaks.append(lines.split('\t')[11])
            TTS_peaks.append(lines.split('\t')[12])

        shuffle_TSS_ocurrence = TSS_peaks.count("1")
        shuffle_UTR5_ocurrence = UTR5_peaks.count("1")
        shuffle_Exon_ocurrence = Exon_peaks.count("1")
        shuffle_Intron_ocurrence = Intron_peaks.count("1")
        shuffle_UTR3_ocurrence = UTR3_peaks.count("1")
        shuffle_TTS_ocurrence = TTS_peaks.count("1")
        data = open(name_output+"shuffle"+str(x)+"/"+name_output+"shuffle"+str(x)+"_intergenic.bed")
        for lines in data.readlines():
            intergenic_peaks.append(lines.split('\t')[3])
        shuffle_quantity_intergenic = len(intergenic_peaks)
        Total_peaks = shuffle_quantity_intergenic + len(TTS_peaks)

        enrichment_TSS_ocurrence.append(Original_TSS_ocurrence / shuffle_TSS_ocurrence)
        enrichment_UTR5_ocurrence.append(Original_UTR5_ocurrence / shuffle_UTR5_ocurrence)
        enrichment_Exon_ocurrence.append(Original_Exon_ocurrence / shuffle_Exon_ocurrence)
        enrichment_Intron_ocurrence.append(Original_Intron_ocurrence / shuffle_Intron_ocurrence)
        enrichment_UTR3_ocurrence.append(Original_UTR3_ocurrence / shuffle_UTR3_ocurrence)
        enrichment_TTS_ocurrence.append(Original_TTS_ocurrence / shuffle_TTS_ocurrence)
        enrichment_intergenic_ocurrence.append(Quantity_intergenic / shuffle_quantity_intergenic)

        shuffle_TSS.append(shuffle_TSS_ocurrence)
        shuffle_UTR5.append(shuffle_TSS_ocurrence)
        shuffle_Exon.append(shuffle_TSS_ocurrence)
        shuffle_Intron.append(shuffle_TSS_ocurrence)
        shuffle_UTR3.append(shuffle_TSS_ocurrence)
        shuffle_TTS.append(shuffle_TSS_ocurrence)
        shuffle_intergenic.append(shuffle_TSS_ocurrence)

    average_shuffle_TSS = np.mean(shuffle_TSS)
    average_shuffle_UTR5 = np.mean(shuffle_UTR5)
    average_shuffle_Exon = np.mean(shuffle_Exon)
    average_shuffle_Intron = np.mean(shuffle_Intron)
    average_shuffle_UTR3 = np.mean(shuffle_UTR3)
    average_shuffle_TTS = np.mean(shuffle_TTS)
    average_shuffle_intergenic = np.mean(shuffle_intergenic)

    average_enrichment_TSS = round(sum(enrichment_TSS_ocurrence)/len(enrichment_TSS_ocurrence), 2)
    average_enrichment_UTR5 = round(sum(enrichment_UTR5_ocurrence)/len(enrichment_UTR5_ocurrence), 2)
    average_enrichment_Exon = round(sum(enrichment_Exon_ocurrence)/len(enrichment_Exon_ocurrence), 2)
    average_enrichment_Intron = round(sum(enrichment_Intron_ocurrence)/len(enrichment_Intron_ocurrence), 2)
    average_enrichment_UTR3 = round(sum(enrichment_UTR3_ocurrence)/len(enrichment_UTR3_ocurrence), 2)
    average_enrichment_TTS = round(sum(enrichment_TTS_ocurrence)/len(enrichment_TTS_ocurrence), 2)
    average_enrichment_intergenic = round(sum(enrichment_intergenic_ocurrence)/len(enrichment_intergenic_ocurrence), 2)

    no_TSS_shuffle = int(Total_peaks - average_shuffle_TSS)
    no_UTR5_shuffle = int(Total_peaks - average_shuffle_UTR5)
    no_Exon_shuffle = int(Total_peaks - average_shuffle_Exon)
    no_Intron_shuffle = int(Total_peaks - average_shuffle_Intron)
    no_UTR3_shuffle = int(Total_peaks - average_shuffle_UTR3)
    no_TTS_shuffle = int(Total_peaks - average_shuffle_TTS)
    no_intergenic_shuffle = int(Total_peaks - average_shuffle_intergenic)

    no_TSS_Original= Total_peaks - Original_TSS_ocurrence
    no_UTR5_Original = Total_peaks - Original_UTR5_ocurrence
    no_Exon_Original = Total_peaks - Original_Exon_ocurrence
    no_Intron_Original = Total_peaks - Original_Intron_ocurrence
    no_UTR3_Original = Total_peaks - Original_UTR3_ocurrence
    no_TTS_Original = Total_peaks - Original_TTS_ocurrence
    no_intergenic_Original = Total_peaks - Quantity_intergenic

    TSS_value = np.array([[average_shuffle_TSS,no_TSS_shuffle], [Original_TSS_ocurrence, no_TSS_Original]])
    UTR5_value = np.array([[average_shuffle_UTR5, no_UTR5_shuffle], [Original_UTR5_ocurrence, no_UTR5_Original]])
    Exon_value = np.array([[average_shuffle_Exon, no_Exon_shuffle], [Original_Exon_ocurrence, no_Exon_Original]])
    Intron_value = np.array([[average_shuffle_Intron, no_Intron_shuffle], [Original_Intron_ocurrence, no_Intron_Original]])
    UTR3_value = np.array([[average_shuffle_UTR3, no_UTR3_shuffle], [Original_UTR3_ocurrence, no_UTR3_Original]])
    TTS_value = np.array([[average_shuffle_TTS, no_TTS_shuffle], [Original_TTS_ocurrence, no_TTS_Original]])
    intergenic_value = np.array([[average_shuffle_intergenic, no_intergenic_shuffle], [Quantity_intergenic, no_intergenic_Original]])

    TSS_pvalue = round(stats.chi2_contingency(TSS_value)[1],4)
    UTR5_pvalue = round(stats.chi2_contingency(UTR5_value)[1],4)
    Exon_pvalue = round(stats.chi2_contingency(Exon_value)[1],4)
    Intron_pvalue = round(stats.chi2_contingency(Intron_value)[1],4)
    UTR3_pvalue = round(stats.chi2_contingency(UTR3_value)[1],4)
    TTS_pvalue = round(stats.chi2_contingency(TTS_value)[1],4)
    intergenic_pvalue = round(stats.chi2_contingency(intergenic_value)[1],4)

    std_enrichment_TSS = np.std(enrichment_TSS_ocurrence)
    std_enrichment_UTR5 = np.std(enrichment_UTR5_ocurrence)
    std_enrichment_Exon = np.std(enrichment_Exon_ocurrence)
    std_enrichment_Intron = np.std(enrichment_Intron_ocurrence)
    std_enrichment_UTR3 = np.std(enrichment_UTR3_ocurrence)
    std_enrichment_TTS = np.std(enrichment_TTS_ocurrence)
    std_enrichment_intergenic = np.std(enrichment_intergenic_ocurrence)

    std_dev = [std_enrichment_intergenic,std_enrichment_TSS,std_enrichment_UTR5,std_enrichment_Exon,std_enrichment_Intron,std_enrichment_UTR3,std_enrichment_TTS]

    x = np.arange(7)
    all_array = average_enrichment_intergenic, average_enrichment_TSS, average_enrichment_UTR5, average_enrichment_Exon, average_enrichment_Intron, average_enrichment_UTR3, average_enrichment_TTS
    labels = "Intergenic", "Upstream", "5`UTR", "Exon", "Intron","3`UTR", "Downstream"
    fig1, ax1 = pyplot.subplots()
    fig1.set_size_inches(8, 7)
    fig1.suptitle('Enrichment over expected of peaks over gene features')
    bar_list= pyplot.bar(x,all_array, yerr=std_dev, align='center', alpha = 0.5, capsize=10)
    bar_list[0].set_facecolor('b')
    bar_list[1].set_facecolor('g')
    bar_list[2].set_facecolor('r')
    bar_list[3].set_facecolor('c')
    bar_list[4].set_facecolor('m')
    bar_list[5].set_facecolor('y')
    bar_list[6].set_facecolor('b')
    pyplot.legend(bar_list, [labels[0]+"="+str(average_enrichment_intergenic),labels[1]+"="+ str(average_enrichment_TSS),labels[2]+"="+ str(average_enrichment_UTR5),labels[3]+"="+ str(average_enrichment_Exon),labels[4]+"="+ str(average_enrichment_Intron),labels[5]+"="+ str(average_enrichment_UTR3) , labels[6]+"="+ str(average_enrichment_TTS)] , loc="upper right",bbox_to_anchor=(1.1, 1.11))
    pyplot.xticks(x, labels)
    pyplot.text(x=0-0.4, y=average_enrichment_intergenic + std_enrichment_intergenic + 0.1, s='p_val=' + str(intergenic_pvalue), size=8)
    pyplot.text(x=0.6, y=average_enrichment_TSS +std_enrichment_TSS+ 0.1, s='p_val=' + str(TSS_pvalue), size=8)
    pyplot.text(x=1.6, y=average_enrichment_UTR5 +std_enrichment_UTR5+ 0.1, s='p_val=' + str(UTR5_pvalue), size=8)
    pyplot.text(x=2.6, y=average_enrichment_Exon +std_enrichment_Exon+ 0.1, s='p_val=' + str(Exon_pvalue), size=8)
    pyplot.text(x=3.6, y=average_enrichment_Intron +std_enrichment_Intron+ 0.1, s='p_val=' + str(Intron_pvalue), size=8)
    pyplot.text(x=4.6, y=average_enrichment_UTR3 +std_enrichment_UTR3+ 0.1, s='p_val=' + str(UTR3_pvalue), size=8)
    pyplot.text(x=5.6, y=average_enrichment_TTS +std_enrichment_TTS+ 0.1, s='p_val=' + str(TTS_pvalue), size=8)
    pyplot.ylabel("Enrichment over expected (Fold Change)")
    pyplot.savefig(folder+"/"+name_output+"_Enrichment.pdf")

############################################################################################################################
