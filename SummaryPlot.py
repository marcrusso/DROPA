
from upsetplot import plot
from matplotlib import pyplot
import pandas as pd
import numpy as np



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
    fig1.suptitle('Percent of Intergenic and Intragenic peaks')
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
        fig1.suptitle('Proportion of peaks in parts of the genes')
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
    fig2.suptitle('Proportion of peaks in parts of the genes')
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
    fig1.suptitle('Frequency of peaks in parts of the genes')
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

    shuffle_quantity_intergenic, shuffle_TSS_ocurrence, shuffle_UTR5_ocurrence, shuffle_Exon_ocurrence, shuffle_Intron_ocurrence, shuffle_UTR3_ocurrence, shuffle_TTS_ocurrence = [],[],[],[],[],[],[]

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

        shuffle_TSS_ocurrence.append(TSS_peaks.count("1"))
        shuffle_UTR5_ocurrence.append(UTR5_peaks.count("1"))
        shuffle_Exon_ocurrence.append(Exon_peaks.count("1"))
        shuffle_Intron_ocurrence.append(Intron_peaks.count("1"))
        shuffle_UTR3_ocurrence.append(UTR3_peaks.count("1"))
        shuffle_TTS_ocurrence.append(TTS_peaks.count("1"))
        data = open(name_output+"shuffle"+str(x)+"/"+name_output+"shuffle"+str(x)+"_intergenic.bed")
        for lines in data.readlines():
            intergenic_peaks.append(lines.split('\t')[3])
        shuffle_quantity_intergenic.append(len(intergenic_peaks))

    average_TSS_ocurrence = sum(shuffle_TSS_ocurrence)/len(shuffle_TSS_ocurrence)
    average_UTR5_ocurrence = sum(shuffle_UTR5_ocurrence)/len(shuffle_UTR5_ocurrence)
    average_Exon_ocurrence = sum(shuffle_Exon_ocurrence)/len(shuffle_Exon_ocurrence)
    average_Intron_ocurrence = sum(shuffle_Intron_ocurrence)/len(shuffle_Intron_ocurrence)
    average_UTR3_ocurrence = sum(shuffle_UTR3_ocurrence)/len(shuffle_UTR3_ocurrence)
    average_TTS_ocurrence = sum(shuffle_TTS_ocurrence)/len(shuffle_TTS_ocurrence)
    average_quantity_intergenic = sum(shuffle_quantity_intergenic)/len(shuffle_quantity_intergenic)

    enrichment_TSS_ocurrence = round(Original_TSS_ocurrence /average_TSS_ocurrence, 2)
    enrichment_UTR5_ocurrence = round(Original_UTR5_ocurrence /average_UTR5_ocurrence, 2)
    enrichment_Exon_ocurrence = round(Original_Exon_ocurrence / average_Exon_ocurrence, 2)
    enrichment_Intron_ocurrence = round(Original_Intron_ocurrence / average_Intron_ocurrence, 2)
    enrichment_UTR3_ocurrence = round(Original_UTR3_ocurrence / average_UTR3_ocurrence, 2)
    enrichment_TTS_ocurrence = round(Original_TTS_ocurrence / average_TTS_ocurrence, 2)
    enrichment_intergenic_ocurrence = round(Quantity_intergenic / average_quantity_intergenic, 2)

    x = np.arange(7)
    all_array = enrichment_intergenic_ocurrence, enrichment_TSS_ocurrence, enrichment_UTR5_ocurrence, enrichment_Exon_ocurrence, enrichment_Intron_ocurrence, enrichment_UTR3_ocurrence, enrichment_TTS_ocurrence
    labels = "Intergenic", "Upstream", "5`UTR", "Exon", "Intron","3`UTR", "Downstream"
    fig1, ax1 = pyplot.subplots()
    fig1.set_size_inches(8, 7)
    fig1.suptitle('Enrichment over expected of peaks in gene regions')
    bar_list= pyplot.bar(x,all_array)
    bar_list[0].set_facecolor('b')
    bar_list[1].set_facecolor('g')
    bar_list[2].set_facecolor('r')
    bar_list[3].set_facecolor('c')
    bar_list[4].set_facecolor('m')
    bar_list[5].set_facecolor('y')
    bar_list[6].set_facecolor('b')
    pyplot.legend(bar_list, [labels[0]+"="+str(enrichment_intergenic_ocurrence),labels[1]+"="+ str(enrichment_TSS_ocurrence),labels[2]+"="+ str(enrichment_UTR5_ocurrence),labels[3]+"="+ str(enrichment_Exon_ocurrence),labels[4]+"="+ str(enrichment_Intron_ocurrence),labels[5]+"="+ str(enrichment_UTR3_ocurrence) , labels[6]+"="+ str(enrichment_TTS_ocurrence)] , loc="upper right",bbox_to_anchor=(1.1, 1.11))
    pyplot.xticks(x, labels)
    pyplot.savefig(folder+"/"+name_output+"_Enrichment.pdf")