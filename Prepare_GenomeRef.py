
import argparse

# The program takes as an Input the folder where the genome reference is

############ INPUT OF VARIABLES ###################
parser = argparse.ArgumentParser(description='Input program')
parser.add_argument('namefolder', type=str, help='Folder with the genome reference')
folder = parser.parse_args().namefolder


#################################################### 3UTR ############################################################

result_name,line0,line1,line2,line4,line5 = [],[],[],[],[],[]
data = open(folder+"/3UTR.bed", "r")
for line in data.readlines():
    complete_name = (line.split('\t')[3])
    partial_name = complete_name.split('_utr3')[0]
    result_name.append(partial_name)
    line0.append(line.split('\t')[0])
    line1.append(line.split('\t')[1])
    line2.append(line.split('\t')[2])
    line4.append(line.split('\t')[4])
    line5.append(line.split('\t')[5])

data.close()
data = open(folder+"/3UTR.bed", "w")

for x in range(len(line0)):
    data.write(line0[x] + '\t' + line1[x] + '\t' +line2[x] + '\t' + result_name[x] + '\t' + line4[x] + '\t' + line5[x])
data.close()

##################################################### 5UTR ############################################################

result_name,line0,line1,line2,line4, line5 = [],[],[],[],[],[]
data = open(folder+"/5UTR.bed", "r")
for line in data.readlines():
    complete_name = (line.split('\t')[3])
    partial_name = complete_name.split('_utr3')[0]
    result_name.append(partial_name)
    line0.append(line.split('\t')[0])
    line1.append(line.split('\t')[1])
    line2.append(line.split('\t')[2])
    line4.append(line.split('\t')[4])
    line5.append(line.split('\t')[5])

data.close()
data = open(folder+"/5UTR.bed", "w")

for x in range(len(line0)):
    data.write(line0[x] + '\t' + line1[x] + '\t' +line2[x] + '\t' + result_name[x] + '\t' + line4[x] + '\t' + line5[x])
data.close()

##################################################### CDE ############################################################

result_name,line0,line1,line2,line4, line5 = [],[],[],[],[],[]
data = open(folder+"/CDE.bed", "r")
for line in data.readlines():
    complete_name = (line.split('\t')[3])
    partial_name = complete_name.split('_cds')[0]
    result_name.append(partial_name)
    line0.append(line.split('\t')[0])
    line1.append(line.split('\t')[1])
    line2.append(line.split('\t')[2])
    line4.append(line.split('\t')[4])
    line5.append(line.split('\t')[5])

data.close()
data = open(folder+"/CDE.bed", "w")

for x in range(len(line0)):
    data.write(line0[x] + '\t' + line1[x] + '\t' +line2[x] + '\t' + result_name[x] + '\t' + line4[x] + '\t' + line5[x])
data.close()