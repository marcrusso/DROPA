# DROPA v1.0.0

# Installation
1) Python:

DROPA core is written in Python 3.5 You can install Python directly from the Console with the command:

	sudo apt-get install python3
	sudo apt-get install python-dev

Install Pip to get the packages that DROPA needs for running.

	sudo apt-get install python3-pip

Here is the list of required packages:

  - numpy(v. 1.16.1)
  - pandas(v. 0.24.1)
  - tqdm(v. 4.31.1)
  - intervaltree(v. 3.0.2)
  - matplotlib(v. 3.0.3)
  - upsetplot(v. 0.2.1)
  - argparse(v. 1.4.0)

To install them:

	pip3 install numpy pandas tqdm intervaltree matplotlib argparse
	pip3 install upsetplot

2) Bedtools

Bedtools is required specifically for the creation of random data. To install bedtools:

	sudo apt-get install bedtools

Bedtools shuffle requires a genome size file structured as follow:

	<ChromosomeName><TAB><ChromosomeSize>

For hg19, mm9 and mm10 this file is provided. Check the Wiki for other genome size files creation.

# Usage
DROPA was tested on a machine running Ubuntu OS (vers. 16.04 LTS and 18.04 LTS)

As input, DROPA requires three files, which are:

●	A file containing query peak locations in Browser Extensible Data (BED) format (BED6 and BED12 are supported); BED file must be sorted by chromosome.

●	A reference set containing information about genes features (5’UTR, 3’UTR, Coding Exons) in BED format and a gene reference in BED12 format; Reference sets for hg19 (Ensembl, UCSCgenes and RefSeq) and RefSeq for mm9 and mm10 are provided. Check the Wiki for custom reference sets creation.

●	A gene expression table  containing the name of each gene in the reference set and its expression value (FPKM, TPM, microarray expression value, etc.), structured as follow:

	<GeneName><TAB><ExprValue>

If this file is not provided, then DROPA skips the gene expression evaluation and annotates each query peak to the gene with the largest overlap.


To test the program (from installation folder):

	python3 DROPA_v0.1.py -ref GeneReference/hg19_RefSeq/ -o TEST -ex Test_hg19_RefSeq_Expression -shuffle 2 -gsize GeneReference/hg19.genome Test_hg19_DRIP_peaks.bed 
