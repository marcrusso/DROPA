# DROPA


DROPA uses Python3 (vers. Python 3.6.4) and R (vers 3.5.1).
Python modules required:
•	numpy (vers 1.14.0)
•	tqdm (vers 4.23.4)
•	pandas (vers 0.22.0)
•	intervaltree (vers 2.1.0)
R libraries required: 
•	UpSetR (vers. 1.3.2)

DROPA was tested on a machine running Ubuntu OS (vers. 16.04 LTS)

As input, DROPA requires three files, which are:

●	A file containing query peak locations in Browser Extensible Data (BED) format (BED6 and BED12 are supported);

●	A reference set containing information about genes features (5’UTR, 3’UTR, exon, intron, etc.) in BED format and a gene reference in Gene Transfer Format (GTF);

●	A 2-column gene expression table  containing the name of each gene and its expression value (FPKM, TPM, microarray expression value, etc.). If this file is not provided, then DROPA skips the gene expression evaluation and annotates each query peak to the gene with the largest overlap.
