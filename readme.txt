dpGMM, Version 1.0
Installation and Usage

INSTALLTION
NONE, immediately use it by python

PREPARATIONS:
linux system
python 3.7
some python modules: apt-get install numpy, pysam, math, sys, matplotlib, rpy2, sklearn, scipy
R version 3.4.4 and later
cghFLasso (R packages)


USAGE

dpGMM is intended to be esay to use. For typical usage with a human genome, dpGMM requires 3 input parameters:
1) BAM file (sorted)
2) window size
3) Reference file (in fasta format)

Note, the BAM file must be sorted.

Example:
python dpGMM.py bam.file windowsize ref.fa


Example using the test data:

python dpGMM.py sim1_6_6100_read.sort.bam 1000 chr21.fa

The output in terminal:
Read reference file: chr21.fa
Read bam file: sim1_6_6100_read.sort.bam
[1] 1
prepare the readdepth data to use dpcluster
RD_mean:  5.940845849724413
write the CNV result to outfile
Well done!

The output file "sim1_6_6100_read.sort.bamCNVfinal_result" including CNVs is as follows:
chr21	start	end	cn	CNtype
chr21	32001001	32003000	5.0	gain
chr21	32003001	32011000	6.0	gain
chr21	32501001	32511000	6.0	gain
chr21	36001001	36021000	4.0	gain
chr21	36501001	36521000	4.0	gain
chr21	37001001	37051000	3.0	gain
chr21	37501001	37551000	3.0	gain
chr21	41000001	41021000	1.0	loss
chr21	41501001	41521000	1.0	loss
chr21	42001001	42053000	1.0	loss
chr21	42499001	42551000	1.0	loss
chr21	46001001	46021000	1.0	loss
chr21	46501001	46521000	1.0	loss
chr21	47001001	47010000	1.0	loss
chr21	47010001	47046000	0.0	loss
chr21	47046001	47051000	1.0	loss
chr21	47551001	47552000	1.0	loss

In addition, users can also set other parameters in the source code "dpGMM.py" file.
