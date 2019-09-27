####################################################
"""
BAM to text file : correct GC and generate RD
parameters:  binSize = 1000
input file: bam = "sim1_6_6100_read.sort.bam"
reference file: 



"""
#####################################################
import numpy as np
import pysam
import math
import sys
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import os
import pandas as pd
import itertools
from sklearn import mixture
from scipy import linalg
import matplotlib as mpl
import copy


def get_chrlist(filename):
    samfile = pysam.AlignmentFile(filename, "rb")
    List = samfile.references
    chrList = np.full(len(List), 0)
    for i in range(len(List)):
        chr = str(List[i]).strip('chr')
        if chr.isdigit():
            chrList[i] = int(chr)
    return chrList
    ####chrList: array([21])

def read_ref_file(filename, ref, num):
    # read reference file
    if os.path.exists(filename):
        print("Read reference file: " + str(filename))
        with open(filename, 'r') as f:
            line = f.readline()
            for line in f:
                linestr = line.strip()
                ref[num] += linestr
    else:
        print("Warning: can not open " + str(filename) + '\n')
    return ref
    ### ref is the refernece sequence

def get_RC(filename, chrList, ReadCount):
    samfile = pysam.AlignmentFile(filename, "rb")
    for line in samfile:
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            if chr.isdigit():
                num = np.argwhere(chrList == int(chr))[0][0]
                posList = line.positions
                ReadCount[num][posList] += 1
    return ReadCount


def ReadDepth(ReadCount, binNum, ref):
	RD = np.full(binNum, 0.0)
	GC = np.full(binNum, 0)
	pos = np.arange(1, binNum+1)
	for i in range(binNum):
		RD[i] = np.mean(ReadCount[i*binSize:(i+1)*binSize])
		cur_ref = ref[i*binSize:(i+1)*binSize]
		N_count = cur_ref.count('N') + cur_ref.count('n')
		if N_count == 0:
			gc_count = cur_ref.count('C') + cur_ref.count('c') + cur_ref.count('G') + cur_ref.count('g')
		else:
			RD[i] = -10000
			gc_count = 0
		GC[i] = int(round(gc_count / binSize, 3) * 1000)
	index = RD > 0
	RD = RD[index]
	GC = GC[index]
	pos = pos[index]
	RD = gc_correct(RD, GC)
	return pos, RD, GC/1000


def gc_correct(RD, GC):
    # correcting gc bias
    bincount = np.bincount(GC)
    global_rd_ave = np.mean(RD)
    for i in range(len(RD)):
        if bincount[GC[i]] < 2:
            continue
        mean = np.mean(RD[GC == GC[i]])
        RD[i] = global_rd_ave * RD[i] / mean
    return RD

def plot(pos, data):
    plt.scatter(pos, data, s=3, c="black")
    #plt.scatter(pos1, data1, s=3, c="red")
    plt.xlabel("ds")
    plt.ylabel("rd")
    plt.show()

def plotgc(pos, data):
    plt.scatter(pos, data, s=3, c="blue")
    #plt.scatter(pos1, data1, s=3, c="red")
    plt.xlabel("ds")
    plt.ylabel("gc")
    plt.show()

def plotGCscatter(GC, RD):
    plt.scatter(GC, RD, s=1,linewidths=0.1)
    plt.xlabel("GC_content")
    plt.ylabel("rd")
    plt.show()

def write_RD(filename, chr, pos, RD, GC):
    output = open(filename, "w")
    for i in range(len(pos)):
        output.write('chr'+str(chr) +"\t" + str(pos[i] * binSize + 1) + '\t' + str((pos[i]+1) * binSize) + '\t' + str(RD[i]) +'\t' +str(GC[i]) +'\n')

def segment(pos, tv_rd):
    start = []
    end = []
    seg_rd = []
    i = 0
    j = 1
    while j < len(tv_rd):
        if j == len(tv_rd) - 1:
            start.append(int(pos[i]))
            end.append(int(pos[j]))
            seg_rd.append(float(tv_rd[i]))
            j += 1
        else:
            if tv_rd[i] == tv_rd[j]:
                j += 1
            else:
                start.append(int(pos[i]))
                end.append(int(pos[j]))
                seg_rd.append(float(tv_rd[i]))
                i = j
    return start, end, seg_rd

def rd_matrix(RD):
    # convert the RD data to two dim matrix (for seg_rd)
    RD = RD.astype(np.float)
    pos = np.array(range(1, len(RD)+1))
    nr_min = np.min(RD)
    nr_max = np.max(RD)
    newpos = (pos - min(pos)) / (max(pos) - min(pos)) * (nr_max - nr_min) + nr_min          ###(add pos a charactor)
    newpos = newpos.astype(np.float)
    rd = np.c_[newpos, RD]
    return rd

def plot_dpGMM_result(dpGMM, rd):
	##plot the dpGMM result 
	color_iter = itertools.cycle(['r','g','b','c','m','k','purple','pink','cyan'])
	Y_ = dpGMM.fit_predict(rd)
	for i, (mean, covar, color) in enumerate(zip(dpGMM.means_,dpGMM.covariances_,color_iter)):
		if not np.any(Y_ == i):
			continue
		plt.scatter(rd[Y_==i,0],rd[Y_==i,1],.8,color = color)
		plt.title('Dirichlet process mixture model')
	plt.show()

def compute_CNV(filename, dpGMM, chr, label, start, end, seg_rd):
	RD_mean = np.dot(dpGMM.means_[:,1], dpGMM.weights_)
	print("RD_mean: ", RD_mean)
	## determine the boundary
	lowerRD = 0.75 * RD_mean
	upperRD =  1.25 * RD_mean
	la = np.where(dpGMM.means_[:,1] < lowerRD)[0]
	ga = np.where(dpGMM.means_[:,1] > upperRD)[0]
	if la.size == 0 and ga.size != 0:
		print("There is no loss")
		gainlabel = ga[0]
		gainindexseg = np.where(label == gainlabel)
		cnvindex = gainindexseg[0]
	elif ga.size == 0 and la.size != 0:
		print("There is no gain")
		losslabel = la[0]
		lossindexseg = np.where(label == losslabel)
		cnvindex = lossindexseg[0]
	elif la.size == 0 and ga.size ==0:
		print("Sorry, there is no variation")
		sys.exit(0)
	else:
		losslabel = la[0]
		gainlabel = ga[0]
		lossindexseg = np.where(label == losslabel)
		gainindexseg = np.where(label == gainlabel)
		cnvindex = np.hstack((lossindexseg[0],gainindexseg[0]))       ## the index of cnv
		cnvindex.sort()
	CNV = np.zeros(len(cnvindex))
	for i in range(len(cnvindex)):
		CNV[i] = round(2 * seg_rd[cnvindex][i] / RD_mean)
	type = ""
	output = open(filename, "w")
	output.write('chr'+str(chr) +"\t" + 'start' + '\t' + 'end' + '\t' + 'seg_rd' + '\t' + 'cn' + '\t' + 'type' + '\n')
	CNV_start = []
	CNV_end = []
	cn = []
	for i in range(len(cnvindex)):
		if CNV[i] < 2:
			CNV_start.append(start[cnvindex[i]])
			CNV_end.append(end[cnvindex[i]])
			cn.append(CNV[i])
			type = 'loss'
			output.write('chr'+str(chr) +"\t" + str(start[cnvindex[i]] * binSize + 1) + '\t' + str(end[cnvindex[i]] * binSize) + '\t' + str(seg_rd[cnvindex[i]]) + '\t' + str(CNV[i]) + '\t' + type + '\n')
			
		if CNV[i] > 2:
			CNV_start.append(start[cnvindex[i]])
			CNV_end.append(end[cnvindex[i]])
			cn.append(CNV[i])
			type = 'gain'
			output.write('chr'+str(chr) +"\t" + str(start[cnvindex[i]] * binSize + 1) + '\t' + str(end[cnvindex[i]] * binSize) + '\t' + str(seg_rd[cnvindex[i]]) + '\t' + str(CNV[i]) + '\t' + type + '\n')
			
	return CNV_start, CNV_end, cn


def CNV_writer(filename, chr, CNV_start, CNV_end, cn):
	start11 = []
	end11 = []
	CNV11 = []
	count1 = 0
	for i in range(len(CNV_start)-1):
		if CNV_end[i] == CNV_start[i+1] and cn[i] == cn[i+1]:
			count1 += 1
			CNV_end[i] = CNV_end[i+1]
		else:
			start11.append(CNV_start[(i-count1)])
			end11.append(CNV_end[i])
			CNV11.append(cn[i])
			count1 = 0
	lastIndex = len(CNV_start) - 2
	if CNV_end[lastIndex] == CNV_start[lastIndex + 1] and cn[lastIndex] == cn[lastIndex + 1]:
		start11.append(CNV_start[lastIndex - count1])
		end11.append(CNV_end[lastIndex])
		CNV11.append(cn[lastIndex])
	else:
		start11.append(CNV_start[lastIndex + 1])
		end11.append(CNV_end[lastIndex + 1])
		CNV11.append(cn[lastIndex + 1])
	#return start11, end11, CNV11	
	CNtype = ""
	output = open(filename, "w")
	output.write('chr'+str(chr) +"\t" + 'start' + '\t' + 'end' + '\t'  + 'cn' + '\t' + 'CNtype' + '\n')
	for y in range(len(start11)):
		if CNV11[y] > 2.0:
			CNtype = 'gain'
		if CNV11[y] < 2.0:
			CNtype = 'loss'
		output.write('chr'+str(chr) +"\t" + str(start11[y]* binSize + 1) + '\t' + str(end11[y]* binSize) + '\t'  + str(CNV11[y]) + '\t' + CNtype + '\n')



##get parameters
#bam = "sim1_6_6100_read.sort.bam"
bam = sys.argv[1]
binSize = 1000
#binSize = int(sys.argv[2])   ###1000
CNVfile = bam + '_result.txt'
#print(chrList)
#print(chrNum)
#reference = "chr21.fa"		
reference = sys.argv[3]
if reference.split("chr"):
	refList = [[int(reference.split("chr")[1].split(".")[0])]]
	a = copy.deepcopy(refList)
	chrList = np.array(a[0])
	chrNum = len(chrList)
else:
	chrList = get_chrlist(bam)
	chrList = chrList[0:22]
	chrNum = len(chrList)
	refList =  [[] for i in range(chrNum)]


for i in range(chrNum):
    refList = read_ref_file(reference, refList, i)

chrLen = np.full(chrNum,0)

for i in range(chrNum):
    chrLen[i] = len(refList[i])    ##len(refList[i])=48129895

print("Read bam file:", bam)

ReadCount = np.full((chrNum, np.max(chrLen)), 0)       ####ReadCount.shape=(1, 83257442)(chrNum, np.max(chrLen)),qie all is 0
ReadCount = get_RC(bam, chrList, ReadCount)			###type(ReadCount)=<class 'numpy.ndarray'>

for i in range(chrNum):
    binNum = int(chrLen[i]/binSize)+1			###48130
    pos, RD, GC = ReadDepth(ReadCount[i], binNum, refList[i])           ####refList[i] is the ref21.fa sequence, ReadCount[0] is the read count of bam
    #plot(pos, RD)
    #plotGCscatter(GC, RD)
    ##write_RD(bam+"_RD", chrList[i], pos, RD, GC)
    ### load R package to segment
    v = robjects.FloatVector(RD)
    m = robjects.r['matrix'](v, ncol=1)
    importr('cghFLasso')
    aa1=robjects.r.cghFLasso(m)
    #robjects.r.png("seg.png")
    #robjects.r.plot(aa1,index=1,type="Lines")    ###segment smoothing figure after cghFLasso
    #robjects.r.title("the smoothed by cghFLasso")
    #robjects.r.dev.off()
    segRD = np.array(aa1[0])         ##aa1[0] = aa$Esti.CopyN in R script
    start, end, seg_rd = segment(pos, segRD)
    start = np.array(start)
    end = np.array(end)
    seg_rd = np.array(seg_rd)
    
    print("prepare the readdepth data to use dpcluster")
    rd = rd_matrix(seg_rd) 
    dpGMM = mixture.BayesianGaussianMixture(n_components = 4,covariance_type='full',max_iter = 1000)
    dpGMM.fit_predict(rd)
    #plot_dpGMM_result(dpGMM, rd)
    label = dpGMM.predict(rd)  
    CNV_start, CNV_end, cn = compute_CNV(CNVfile, dpGMM, chrList[i], label, start, end, seg_rd)
    print("write the CNV result to outfile")
    filename1 = bam + "CNVfinal_result"
    CNV_writer(filename1, chrList[i], CNV_start, CNV_end, cn)
    print("Well done!")
