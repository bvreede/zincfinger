#!/usr/bin/python

'''
Script to run findmotif.py twice: once on each species individually,
and only to identify which motifs pass the determined threshold
and will be used for further analysis; then a second run to actually
make the stats given these motifs.
'''
import os,config
import seaborn as sns
import pylab as pl
import numpy as np

#CUSTOMIZE: selection criteria. A motif will be chosen if:
# it is found at least once in XXX or more individual species;
a=15
#OR
# it is found more than XXX times in a single species
b=5

bargraphfig = "%s/%s/%s_perspecies" %(config.mainfolder,config.imgfolder,config.idr)
heatmapfig = "%s/%s/%s_perspecies" %(config.mainfolder,config.imgfolder,config.idr)


#set figure aesthetics
sns.set_style("white")
sns.set_context("poster")



# run findmotif on all files in the input folder, and concatenate all files to a single input file
confile = "%s/%s/%s-alli_seq.fa" %(config.mainfolder,config.dbfolder,config.idr)
concat = open(confile, "w")
for f in os.listdir("%s/%s" %(config.mainfolder,config.seqpfolder)):
	os.system("findmotif.py %s/%s/%s screen" %(config.mainfolder,config.seqpfolder,f))
	# add file content to concatenated sequence file
	fo  = open("%s/%s/%s" %(config.mainfolder,config.seqpfolder,f))
	for line in fo:
		concat.write(line)
	fo.close()
concat.close()

# concatenate the hmm files
confilehmm = "%s/%s/%s-alli_hmmsearch.txt" %(config.mainfolder,config.hmmfolder,config.idr)
concathmm = open(confilehmm, "w")
for f in os.listdir("%s/%s" %(config.mainfolder,config.hmmfolder)):
	if "alli_hmmsearch.txt" in f:
		continue
	# add file content to concatenated sequence file
	fo  = open("%s/%s/%s" %(config.mainfolder,config.hmmfolder,f))
	for line in fo:
		concathmm.write(line)
	fo.close()
concathmm.close()


#choose the motifs: first read the outputfiles from the first findmotif screen
motifdict,totaldict,nadict = {},{},{}
for f in os.listdir("%s/%s" %(config.mainfolder,config.resfolder)):
	if "_motifcount.txt" in f:
		spp = f.split('_')[0][-4:]
		totaldict[spp],nadict[spp] = {},{}
		fo = open("%s/%s/%s" %(config.mainfolder,config.resfolder,f))
		for line in fo:
			mo = line.split('_na: ')
			if len(mo) < 2: #total hits
				mo = line.split('_total: ')
				totaldict[spp][mo[0]] = mo[1].strip()
				continue
			# only non-ambiguous hits
			nadict[spp][mo[0]] = mo[1].strip()
			if mo[0] in motifdict:
				motifdict[mo[0]] = motifdict[mo[0]] + [int(mo[1].strip())]
			else:
				motifdict[mo[0]] = [int(mo[1].strip())]

# choose motifs from the dictionary
selected_motifs=open(config.semot, "w")
motifs = []
for m in motifdict:
	cs = motifdict[m]
	nnz = len(cs) - cs.count(0) #the number of entries with non-zero hits
	if max(cs) >= b:	# more than 5 hits in a single species
		motifs.append(m)
	elif nnz >= a: # hits in more than 10 species
		motifs.append(m)

# sort the motifs
morder = []
for m in motifs:
	ms = m.split('_')
	# enter first item in morder:
	if len(morder) == 0:
		morder.append(m)
		continue
	# now, add items to morder
	for i,k in enumerate(morder):
		mk = k.split('_')
		if int(ms[0]) > int(mk[0]):
			if i==len(morder)-1:
				morder.append(m)
				break
			continue
		elif int(ms[0]) == int(mk[0]):
			if int(ms[1]) > int(mk[1]):
				if i==len(morder)-1:
					morder.append(m)
					break
				continue
			elif int(ms[1]) == int(mk[1]):
				if int(ms[2]) > int(mk[2]):
					if i==len(morder)-1:
						morder.append(m)
						break
					continue
				else:
					morder[i:i] = [m]
					break
			else:
				morder[i:i] = [m]
				break
		else:
			morder[i:i] = [m]
			break

for m in morder:
	selected_motifs.write("%s\n" %m)
selected_motifs.close()
		

# make heatmap and bargraph
def makestackedbargraph(values1,values2,labels,name):
	'''
	Makes a stacked bar chart with two sets of values on y and labels on x.
	'''
	fig = pl.figure()
	ind = np.arange(len(values1))
	p1 = pl.bar(ind,values1,color="navy")
	p2 = pl.bar(ind,values2,color="darkkhaki",bottom=values1)
	pl.xticks(ind + 0.5, labels, rotation=90)
	#pl.yticks(np.arange(30000,60000,5000)) #only used to restrict plot values
	#pl.ylim((30000,60000)) #only used to restrict plot values
	#pl.yticks(np.arange(0,9000,1000)) #only used to restrict plot values
	#pl.ylim((0,9000)) #only used to restrict plot values
	pl.savefig("%s-%s.svg" %(bargraphfig,name))
	pl.clf()
	pl.close()

def makeheatmap(matrix,labelsx,labelsy,name):
	'''
	Makes a heatmap of a matrix, using Seaborn.
	'''
	sns.heatmap(matrix, xticklabels=labelsx,yticklabels=labelsy,linewidths=.5)#,cmap="YlOrBr")
	pl.savefig("%s-%s.svg" %(heatmapfig,name))
	pl.clf()
	pl.close()

# data for bargraph
namot,ambmot,totalmot = [],[],[]
for s in config.sppall:
	na = 0
	tot = 0
	for m in morder:
		na += int(nadict[s][m])
		tot += int(totaldict[s][m])
	namot.append(na)
	ambmot.append(tot-na)
	totalmot.append(tot)
	totaldict[s]["total"] = tot


# data for heatmap
heatmatrix = []
for s in config.sppall:
	heatrow = []
	for m in morder:
		luton = int(totaldict[s][m])/float(totaldict[s]["total"])
		heatrow.append(luton)
	heatmatrix.append(heatrow)

#make the images
makeheatmap(heatmatrix,morder,config.sppall,"motifs-proportion")
makestackedbargraph(namot,ambmot,config.sppall,"motifs-stacked")


# run findmotif.py on the concatenated file, without the screen option
os.system("findmotif.py %s" %confile)


