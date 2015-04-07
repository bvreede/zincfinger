#!/usr/bin/python

'''
This script can be used to take a csv file where clusters are indicated (in column 4 as a number) and
use it to make individual fasta files and visualizations of the genes in this cluster.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 1 April 2015
'''

import csv, sys

maindb = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters"
infile = "nameclusters6.csv" #csv file with [geneID,genename,proteinID,clusternumber]

#fasta file with sequences
fasta = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-allz_seq.fasta"

'''
Read the input from the cluster file (= gene names and their assigned cluster)
and save to local memory.
'''
infileo = csv.reader(open("%s/%s" %(maindb,infile)))

clusternames,clusterinfo,singlecluster = [],[],[]
curcluster = ''
#singlecluster will collect information per cluster
#clusternames will save all cluster names, in order of appearance
#clusterinfo is the total 3d matrix with clusters in D1 (0), genes in D2 (1), and ID,name,protID in D3 (2)
for line in infileo:
	if line[3] in clusternames:
		if curcluster != line[3]:
			sys.exit("Script is not robust to clusters that are not in order. Sort and try again!")
		else:
			singlecluster.append(line[0:3])
	else:
		clusternames.append(line[3])
		if len(singlecluster) >= 1:
			clusterinfo.append(singlecluster)
		singlecluster = []
		singlecluster.append(line[0:3])
		curcluster = line[3]
clusternames.append(line[3])
if len(singlecluster) >= 1:
	clusterinfo.append(singlecluster)

'''
Save all fasta files to memory so sequences can be recalled per cluster.
'''
fastadict = {}
def fastadicter(fastadb):
	'''
	Takes an open fastafile and adds the sequences to a general sequences dictionary, with the
	headers as keys.
	'''
	global fastadict
	sequence = ""
	header = ""
	for line in fastadb:
		if line[0] == ">":
			if header != "":
				fastadict[header] = sequence
			header = line[1:].strip().replace(':','-')
			sequence = ""
		else:
			sequence += line.strip()
		fastadict[header] = sequence

fo = open(fasta)
fastadicter(fo)


'''
loutput:
CHECK per cluster: fasta file with sequences (send to Ariel and have him make sequence trees)
- per cluster: motif scores in separate file (and run the viz script on this)
- Drosophila genes only; gorilla input list
- Drosophila genes only; with gorilla output: make visualization (heatmap)

'''

for i,c in enumerate(clusterinfo):
	clusterid = clusternames[i]
	#generate a fastafile per cluster with only the genes from this cluster
	cfile = "%s/cluster-%s.fa" %(maindb,clusterid)
	out = open(cfile, "w")
	for n in c:
		name = "%s|%s|%s" %(n[0],n[1],n[2])
		out.write(">%s\n" %name)
		out.write(fastadict[name])
		out.write("\n\n")
	out.close()
	


