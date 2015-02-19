#!/usr/bin/python
'''
This script can be used to generate a heatmap to visualize GO terms associated
to a group of genes. It does NOT calculate significance or any other measures;
it needs to work with an input list of genes and GO terms, and find the corresponding
hits in a csv database (downloaded from Ensembl BioMart and prepared with tableprep.py).
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 18 February 2015
'''
import csv, sys
from os import path
import matplotlib.pyplot as plt
import numpy as np

#specify folder for inputfiles, and an errormessage given on every usage abort to the user
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/"
datafolder = "%sresults/singlemotif/" %dbfolder
gosource = "%sdatabases/140720-SM00355-dmel2.csv" %dbfolder
errormess = "USAGE: goheat.py motif_incl/excl (e.g.: goheat.py 2_12_4_excl)\nThe program will automatically find the correct input files in %s" %datafolder


if len(sys.argv) <= 1:
	sys.exit(errormess)

#generate names of input files: the user only indicates the prefix
genefile = datafolder + sys.argv[1] + ".csv"
gofile = datafolder + sys.argv[1] + "_results.csv"

#check if all input files are correct (if they don't exist: abort)
if not path.exists(genefile):
	sys.exit("ABORT: at least one of the input files was not found.\n%s" %errormess)
if not path.exists(gofile):
	sys.exit("ABORT: at least one of the input files was not found.\n%s" %errormess)
if not path.exists(gosource):
	sys.exit("ABORT: could not locate the GO database. Correct the path in the script.")

#open the files
genes = csv.reader(open(genefile))
goterms = csv.reader(open(gofile))
godb = csv.reader(open(gosource))

#collect genes and GO terms
goli, gnli, pidli = [],[],[] #lists of (unique) go IDs, gene names, protein IDs
godict,genedict = {},{} #translation dictionaries to find gene names with protein IDs, and go description with go ID
for n,line in enumerate(genes):
	if n == 0:
		continue
	genedict[line[2]] = line[1]
	gnli.append(line[1])
	pidli.append(line[2])
gnli = list(set(gnli))
pidli = list(set(pidli))

for n,line in enumerate(goterms):
	if n == 0:
		continue
	goli.append(line[0])
	godict[line[0]] = line[1]


#from the GO-file: get data on which genes are associated with which GO terms
#make a table (gene x GO)
gohitdb = []
gohitdb2 = []
for n,line in enumerate(godb):
	#if n == 0:
	#	continue
	#if n == 10:
	#	break
	#print line
	if line[2] in pidli and line[6] in goli:
		gohit = frozenset([line[1],line[6]])
		gohit2 = [line[1],line[6]]
		gohitdb.append(gohit)
		gohitdb2.append(gohit2)
gohitdb = set(gohitdb)


for hit in gohitdb:
	print hit

print gohitdb2

data = []
for gene in gnli:
	row = []
	for go in goli:
		check = [gene,go]
		if check in gohitdb2:
			row.append(1)
		else:
			row.append(0)
	data.append(row)

data = np.array(data)


#cluster genes by similarity? and GO terms too?

#make heatmap
#data = np.random.rand(10,4)
#print data

#data = np.array([[3,4],[5,1],[9,0]])
print data


#specify plot
fig, ax = plt.subplots()
heatmap = ax.pcolor(data, cmap=plt.cm.YlOrBr)
ax.set_xticks(np.arange(data.shape[1])+0.5,minor=False)
ax.set_yticks(np.arange(data.shape[0])+0.5,minor=False)
ax.invert_yaxis()
ax.xaxis.tick_top()

ax.set_xticklabels(goli, minor=False, rotation='vertical')
ax.set_yticklabels(gnli, minor=False)

plt.show()

