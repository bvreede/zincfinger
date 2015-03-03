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
import csv, sys, urllib2
from os import path
import matplotlib.pyplot as plt
import numpy as np

#specify folder for inputfiles, and an errormessage given on every usage abort to the user
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/"
datafolder = "%sresults/singlemotif/" %dbfolder
#gosource = "%sdatabases/140720-SM00355-dmel2.csv" %dbfolder
gosource = "%sdatabases/150219-SM00355-dmel_corr.csv" %dbfolder
errormess = "USAGE: goheat.py motif_incl/excl (e.g.: goheat.py 2_12_4_excl)\nThe program will automatically find the correct input files in %s" %datafolder

#two parts of the AMIGO link, they are separated by the GO term.
golink1 = "http://golr.geneontology.org/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=10000&start=0&fl=bioentity&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=|&fq=document_category:%22bioentity%22&fq=taxon_closure_label:%22Drosophila%20%3Cfruit%20fly,%20genus%3E%22&facet.field=source&facet.field=type&facet.field=panther_family_label&facet.field=taxon_closure_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q="
golink2 = "&qf=bioentity^2&qf=bioentity_label_searchable^2&qf=bioentity_name_searchable^1&qf=bioentity_internal_id^1&qf=synonym^1&qf=isa_partof_closure_label_searchable^1&qf=regulates_closure^1&qf=regulates_closure_label_searchable^1&qf=panther_family_searchable^1&qf=panther_family_label_searchable^1&qf=taxon_closure_label_searchable^1"


#use AMIGO or BioMart inputfile: a for AMIGO, b for BioMart
#source = 'a'
source = 'b'


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
gohitdb = []
for n,line in enumerate(godb):
	if line[5] == 'biological_process' and line[6] in goli:
		print line
	if line[2] in pidli and line[6] in goli:
		gohit = [line[1],line[6]]
		#print line
		gohitdb.append(gohit)

#make an array for all genes and GO terms: 0 if they don't associate; 1 if they do
data = []
for gene in gnli:
	row = []
	for go in goli:
		check = [gene,go]
		if check in gohitdb:
			row.append(1)
		else:
			row.append(0)
	data.append(row)
data = np.array(data)


#cluster genes by similarity? and GO terms too?

#make a heatmap:
#specify plot
fig, ax = plt.subplots()
heatmap = ax.pcolor(data, cmap=plt.cm.YlOrBr)
#put labels halfway each column/row
ax.set_xticks(np.arange(data.shape[1])+0.5,minor=False)
ax.set_yticks(np.arange(data.shape[0])+0.5,minor=False)
plt.axis('tight') #remove the white bar
ax.invert_yaxis() #start from the top
ax.xaxis.tick_top() #labels on top

#set the labels
ax.set_xticklabels(goli, minor=False, rotation='vertical')
ax.set_yticklabels(gnli, minor=False)

plt.tight_layout()#prevents axis labels from being cut off

#show or save
plt.show()
# plt.savefig("%s%s_heatmap.png" %(datafolder,sys.argv[1]), dpi=300)

