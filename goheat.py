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

datafolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/singlemotif/"
errormess = "USAGE: goheat.py motif_incl/excl (e.g.: goheat.py 2_12_4_excl)\nThe program will automatically find the correct input files in %s" %datafolder

if len(sys.argv) <= 1:
	sys.exit(errormess)

genefile = datafolder + sys.argv[1] + ".csv"
gofile = datafolder + sys.argv[1] + "_results.csv"

if not path.exists(genefile):
	sys.exit("ABORT: one of the input files was not found.\n%s" %errormess)
if not path.exists(gofile):
	sys.exit("ABORT: one of the input files was not found.\n%s" %errormess)

genes = csv.reader(open(genefile))
goterms = csv.reader(open(gofile))

#collect genes and GO terms
geneli, goli = [],[]
for line in genes:
	geneli.append(line[0:3])

for line in goterms:
	goli.append(line[0:2])

print geneli,goli


#from the GO-file: get data on which genes are associated with which GO terms
#make a table (gene x GO)

#cluster genes by similarity? and GO terms too?

#make heatmap
