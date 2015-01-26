#!/usr/bin/python
'''
This script can be used to find and group all genes with either
(1) domain X exclusively; or
(2) containing domain X, irrespective of others.
The input required for this script is the output database of findmotif.py,
consisting of three header columns (gene ID, gene name, protein ID), a column
for sequence lengths, and per motif the locations of hits (each motif in a separate
column).
Possible option for improvement: add an escape to the 'exclusive' hit finder that
allows the existence of other motifs so long as they overlap.

Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 25 January 2015
'''
import csv,sys

if len(sys.argv) <= 1:
	sys.exit("USAGE: select4gorilla.py path/to/inputfile \n(inputfile is database of motifhits -- output of findmotif.py -- in csv)")
inputdb = csv.reader(open(sys.argv[1])) # input file
outlist = sys.argv[1].split('/')[:-1]
outfolder = '/'.join(outlist) + '/singlemotif' # all folders from input path minus the file; used as folder for any output

#read input and put in memory
db = []
for k,line in enumerate(inputdb):
	if k == 0:
		header = line
	else:
		db.append(line)

'''
Go through database with a given column in mind
and print the outputfiles for this column.
'''
def readdb(n,outi,oute):
	for line in db:
		pline = ','.join(line) #the line as it will be on a csv resultsfile
		if line[n] == '': #there is a hit on this protein for the motif
			# (1) motif inclusive:
			# put those lines in the database that have hits on this motif.
			outi.write("%s\n" %pline)
			#check if there are other motif hits in this protein
			e = 0
			for k in range(4,len(line)):
				if line[k] != '':
					if k == n:
						continue
					else:
						e = 1 #other hits found; no longer interested in this motif
						break
			# (2) motif exclusive:
			# put those lines in the database that have hits on this motif, and nothing in others
			if e == 0:
				oute.write("%s\n" %pline)

#per motif: open outputfile, search through the database
for n,m in enumerate(header):
	if n > 3 and len(m) > 0: #these are the motifs
		outi = open("%s/%s_incl.csv" %(outfolder,m), 'w')
		oute = open("%s/%s_excl.csv" %(outfolder,m), 'w')
		readdb(n,outi,oute)
		outi.close()
		oute.close()

