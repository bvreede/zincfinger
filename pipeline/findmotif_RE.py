'''
This script can be used to detect distinct C2H2 zinc finger motifs
in a protein sequence.
The output is:
(1) a fasta file with domains in order (+ a readme);
(2) a csv file with all info that can be used for visualization;
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 10 October 2014
'''
import re

customout = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/playground/2_8_3.txt", "w")

### SPECIFY INFORMATION: DATA TO USE###
### input consists of: dbfolder/seqfolder/prefix-species_suffix
### output consists of: dbfolder/resfolder/
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
seqfolder = "sequences"
resfolder = "results"
prefix = "140720-SM00355"
suffix = "seq.fasta"
species = ["allz"] #["dmel","tcas","dpul","isca","smar"]
motiflist = ['2_8_3','2_8_4','2_8_5','4_8_3','4_8_4','4_8_5','2_12_3','2_12_4','2_12_5','4_12_3','4_12_4','4_12_5','2_15_3','2_15_4','2_15_5','4_15_3','4_15_4','4_15_5'] # use only numbers, indicating the distances between C-C-H-H (separated by _)

### END OF CUSTOM INFO. DON'T FUCK WITH THE FOLLOWING ###

'''
Motif to cluster-able sequence: make a dictionary to translate
the zf-domains to a single letter
'''
alphabet = "ABCDEFGHIJKLMNPQRSTUVWXYZ"
domaindict = {}
for i in range(len(motiflist)):
		domaindict[motiflist[i]] = alphabet[i]

'''
Define regular expressions for all zf-domains (by the C-H distances), and save in a
dictionary. Each domain is only represented as a forward domain.
'''
motifdict = {}
motiflength = {}
for m in motiflist:
	cc,ch,hh = m.split('_')
	motif = '[A-Z]{2}C[A-Z]{%s}C[A-Z]{%s}H[A-Z]{%s}H' %(cc,ch,hh) # construct the regular expression
	remotif = re.compile(motif)
	motifdict[m] = remotif
	d = [int(x) for x in cc,ch,hh]
	dl = sum(d)
	motiflength[m] = dl

'''
From an open fasta file, generate a dictionary containing the header as key, and the
sequence as value.
'''
def fastadicter(fastadb):
	fastadict = {}
	sequence = ""
	header = ""
	for line in fastadb:
		if line[0] == ">":
			if header != "":
				fastadict[header] = sequence
			header = line[1:].strip()
			sequence = ""
		else:
			sequence += line.strip()
	fastadict[header] = sequence
	return fastadict



'''
INFO FOR THE FASTA README! Collects all motifs and combinations of ~
to give a 'legend' to the information
'''
def makereadme(motiflist):
	domaindict = {"space (undefined length)": 'O'}
	domains = "ABCDEFGHIJKLMNPQRSTUVWXYZabcdefghijklmnpqrstuvwxyz"
	motiflist = list(set(motiflist))
	for i in range(len(motiflist)):
		domaindict[motiflist[i]] = domains[i]
	readmetxt = "The attached fasta file contains the sequence of different C2H2 zinc finger \
	domains in the provided source file. They are encoded as follows:\n"
	for k in domaindict:
		rm = "%s = %s\n" %(domaindict[k],k)
		readmetxt += rm
	return readmetxt,domaindict

'''
Per species, read the fasta file, open an output file, and scan all sequences for the presence of
motifs.
'''
cntr = 0
for sp in species:
	# write the header for the database
	fastadb = open("%s/%s/%s-%s_%s" %(dbfolder,seqfolder,prefix,sp,suffix))
	outputdb = open("%s/%s/motifhits_%s.csv" %(dbfolder,resfolder,sp), "w")
	outfasta = open("%s/%s/motifseq_%s.fasta" %(dbfolder,resfolder,sp), "w")
	outputdb.write("Gene_stable_ID,Gene_name,Protein_stable_ID,Sequence_length,")
	for m in motifdict:
		outputdb.write("%s," %m)
	outputdb.write("\n")
	# database content follows below
	fastadict = fastadicter(fastadb) # translate the fasta file into a dictionary
	seqdict = {}
	for key in fastadict:
		# get info for the first columns (ID and sequence length)
		ids = key.split('|')
		seqlen = len(fastadict[key]) #length of the sequence
		outputdb.write("%s,%s,%s,%s," %(ids[0],ids[1],ids[2],seqlen)) #turn the header name into gene ID/name/prot ID

		# screen the sequence for motifs
		seqdict.clear() #for the fasta file: collect positions as key and the corresponding domain at that position as value
		seqlist = [] #for the fasta file: collect all positions to put them in order later on
		for m in motifdict: #go through each motif and find all instances in the sequence
			domain = motifdict[m]
			positions = ''
			for i in domain.finditer(fastadict[key]):
				strt = i.start()
				positions += str(strt)
				positions += '|'
				customseq = i.group()
				if m == '2_8_3':
					cntr +=1
					customout.write(">seq%s\n%s\n\n" %(cntr,customseq))
				if strt in seqdict:
					nm = seqdict[strt] + "/" + m
					motiflist.append(nm)
					seqdict[strt] = nm
				else:
					seqdict[strt] = m
					seqlist.append(strt)
			outputdb.write('%s,' %(positions[:-1])) #remove final pipe from total positions
		outputdb.write("\n")
		outfasta.write(">%s\n" %key) #start collecting the info for the fasta file
		seqlist.sort()
		cl = 0
		readmetxt,domaindict = makereadme(motiflist)
		for n in seqlist:
			if cl != 0 and int(n) - cl > 11:
				outfasta.write("O")
			if len(seqdict[n].split('/')) > 1: #if we're dealing with a double or triple motif
				lmot = seqdict[n].split('/')
				lencol = []
				for lm in lmot:
					lencol.append(motiflength[lm])
				l3 = max(lencol)
				cl = int(n) + l3
			elif cl < int(n) + motiflength[seqdict[n]]:
				cl = int(n) + motiflength[seqdict[n]]
			outfasta.write(domaindict[seqdict[n]])
		outfasta.write("\n\n")
	#MAKE README
	readme = open("%s/%s/motifseq_%s-README.txt" %(dbfolder,resfolder,sp), "w")
	readme.write(readmetxt)
	readme.close()
	#END README
	outputdb.close()
	outfasta.close()
	#matrix(sp)



"""
'''
Add-on made for clustering:
Turn the csv file into a matrix with single genes in the rows, and a motif count
in the columns. Only takes into account the first protein per gene.
'''
def matrix(sp):
	prev_output = open("%s/%s/motifhits_%s.csv" %(dbfolder,resfolder,sp)) # outputfile that saves motif locations and sequence length
	new_output = open("%s/%s/motifhits_%s-count.csv" %(dbfolder,resfolder,sp), "w") #outputfile that saves n motifs
	testID = "Gene_stable_ID"
	for line in prev_output:
		l = line.strip().split(',')
		if l[0] == testID:
			if l[0] == "Gene_stable_ID":
				lr = line.replace('Sequence_length,','')
				new_output.write(lr)
			continue
		else:
			new_output.write("%s,%s,%s," %(l[0],l[1],l[2]))
			for i in range(4,len(l)-1):
				if len(l[i]) != 0:
					new_output.write(str(len(l[i].split('|'))))
					new_output.write(",")
				else:
					new_output.write("0,")
			new_output.write("\n")
"""

