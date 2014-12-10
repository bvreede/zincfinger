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
species = ["dmel","tcas","dpul","isca","smar"]
motiflist = ['2_8_3','2_8_4','2_8_5','4_8_3','4_8_4','4_8_5','2_12_3','2_12_4','2_12_5','4_12_3','4_12_4','4_12_5','2_15_3','2_15_4','2_15_5','4_15_3','4_15_4','4_15_5'] # use only numbers, indicating the distances between C-C-H-H (separated by _)

HFresidues = ['V','I','L','M','F','W','C','A','Y','H','T','S','P','G','R','K'] #hydrophobic residues. Not used in the script, as there does not seem to be a reliable way to check for false positives with this information. 



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
plink = 7
alink = 7
for m in motiflist:
	cc,ch,hh = m.split('_')
	motif = '[A-Z]{%s}C[A-Z]{%s}C[A-Z]{%s}H[A-Z]{%s}H[A-Z]{%s}' %(plink,cc,ch,hh,alink) # construct the regular expression
	remotif = re.compile(motif)
	motifdict[m] = remotif
	d = [int(x) for x in cc,ch,hh]
	dl = sum(d) + 4
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
translate a set of possible zf hits at a range of positions to a single
list of positions + the correct motif at that position.
Yeah, I know. Not possible.
'''
def resolvemotifs(sequences,motifs):
	#score the motifs on how good they are: L on H-2, other hydrophobic on L-5, R in final H-H loop
	scores = []
	for n in range(len(sequences)):
		score = 0
		CC,CH,HH = motifs[n].split('_')
		Lpos = plink+int(CC)+int(CH)-1
		if sequences[n][Lpos] == 'L':
			score += 2
		elif sequences[n][Lpos] in HFresidues:
			score += 1
		Fpos = Lpos-6
		if sequences[n][Fpos] == 'F':
			score += 2
		elif sequences[n][Fpos] in HFresidues:
			score += 1
		Hs = plink+int(CC)+int(CH)+4
		He = plink+int(CC)+int(CH)+int(HH)+2
		H2H = sequences[n][Hs:He]
		for r in H2H:
			if r == 'R':
				score += 1
		scores.append(score)
	return scores
		

def resolvematrix(posmatrix,seqdict,motdict):
	finalpos = []
	finaldict = {}
	for line in posmatrix: # lines of positions that need to be assessed together
		allpos,allmots,allseqs,allend = [],[],[],[]
		for item in line:
			mots = motdict[item].split('/') # get all motifs on this position
			seqs = seqdict[item].split('/') # get all sequences on this position
			for n in range(len(mots)):
				allpos.append(item)
				allmots.append(mots[n])
				allseqs.append(seqs[n])
		allscores = resolvemotifs(allseqs,allmots)
		for n in range(len(allpos)):
			end = int(allpos[n])+int(motiflength[allmots[n]])+7
			allend.append(end)
		if len(allpos) > 1:
			print allpos,allmots,allseqs,allscores,allend
	return finalpos,finaldict


'''
Per species, read the fasta file, open an output file, and scan all sequences for the presence of
motifs.
'''
for sp in species:
	# write the header for the database
	fastadb = open("%s/%s/%s-%s_%s" %(dbfolder,seqfolder,prefix,sp,suffix))
	outputdb = open("%s/%s/motifhits_%s.csv" %(dbfolder,resfolder,sp), "w")
	outfasta = open("%s/%s/motifseq_%s.fasta" %(dbfolder,resfolder,sp), "w")
	outputdb.write("Gene_stable_ID,Gene_name,Protein_stable_ID,Sequence_length,")
	for m in motifdict:
		outputdb.write("%s," %m)
	outputdb.write("\n")

	# read the fasta file into memory
	fastadict = fastadicter(fastadb) # translate the fasta file into a dictionary
	seqdict = {} # dictionary for [start position]: sequence
	motdict = {} # dictionary for [start position]: motif type

	# go through the sequences
	for key in fastadict:
		print key
		# get info for the first columns (ID and sequence length)
		ids = key.split('|')
		seqlen = len(fastadict[key]) #length of the sequence
		outputdb.write("%s,%s,%s,%s," %(ids[0],ids[1],ids[2],seqlen)) #turn the header name into gene ID/name/prot ID

		# screen the sequence for motifs
		seqdict.clear() #for the fasta file: collect positions as key and the corresponding sequence at that position as value
		motdict.clear() #as seqdict, but with the motif name instead of sequence
		poslist = [] #for the fasta file: collect all positions to put them in order later on

		# check each motif individually
		for m in motifdict: #go through each motif and find all instances in the sequence
			domain = motifdict[m]
			positions = ''
			CC,CH,HH = m.split('_')
			for i in domain.finditer(fastadict[key]):
				mseq = i.group() # the sequence picked up by the RE
				
				LL = mseq[plink+int(CC)+int(CH)-1]
				FF = mseq[plink+int(CC)+int(CH)-1-6]
				#if LL not in HFresidues or FF not in HFresidues:
				#	print LL, FF
				#print mseq[plink+int(CC)+int(CH)+2:plink+int(CC)+int(CH)+int(HH)+4]
				strt = i.start()
				positions += str(strt)
				positions += '|'
				if strt in seqdict:
					ns = seqdict[strt] + "/" + mseq
					nm = motdict[strt] + "/" + m
					#motiflist.append(nm) # possibly superfluous if no longer working with double motifs for cluster
					seqdict[strt] = ns
					motdict[strt] = nm
				else:
					seqdict[strt] = mseq
					motdict[strt] = m
					poslist.append(strt)
			# consider moving the following 3 lines because they need to be after the motif screen
			#outputdb.write('%s,' %(positions[:-1])) #remove final pipe from total positions
		#outputdb.write("\n")
		#outfasta.write(">%s\n" %key) #start collecting the info for the fasta file

		# after screen, filter the domains that are conflicting.
		poslist = list(set(poslist))
		poslist.sort() #sorts all starting positions to make sure they are in order

		#cl = 0 
		#readmetxt,domaindict = makereadme(motiflist)

		posdone = [] # to add positions that have already been assessed
		posmatrix = [] # to collect all sequences that should be assessed together
		pos4matrix = [] # to collect the lines of the matrix
		for n in range(len(poslist)):
			pos = poslist[n]
			if pos in posdone:
				continue
	
			pos4matrix.append(pos)
			for p in range(1,6):
				pp = p + pos
				if pp in poslist:
					pos4matrix.append(pp)
					posdone.append(pp)
			posmatrix.append(pos4matrix)
			pos4matrix = []
	

		# now I have a posmatrix of all positions for this gene that need to be assessed together
		# to generate A SINGLE OUTPUT of position - key.
		# how to fucking do this?
		
		finalpos,finaldict = resolvematrix(posmatrix,seqdict,motdict)




"""

			# determine here whether there is a chain or not. Perhaps have a parameter that says whether the motif is part of a chain, and if so: is it the first, one of the middle, or last? if it's the first/last: put a 'space' sign in the translation!

			# for those that conflict: determine the right one or if there is a possible ambiguous zf:
			# does it have the two hydrophobic residues in between C and H?
			# if it is in a chain, are there 7 residues between the last H and the next C? [linker sequence]
			# if these qualities have been met, and there is still a conflict: then label the zf with the appropriate combo term. If there is a triple or more, report it for manual curation (perhaps make a file with exceptions so that once manually curated zfs have been identified and assessed, they won't need to be done again?)
"""


