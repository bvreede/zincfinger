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

HFresidues = ['V','I','L','M','F','W','C','A','Y','H','T','S','P','G','R','K'] #hydrophobic residues.

domains = set() #collect all different domains and combinations of domains that are *actually* found & used.

'''
Define regular expressions for all zf-domains (by the C-H distances), and save in a
dictionary. Each domain is only represented as a forward domain.
'''
motifdict = {} # dictionary of motifs and their regular expression
motiflength = {} # dictionary of motifs and their lengths
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
Motif to cluster-able sequence: make a dictionary to translate
the zf-domains to a single letter
'''
alphabet = "ABCDEFGHIJKLMNPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz!@^&*()-=+{}:;|<>?/,."
translationdict = {} # dictionary of motifs and the corresponding string element
alphacount = 0

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
translate a set of possible zf hits at a range of positions to a single
list of positions + the correct motif at that position (or the best possible double).
'''
def resolvemotifs(positions,sequences,motifs,prevend,nextst): # called PER (CONFLICTING SET OF) MOTIFS
	#score the motifs on how good they are: L on H-2, other hydrophobic on L-5, R in final H-H loop
	scores = []
	for n,seq in enumerate(sequences):
		score = 0
		CC,CH,HH = motifs[n].split('_')
		#check the presence of hydrophobic residues in specific locations
		Lpos = plink+int(CC)+int(CH)-1
		if seq[Lpos] == 'L':
			score += 2
		elif seq[Lpos] in HFresidues:
			score += 1
		Fpos = Lpos-6
		if seq[Fpos] == 'F':
			score += 2
		elif seq[Fpos] in HFresidues:
			score += 1
		#check the presence of an R in the H-H loop
		Hs = plink+int(CC)+int(CH)+4
		He = plink+int(CC)+int(CH)+int(HH)+2
		H2H = seq[Hs:He]
		for r in H2H:
			if r == 'R':
				score += 1
				break #to prevent double scoring in event of RR
		#check the location respective to previous and next domain
		st_motif = int(positions[n])
		end_motif = st_motif + int(motiflength[motifs[n]]) + 7
		if st_motif < min(prevend):
			score -= 5
		elif st_motif in prevend:
			score += 5
		if end_motif > max(nextst):
			score -= 5
		elif end_motif in nextst:
			score += 5
		scores.append(score)
	return scores
		

def resolvematrix(posmatrix,seqdict,motdict): #called PER GENE
	allstarts,allmotifs,alllength = [],[],[]
	prevend = [0]
	for k,line in enumerate(posmatrix): # lines of positions that need to be assessed together
		allpos,allmots,allseqs = [],[],[]
		finalpos,finalmot,finalseq = [],[],[]
		for item in line:
			mots = motdict[item].split('/') # get all motifs on this position
			seqs = seqdict[item].split('/') # get all sequences on this position
			for n,mot in enumerate(mots):
				allpos.append(item)
				allmots.append(mot)
				allseqs.append(seqs[n])
		try:
			allscores = resolvemotifs(allpos,allseqs,allmots,prevend,posmatrix[k+1])
		except IndexError:
			allscores = resolvemotifs(allpos,allseqs,allmots,prevend,[0])
		# prevend: the previous end. Get the indices of the HIGHEST scoring motifs, and stuff the END SITE in a list
		prevend = []
		for i,score in enumerate(allscores):
			if score == max(allscores): # this/ese is/are our top hit(s) from the current round!
				prevend.append(allpos[i] + motiflength[allmots[i]] + 7)
				# pick the best scoring motif
				finalpos.append(allpos[i])
				finalmot.append(allmots[i])
				finalseq.append(allseqs[i])
		allstarts.append(min(allpos))
		alllength.append(max(allpos) + max([motiflength[i] for i in finalmot])-min(allpos)) #max length of the motif WITHOUT LINKER SEQ
		allmotifs.append(frozenset(finalmot))
		domains.add(frozenset(finalmot))
	return allstarts,allmotifs,alllength

'''
translate a set of motifs + corresponding start sites and lengths
to a single string that contains both translations for the domains
(or domain combinations) used, and for the spaces in between.
'''
def translation(starts,motifs,lengths):
	transl = ''
	for n,start in enumerate(starts):
		if n != 0:
			if start - (starts[n-1] + lengths[n-1] + 7) > 0:
				transl += '_'
		if motifs[n] not in translationdict:
			global alphacount
			translationdict[motifs[n]] = alphabet[alphacount]
			alphacount += 1
		transl += translationdict[motifs[n]]
	return transl


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
				strt = i.start()
				positions += str(strt)
				positions += '|'
				if strt in seqdict:
					ns = seqdict[strt] + "/" + mseq
					nm = motdict[strt] + "/" + m
					seqdict[strt] = ns
					motdict[strt] = nm
				else:
					seqdict[strt] = mseq
					motdict[strt] = m
					poslist.append(strt)

		# screen to filter the domains that are conflicting.
		poslist = list(set(poslist))
		poslist.sort() #sorts all starting positions to make sure they are in order
		posdone = [] # to add positions that have already been assessed
		posmatrix = [] # to collect all sequences that should be assessed together
		pos4matrix = [] # to collect the lines of the matrix
		for pos in poslist:
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
		starts,motifs,lengths = resolvematrix(posmatrix,seqdict,motdict)
		transl = translation(starts,motifs,lengths)
		outfasta.write(">%s\n%s\n\n" %(key,transl))
			#writeimage(key,starts,lengths)
			# write image interpretation file. Info needed:
				# key
				# per motif:	- start site (add plink to this!
				#		- length
	outfasta.close()
