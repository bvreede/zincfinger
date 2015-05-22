#!/usr/bin/python

'''
This script can be used to detect distinct C2H2 zinc finger motifs
in files with protein sequences.
The output is:
(1) a fasta file with domains in order;
(2) a csv file with all info that can be used for visualization.
(3) a stats file with numbers of motifs and the amount of double motifs, and a heatmap of these stats
(4) fasta files for all motifs with sequences that were found
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 10 October 2014
'''
import re, sys
import pylab as pl
import numpy as np

if len(sys.argv) <= 1:
	sys.exit("USAGE: findmotif.py path/to/inputfile (input needs to be a protein fasta file).\nOutputfolders are indicated in the script; edit the script if you want to alter them.")

infile = sys.argv[1]


### SPECIFY INFORMATION: DATA TO USE###
### output consists of: dbfolder/resfolder/ or dbfolder/imgfolder or dbfolder/seqfolder, depending on output type
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
seqfolder = "sequences"
resfolder = "results"
imgfolder = "images"
dbfolder2 = "databases"

#define motifs: first a complete dataset
moCC = [2,3,4] #distances between CC
moCH = [7,8,9,10,11,12,13,14,15] #distances between CH
moHH = [2,3,4,5,6] #distances between HH
# turn this on if you want to search for all possible combinations of the above: (and don't forget to turn the custom list off!)
#motiflist = ['%s_%s_%s' %(m,n,o) for m in moCC for n in moCH for o in moHH] 

# turn this on for custom motif list
motiflist = ['2_7_4','2_8_3','2_9_3','2_10_5','2_11_3','2_11_4','2_12_2','2_12_3','2_12_4','2_12_5','2_12_6','2_13_3',
'2_13_4','2_14_3','2_14_4','2_15_4','3_8_3','4_12_3','4_12_4','4_15_3']

# options: what do you want the script to do?
stats = 0 # set to 1 if you want to make stats (how many motifs were found; how many duplicates; etc) and images (heatmap)
saveseq = 0 # set to 1 if you want to make fasta files of all motifs that were found
translate_hits = 1 #set to 1 if you want to generate an output fasta file with the translated motifhits
frequency = 1 #set to 1 if you want to generate the list of all motifs used for frequency-dependent motif sampling (NB only works when translate_hits is also set to 1!

# output file with list of all motifs, to be used for frequency-dependent sampling
if frequency == 1:
	allmotifstxt = open("%s/%s/allmotifs.txt" %(dbfolder,dbfolder2), "w")

'''
Define regular expressions for all zf-domains (by the C-H distances), and save in a
dictionary. Each domain is only represented as a forward domain.
'''
motifdict = {} # dictionary of motifs and their regular expression
motiflength = {} # dictionary of motifs and their lengths
plink = 4
alink = 4
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
alphabet = """ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`"""
translationdict = {motif: alphabet[a] for a,motif in enumerate(motiflist)} # dictionary of motifs and the corresponding string element

'''
Make dictionaries to count all motifs and combinations of motifs.
'''
motifcount = {m: 0 for m in motiflist}
motifdoublecount = {m: 0 for m in motiflist}
combodict = {a: 0 for a in set([frozenset([m,n]) for m in motiflist for n in motiflist])}

'''
Open files to save the sequences of all hits per motif
'''
if saveseq == 1:
	for m in motiflist:
		mfile = open("%s/%s/%s.txt" %(dbfolder,seqfolder,m), "w")
		mfile.close()

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
translate a set of motifs + corresponding start sites and lengths
to a single string that contains both translations for the domains
(or domain combinations) used, and for the spaces in between.
'Translation' in this case means: into the string used for clustering.
'''
def translation(posmatrix,motdict,seqdict):
	# translate posmatrix and motdict into two lists of (lists of) start sites and the corresponding motifs 
	starts, mots,seqs = [],[],[]
	for pos in posmatrix:
		st,mo,sq = [],[],[]
		for p in pos:
			mos = motdict[p].split('/')
			seq = seqdict[p].split('/')
			for i,m in enumerate(mos):
				mo.append(m)
				st.append(p)
				sq.append(seq[i])
		seqs.append(sq)
		starts.append(st)
		mots.append(mo)
	# translate the starts + motifs into a string
	transl = ''
	for n,start in enumerate(starts):
		#add a spacer if the motifs are not connected
		if n != 0:
			ends = [s + motiflength[mots[n-1][k]] for k,s in enumerate(starts[n-1])]
			if min(start) - max(ends) > 9: # 9 is the cut off for domains in a loop
				transl += 'Z'			
		#translate the domains to either a single letter or a regex if there are more
		if len(mots[n]) > 1:
			regex = ''
			for m in mots[n]:
				regex += translationdict[m] + '|'
			allmots = set([frozenset([q,r]) for q in mots[n] for r in mots[n] if len(frozenset([q,r]))>1])
			for a in allmots:
				combodict[a] += 1
			transl += '{' + regex[:-1] + '}'
			allmotifstxt.write('{' + regex[:-1] + '}\n') #add the regex to the 'allmotifs' document for frequency-dependent sampling
			for m in mots[n]:
				motifdoublecount[m] += 1
				if mots[n].count(m) > 1:
					combodict[frozenset([m])] += (1./mots[n].count(m)) #divide by total count, otherwise it will count +2 or more if it passes this point twice
		else:
			transl += translationdict[mots[n][0]]
			allmotifstxt.write("%s\n" %translationdict[mots[n][0]]) #add the (translated) motif to the 'allmotifs' document for frequency-dependent sampling
	return transl
	
'''
Read the fasta file, open an output file, and scan all sequences for the presence of
motifs.
'''
# write the header for the database
sp = infile.split('/')[-1].split('.')[0]
fastadb = open("%s" %(infile))
outputdb = open("%s/%s/%s_motifhits.csv" %(dbfolder,resfolder,sp), "w")
allmotifsfa = open("%s/%s/allmotifs.fasta" %(dbfolder,seqfolder), "w")

if translate_hits == 1:
	outfasta = open("%s/%s/%s_motifseq.fasta" %(dbfolder,resfolder,sp), "w")
outputdb.write("Gene_stable_ID,Gene_name,Protein_stable_ID,Sequence_length,")
for m in motiflist:
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
		thisseqcount = 0 #per motif per seq, to give an index for each aminoacid sequence found
		if saveseq == 1:
			mfile = open("%s/%s/%s.txt" %(dbfolder,seqfolder,m), "a")
		domain = motifdict[m]
		CC,CH,HH = m.split('_')
		for i in domain.finditer(fastadict[key]):
			motifcount[m] += 1 # count the found motif
			mseq = i.group() # the sequence picked up by the RE
			if saveseq == 1:
				mfile.write(">%s\n%s\n\n" %(key,mseq))
			allmotifsfa.write(">%s|%s-%s\n%s\n\n" %(key,translationdict[m],thisseqcount,mseq))
			thisseqcount += 1
			strt = i.start() + plink
			if strt in seqdict:
				ns = seqdict[strt] + "/" + mseq
				nm = motdict[strt] + "/" + m
				seqdict[strt] = ns
				motdict[strt] = nm
			else:
				seqdict[strt] = mseq
				motdict[strt] = m
				poslist.append(strt)
		if saveseq == 1:
			mfile.close()

	# screen to filter the domains that are conflicting.
	poslist = list(set(poslist))
	poslist.sort() # sorts all starting positions to make sure they are in order
	posdone = [] # to add positions that have already been assessed
	posmatrix = [] # to collect all sequences that should be assessed together
	pos4matrix = [] # to collect the lines of the matrix
	for pos in poslist:
		if pos in posdone:
			continue
		pos4matrix.append(pos)
		for p in range(1,6):
			pp = pos + p
			if pp in poslist:
				pos4matrix.append(pp)
				posdone.append(pp)
		posmatrix.append(pos4matrix)
		pos4matrix = []

	# two different lists: one accurate, for visualization
	for ml in motiflist: # the list of motifs as stated above (all motifs to look for)
		outstring = ""
		for pos in motdict:
			posmots = motdict[pos].split('/')
			if len(posmots) > 1:
				for pm in posmots:
					if pm == ml:
						outstring += str(pos) + '|'
			else:
				if motdict[pos] == ml:
					outstring += str(pos) + '|'
		outputdb.write("%s," %outstring[:-1])			
	outputdb.write("\n")

	# another for clustering: space or no space
	if translate_hits == 1:
		transl = translation(posmatrix,motdict,seqdict)
		if len(transl) > 0:
			outfasta.write(">%s\n%s\n\n" %(key,transl))
if translate_hits == 1:
	outfasta.close()
outputdb.close()

'''
make the stats!
in rows: for each motif, count the times it coincides with another motif, and divide this by the total of hits
for that motif.
Make a heatmap with the stats.
'''
if stats == 1:
	#open file and array
	stats = open("%s/%s/motifstats_total.csv" %(dbfolder,resfolder), "w")
	doublematrix = []
	#print headers
	stats.write(",")
	for m in motiflist:
		stats.write("%s," %m)
	stats.write("total_combo,total\n")
	#per motif count combinations
	for m in motiflist:
		stats.write("%s," %m)
		doublelist = []
		for n in motiflist:
			if motifcount[m] == 0:
				norm_combo = 0
			else:
				norm_combo = combodict[frozenset([m,n])]/float(motifcount[m])
			stats.write("%s," %norm_combo)
			doublelist.append(norm_combo)
		doublematrix.append(doublelist)
		stats.write("%s,%s\n" %(motifdoublecount[m],motifcount[m]))
	stats.close()
	
	###HEATMAP###
	data = pl.array(doublematrix)
	colourformap = "YlOrBr"
	fig,ax = pl.subplots()
	heatmap = pl.pcolor(data, cmap=colourformap)
	cbar = pl.colorbar(heatmap)
	
	# following commented code from stackoverflow; probably not necessary, but leaving it just in case
	#cbar.ax.set_yticklabels(['0','1','2','>3'])
	#cbar.set_label('#double motifs / total motifs', rotation=270)
	
	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
	pl.axis('tight') #remove the white bar
	ax.invert_yaxis() #make sure it starts counting from the top
	
	#make the labels
	ax.set_xticklabels(motiflist, minor=False, rotation=90)
	ax.set_yticklabels(motiflist, minor=False)
	
	# save the figure
	pl.savefig("%s/%s/heatmap_%s.svg" %(dbfolder,imgfolder,colourformap), dpi = 300)
