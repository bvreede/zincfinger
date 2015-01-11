'''
This script can be used to detect distinct C2H2 zinc finger motifs
in a protein sequence.
The output is:
(1) a fasta file with domains in order;
(2) a csv file with all info that can be used for visualization.
(3) a stats file with numbers of motifs and the amount of double motifs, and a heatmap
(4) fasta files for all motifs with sequences that were found
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 10 October 2014
'''
import re
import pylab as pl
import numpy as np

### SPECIFY INFORMATION: DATA TO USE###
### input consists of: dbfolder/seqfolder/prefix-species_suffix
### output consists of: dbfolder/resfolder/
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
seqfolder = "sequences"
resfolder = "results"
prefix = "150111-SM00355"
suffix = "seq.fasta"
species = ["dmel","tcas","dpul","isca","smar","turt"]
motiflist = ['2_8_3','2_8_4','2_8_5','2_8_6','4_8_3','4_8_4','4_8_5','4_8_6','2_12_3','2_12_4','2_12_5','2_12_6','4_12_3','4_12_4','4_12_5','4_12_6','2_15_3','2_15_4','2_15_5','2_15_6','4_15_3','4_15_4','4_15_5','4_15_6'] # use only numbers, indicating the distances between C-C-H-H (separated by _)
#HFresidues = ['V','I','L','M','F','W','C','A','Y','H','T','S','P','G','R','K'] #hydrophobic residues.


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
alphabet = "ABCDEFGHIJKLMNPQRSTUVWXY"
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
for m in motiflist:
	mfile = open("%s/%s/%s.txt" %(dbfolder,resfolder,m), "w")
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
			for m in mots[n]:
				motifdoublecount[m] += 1
				if mots[n].count(m) > 1:
					combodict[frozenset([m])] += (1./mots[n].count(m)) #divide by total count, otherwise it will count +2 or more if it passes this point twice
		else:
			transl += translationdict[mots[n][0]]
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
			mfile = open("%s/%s/%s.txt" %(dbfolder,resfolder,m), "a")
			domain = motifdict[m]
			CC,CH,HH = m.split('_')
			for i in domain.finditer(fastadict[key]):
				motifcount[m] += 1 # count the found motif
				mseq = i.group() # the sequence picked up by the RE
				mfile.write(">%s\n%s\n\n" %(key,mseq))
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
		transl = translation(posmatrix,motdict,seqdict)
		outfasta.write(">%s\n%s\n\n" %(key,transl))
	outfasta.close()
	outputdb.close()

'''
make the stats!
in rows: for each motif, count the times it coincides with another motif, and divide this by the total of hits
for that motif.
'''
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

'''
Make a heatmap with the above stats.
'''
data = pl.array(doublematrix)
colourformap = "YlOrBr"
fig,ax = pl.subplots()
heatmap = pl.pcolor(data, cmap=colourformap)
cbar = pl.colorbar(heatmap)

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

#pl.show()
pl.savefig("%s/%s/heatmap_%s.svg" %(dbfolder,resfolder,colourformap), dpi = 300)

