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
import re, sys,config
import pylab as pl
import numpy as np


if len(sys.argv) <= 1:
	sys.exit("USAGE: findmotif.py path/to/inputfile (input needs to be a protein fasta file).\nOutputfolders are indicated in the script; edit the script if you want to alter them.")

infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]
#NB! filenames should always start with an input ID specifier (separate elements in dashes) and end with output ID specifiers.
#e.g.: 150525-dmel_seq.fa or 141212-tcas_heatmap.svg

"""
# options: what do you want the script to do?
stats = 1 # set to 1 if you want to make stats (how many motifs were found; how many duplicates; etc) and images (heatmap)
saveseq = 1 # set to 1 if you want to make fasta files of all motifs that were found
translate_hits = 1 #set to 1 if you want to generate an output fasta file with the translated motifhits
frequency = 1 #set to 1 if you want to generate the list of all motifs used for frequency-dependent motif sampling (NB only works when translate_hits is also set to 1!
"""

### OUTPUT FILES ###
### specify only names, and no extensions, so sub-outputfiles can be specified later ###
# images: bar graph and heatmap
bargraphfig = "%s/%s/%s_motifs-bar" %(config.mainfolder,config.imgfolder,infilebrev)
heatmapfig = "%s/%s/%s_motifs-heat" %(config.mainfolder,config.imgfolder,infilebrev)
# databases: hit results by site
hitsdb = "%s/%s/%s_hitsdb" %(config.mainfolder,config.dbfolder,infilebrev)
# fasta files: translated hits, aa sequence of hits
transdb = "%s/%s/%s_protstring" %(config.mainfolder,config.seqfolder,infilebrev)
motseq = "%s/%s/%s-motseq" %(config.mainfolder,config.seqfolder,infilebrev)
allmotifs = "%s/%s/%s_allmotifs" %(config.mainfolder,config.seqfolder,infilebrev) #for frequency and aa sequence of specific motifs
statsdb = "%s/%s/%s_motifstats" %(config.mainfolder,config.resfolder,infilebrev)


#Make dictionaries to count all motifs and combinations of motifs.
motifcount = {m: 0 for m in config.motiflist}
motifdoublecount = {m: 0 for m in config.motiflist}
combodict = {a: 0 for a in set([frozenset([m,n]) for m in config.motiflist for n in config.motiflist])}


def translation(posmatrix,motdict,seqdict):
	'''
	translate a set of motifs + corresponding start sites and lengths
	to a single string that contains both translations for the domains
	(or domain combinations) used, and for the spaces in between.
	'Translation' in this case means: into the string used for clustering.
	'''
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
			ends = [s + config.motiflength[mots[n-1][k]] for k,s in enumerate(starts[n-1])]
			if min(start) - max(ends) > 9: # 9 is the cut off for domains in a loop
				transl += 'Z'			
		#translate the domains to either a single letter or a regex if there are more
		if len(mots[n]) > 1:
			regex = ''
			for m in mots[n]:
				regex += config.translationdict[m] + '|'
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
			transl += config.translationdict[mots[n][0]]
			allmotifstxt.write("%s\n" %config.translationdict[mots[n][0]]) #add the (translated) motif to the 'allmotifs' document for frequency-dependent sampling
	return transl


def makeheatmap(doublematrix,name):
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
	ax.set_xticklabels(config.motiflist, minor=False, rotation=90)
	ax.set_yticklabels(config.motiflist, minor=False)
	
	# save the figure
	pl.savefig("%s-%s.svg" %(heatmapfig,name), dpi = 300)
	pl.clf()
	pl.close()


# Read the fasta file, open output files
fastadb = open("%s" %(infile))
outputdb = open("%s.csv" %hitsdb, "w")
outfasta = open("%s.fa" %transdb, "w")
allmotifsfa = open("%s.fa" %allmotifs, "w")
allmotifstxt = open("%s.txt" %allmotifs, "w")
for m in config.motiflist:
	mfile = open("%s-%s.fa" %(motseq,m), "w")
	mfile.close()
outfasta = open("%s.fa" %transdb, "w")



# Make the headers for the outputdb
outputdb.write("Gene_stable_ID,Gene_name,Protein_stable_ID,Sequence_length,")
for m in config.motiflist:
	outputdb.write("%s," %m)
outputdb.write("\n")


fastadict = config.fastadicter(fastadb) # translate the fasta file into a dictionary

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
	for m in config.motifdict: #go through each motif and find all instances in the sequence
		thisseqcount = 0 #per motif per seq, to give an index for each aminoacid sequence found
		mfile = open("%s-%s.fa" %(motseq,m), "a")
		domain = config.motifdict[m]
		CC,CH,HH = m.split('_')
		for i in domain.finditer(fastadict[key]):
			motifcount[m] += 1 # count the found motif
			mseq = i.group() # the sequence picked up by the RE
			mfile.write(">%s\n%s\n\n" %(key,mseq))
			allmotifsfa.write(">%s|%s-%s\n%s\n\n" %(key,config.translationdict[m],thisseqcount,mseq))
			thisseqcount += 1
			strt = i.start() + config.plink
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
	for ml in config.motiflist: # the list of motifs as stated above (all motifs to look for)
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
	if len(transl) > 0:
		outfasta.write(">%s\n%s\n\n" %(key,transl))
outfasta.close()
allmotifstxt.close()
outputdb.close()

'''
make the stats!
in rows: for each motif, count the times it coincides with another motif, and divide this by the total of hits
for that motif.
Make a heatmap with the stats.
'''
#open file and array
stats = open("%s.csv" %statsdb, "w")
doublematrix1,doublematrix2 = [],[]
#print headers
stats.write(",")
for m in config.motiflist:
	stats.write("%s," %m)
stats.write("total_combo,total\n")
#per motif count combinations
for m in config.motiflist:
	stats.write("%s," %m)
	doublelist1,doublelist2 = [],[]
	for n in config.motiflist:
		if motifcount[m] == 0:
			norm_combo1 = 0
		else:
			norm_combo1 = combodict[frozenset([m,n])]/float(motifcount[m])
		if motifcount[m] + motifcount[n] == 0:
			norm_combo2 = 0
		else:
			norm_combo2 = 2*combodict[frozenset([m,n])]/float(motifcount[m] + motifcount[n])
		stats.write("%s," %norm_combo1)
		doublelist1.append(norm_combo1)
		doublelist2.append(norm_combo2)
	doublematrix1.append(doublelist1)
	doublematrix2.append(doublelist2)
	stats.write("%s,%s\n" %(motifdoublecount[m],motifcount[m]))
stats.close()

makeheatmap(doublematrix1,"singlenorm")
makeheatmap(doublematrix2,"doublenorm")
