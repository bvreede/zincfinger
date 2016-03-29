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

import re, sys,config,os.path
import pylab as pl
import numpy as np
import seaborn as sns
pl.rcParams['xtick.labelsize']=8
pl.rcParams['ytick.labelsize']=8



errormsg= "USAGE: findmotif.py path/to/inputfile [screen/run]\nInput needs to be a protein fasta file.\n[screen] is an optional argument; if it is used, most outputfiles are not generated and the only purpose of the run is to count how many motifs are found.\nOutputfolders are indicated in the script; edit the script if you want to alter them."

if len(sys.argv) <= 1:
	sys.exit(errormsg)

try:
	rtype = sys.argv[2]
	if rtype != 'screen':	#option that only screens the total number of hits, and doesn't do any further analyses
		sys.exit(errormsg)
except IndexError:
	rtype = 'run'

#what kind of motiflist/motiflength/motifdict/translationdict to use:
if rtype == 'screen':
	motiflist = config.motiflist1
	motiflength = config.motiflength1
	motifdict = config.motifdict1
	translationdict = config.translationdict1
else:
	motiflist = config.motiflist2
	motiflength = config.motiflength2
	motifdict = config.motifdict2
	translationdict = config.translationdict2


#set figure aesthetics
sns.set_style("white")
sns.set_context("poster")


infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]
hmmfile = "%s/%s/%s_hmmsearch.txt" %(config.mainfolder,config.hmmfolder,infilebrev)

#NB! filenames should always start with an input ID specifier (separate elements in dashes) and end with output ID specifiers.
#e.g.: 150525-dmel_seq.fa or 141212-tcas_heatmap.svg


### OUTPUT FILES ###
# images: bar graph and heatmap
bargraphfig = "%s/%s/%s_motifcount" %(config.mainfolder,config.imgfolder,infilebrev)
heatmapfig = "%s/%s/%s_motifoverlap" %(config.mainfolder,config.imgfolder,infilebrev)

# databases: hit results by site
hitsdb = "%s/%s/%s_hitsdb" %(config.mainfolder,config.dbfolder,infilebrev)

motseq = "%s/%s/%s-motseq" %(config.mainfolder,config.dbfolder,infilebrev)
allmotifs = "%s/%s/%s_allmotifs" %(config.mainfolder,config.dbfolder,infilebrev) #for frequency and aa sequence of specific motifs
results = "%s/%s/%s_motifcount" %(config.mainfolder,config.resfolder,infilebrev)
# fasta files: translated hits, aa sequence of hits
transdb = "%s/%s/%s_protstring" %(config.mainfolder,config.seqmfolder,infilebrev)
statsdb = "%s/%s/%s_motifstats" %(config.mainfolder,config.resfolder,infilebrev)


# Read the fasta file, open output files

fastadb = open("%s" %(infile))
hmmdb = open("%s" %hmmfile)

if rtype != 'screen':
	outputdb = open("%s.csv" %hitsdb, "w")
	outfasta = open("%s.fa" %transdb, "w")
	allmotifsfa = open("%s.fa" %allmotifs, "w")
	allmotifstxt = open("%s.txt" %allmotifs, "w")
	for m in motiflist:
		mfile = open("%s-%s.fa" %(motseq,m), "w")
		mfile.close()
	outfasta = open("%s.fa" %transdb, "w")


#Make dictionaries to count all motifs and combinations of motifs.
motifcount = {m: 0 for m in motiflist}
motifdoublecount = {m: 0 for m in motiflist}
nonambcount = {m: 0 for m in motiflist}
combodict = {a: 0 for a in set([frozenset([m,n]) for m in motiflist for n in motiflist])}


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
			if rtype != 'screen':
				allmotifstxt.write('{' + regex[:-1] + '}\n') #add the regex to the 'allmotifs' document for frequency-dependent sampling
			for m in mots[n]:
				motifdoublecount[m] += 1
				if mots[n].count(m) > 1:
					combodict[frozenset([m])] += (1./mots[n].count(m)) #divide by total count, otherwise it will count +2 or more if it passes this point twice
		else:
			nonambcount[mots[n][0]] += 1 # counter only for nonambiguous motifs
			transl += translationdict[mots[n][0]]
			if rtype != 'screen':
				allmotifstxt.write("%s\n" %translationdict[mots[n][0]]) #add the (translated) motif to the 'allmotifs' document for frequency-dependent sampling
	return transl


def makestackedbargraph(values1,values2,labels,name):
	'''
	Makes a stacked bar chart with two sets of values on y and labels on x.
	'''
	fig = pl.figure()
	ind = np.arange(len(values1))
	p1 = pl.bar(ind,values1,color="navy")
	p2 = pl.bar(ind,values2,color="darkkhaki",bottom=values1)
	pl.xticks(ind + 0.5, labels, rotation=90)
	#pl.yticks(np.arange(30000,60000,5000)) #only used to restrict plot values
	#pl.ylim((30000,60000)) #only used to restrict plot values
	pl.yticks(np.arange(0,7000,1000)) #only used to restrict plot values
	pl.ylim((0,7000)) #only used to restrict plot values
	pl.savefig("%s-%s.svg" %(bargraphfig,name))
	pl.clf()
	pl.close()

def makeheatmap(matrix,labelsx,labelsy,name):
	'''
	Makes a heatmap of a matrix, using Seaborn.
	'''
	sns.heatmap(matrix, xticklabels=labelsx,yticklabels=labelsy,linewidths=.5)#,cmap="YlOrBr")
	pl.savefig("%s-%s.svg" %(heatmapfig,name))
	pl.clf()
	pl.close()

def hmmdicter(infile):
	'''
	Takes the result file from hmmsearch and turns it into a dictionary
	which has gene names as key, and hits (start site, sequence, end site)
	as values. All hits are collected as lists in a list, so that
	multiple hits are saved with the same protein as key.
	'''
	hmmdict= {}
	linename = ""
	for line in infile:
		try:
			lineli = line.split()
			linehead = lineli[0]
		except IndexError:
			continue
		if linehead == ">>":
			linename = lineli[1]
		if lineli[0] == linename: # this line contains positional information and sequence of the hmm hit
			if linehead in hmmdict: #protein already has an entry; append information
				# in hmmdict: list of lists of hits
				newvalue = hmmdict[linehead]
				newvalue += [lineli[1:]]
				hmmdict[linehead] = newvalue
			else: #make new entry for the protein
				hmmdict[linehead] = [lineli[1:]] #lineli in 2nd dimension because further hits may follow
	return hmmdict

def test_hmmentry(strt,key):
	'''
	Tests a hit found by the regular expression against the db created
	with the pfam hmm. If it is found, return 1. If not, return 0 (this
	will then reject the regular expression hit).
	'''
	verify = 0
	hits = hmmdict[key]
	for h in hits:
		strt_h, end_h = int(h[0])-8,int(h[2])+8
		if strt >= strt_h and strt <= end_h:
			verify = 1
	return verify


# Make the headers for the outputdb
if rtype != 'screen':
	outputdb.write("Gene_stable_ID,Gene_name,Protein_stable_ID,Sequence_length,")
	for m in motiflist:
		outputdb.write("%s," %m)
	outputdb.write("\n")


fastadict = config.fastadicter(fastadb) # translate the fasta file into a dictionary

hmmdict = hmmdicter(hmmdb) # translate the hmmer output into a dictionary


seqdict = {} # dictionary for [start position]: sequence
motdict = {} # dictionary for [start position]: motif type
# go through the sequences

for key in fastadict:
	# get info for the first columns (ID and sequence length)
	ids = key.split('|')
	seqlen = len(fastadict[key]) #length of the sequence
	# screen the sequence for motifs
	seqdict.clear() #for the fasta file: collect positions as key and the corresponding sequence at that position as value
	motdict.clear() #as seqdict, but with the motif name instead of sequence
	poslist = [] #for the fasta file: collect all positions to put them in order later on
		# check each motif individually
	if key not in hmmdict:
		continue
	if rtype != 'screen':
		outputdb.write("%s,%s,%s,%s," %(ids[0],ids[1],ids[2],seqlen)) #turn the header name into gene ID/name/prot ID
	for m in motifdict: #go through each motif and find all instances in the sequence. NB: m is a regular expression.
		thisseqcount = 0 #per motif per seq, to give an index for each aminoacid sequence found
		if rtype != 'screen':
			mfile = open("%s-%s.fa" %(motseq,m), "a")
		domain = motifdict[m]
		for i in domain.finditer(fastadict[key]):
			mseq = i.group() # the sequence picked up by the RE
			strt = i.start() + config.plink
			# test whether this hit was found also by the pfam screen: hmmdict
			hmmverify = test_hmmentry(strt,key)
			if hmmverify == 0:
				continue
			motifcount[m] += 1 # count the found motif
			if rtype != 'screen':
				mfile.write(">%s\n%s\n\n" %(key,mseq))
				allmotifsfa.write(">%s|%s-%s\n%s\n\n" %(key,translationdict[m],thisseqcount,mseq))
			thisseqcount += 1
			if strt in seqdict:
				ns = seqdict[strt] + "/" + mseq
				nm = motdict[strt] + "/" + m
				seqdict[strt] = ns
				motdict[strt] = nm
			else:
				seqdict[strt] = mseq
				motdict[strt] = m
				poslist.append(strt)
		if rtype != 'screen':
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
		if rtype != 'screen':
			outputdb.write("%s," %outstring[:-1])
	if rtype != 'screen':
		outputdb.write("\n")

	# another for clustering
	transl = translation(posmatrix,motdict,seqdict)
	if len(transl) > 0 and rtype != 'screen':
		outfasta.write(">%s\n%s\n\n" %(key,transl))

if rtype != 'screen':
	outfasta.close()
	allmotifstxt.close()
	outputdb.close()


'''
make the stats!
in rows: for each motif, count the times it coincides with another motif, and divide this by the total of hits
for that motif.
Make a heatmap with the stats.
'''
if rtype != 'screen':
	#open statsdb collecting overlap
	overlapstats = open("%s.csv" %statsdb, "w")
	mcounts,mcounts_na = [],[]
else:
	# open results file for this species/inputfile
	sumout = open("%s.txt" %results, "w")

for m in motiflist:
	if rtype == 'screen':
		sumout.write("%s_total: %s\n" %(m,motifcount[m])) 	# for total counts in results text
		sumout.write("%s_na: %s\n" %(m,nonambcount[m]))		# for non-ambiguous counts in results text
	else:
		overlapstats.write(",%s" %m)				# headers for stats on overlapping motifs
		mcounts.append(motifcount[m])				# for total counts in bar graph
		mcounts_na.append(nonambcount[m])			# for non-ambiguous counts in bar graph
if rtype != 'screen':
	overlapstats.write("\n")
	#per motif count combinations
	doublematrix = []
	for m in motiflist:
		overlapstats.write("%s" %m)
		doublelist = [] #rows for heatmap matrix
		for n in motiflist:
			overlapstats.write(",%s" %combodict[frozenset([m,n])])
			if motifcount[m] == 0:
				norm_combo = 0.0 #prevent dividing by 0
			else:
				norm_combo = combodict[frozenset([m,n])]/float(motifcount[m])
			doublelist.append(norm_combo)
		overlapstats.write("\n")
		doublematrix.append(doublelist)
	#close output files
	overlapstats.close()
else:
	sumout.close()





# stacked bar graph of nonambiguous + ambiguous motifs:
#makebargraph(mcounts,motiflist,"motifs")
#makebargraph(mcounts_na,motiflist,"non-ambiguous_motifs")
if rtype != 'screen':
	mcounts_am = [p-mcounts_na[n] for n,p in enumerate(mcounts)]
	makeheatmap(doublematrix,motiflist,motiflist,"singlenorm")
	makestackedbargraph(mcounts_na,mcounts_am,motiflist,"motifs-stacked")

