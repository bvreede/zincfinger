#!/usr/bin/python

'''
This script takes the selected motif sequences that
resulted from 'orthconservation.py' and outputs which
motifs replace which. Further, it looks in detail at
the amino acid sequences of motifs that have not been
substituted between orthologs, to score substitutions
in the amino acid sequence itself.

Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 10 August 2015
'''


import sys,config,csv,itertools,random
import pylab as pl
import numpy as np

if len(sys.argv) <= 1:
	sys.exit("USAGE: python orthconservation-part2.py path/to/inputfile (input file is the '-detail.csv' output of orthconservation.py).\nOutputfolders are indicated in the script; edit the script if you want to alter them.")


## make stacked bargraphs, splitting the alternative and standard motifs
## give random-model a different colour if possible



### INPUT FILES ###
infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]
orthin = [line for line in csv.reader(open(infile))] #this opens the file with ortholog combinations to investigate
motifaaseq = "%s/%s/%s-alli_allmotifs.fa" %(config.mainfolder,config.dbfolder,config.idr)

### OUTPUT FILES ###
#heatmaptxt = "%s/%s/%s_conservation-heatmap" %(config.mainfolder,config.evfolder,infilebrev)
#bargraphtxt = "%s/%s/%s_conservation-bargraph" %(config.mainfolder,config.evfolder,infilebrev)
seqconshm =  "%s/%s/%s_sequencehm" %(config.mainfolder,config.imgfolder,infilebrev)


### LIST OF SPECIES FOR COMPARISON ###
# amino acid sequence comparisons are only made between one member of group1 and one of group2
group1 = config.anim
group2 = config.plan + config.prot


#make dictionary to store conservation/non-conservation per motif
consdict = {m:[0,0] for m in config.motiflist2} #first digit for conserved, second for not-conserved

# make dictionaries where motif combos can be counted.
ambiguous = {} #dictionary to count appearance of ambiguous combinations
for m in config.motiflist2:
	m_ = config.translationdict2[m]
	for n in config.motiflist2:
		n_ = config.translationdict2[n]
		mn = frozenset([m_,n_])
		ambiguous[mn] = 0

# copy combination dictionaries
orthambi = dict(ambiguous) # same as ambiguous but this one to count ambiguous combinations that are orthologous
substitutions = dict(ambiguous) # same as ambiguous, to count substitutions between orthologs

# make dictionary where motifs individually can be counted
individual = {}
for m in config.motiflist2:
	n = config.translationdict2[m]
	individual[n] = 0

# make dictionary of amino acid sequences for all species
# and rewrite the dictionary to have different headers (proteinID|motifclass-number)
aadict = {}
motifsaa = open(motifaaseq)
tempdict = config.fastadicter(motifsaa)
for key in tempdict:
	newhead = key.split('|')[2] + '|' + key.split('|')[3]
	aadict[newhead] = tempdict[key]


# make dictionary of motifs as key, with all amino acid sequences in a list as value.
# this can be used to pick a random motif later.
motrandomdx = {config.translationdict2[m]: [] for m in config.motiflist2}
for key in aadict:
	motif = key.split('|')[1][0] #the letter indicating motif type
	sequence = aadict[key]
	motrandomdx[motif].append(sequence)

# make dictionary for conservation measurements of amino acid sequences
conservation,conserv_rand = {},{} #dictionary where lists of conservation per site are stored per motif, and same for random comparisons
for m in config.motiflist2:
	mlen = config.motiflength2[m] + config.plink + config.alink + 1
	mli = [0 for n in range(mlen)] #a 0 for each site, and finally a 0 that will count how many times this motif was found conserved
	conservation[m] = list(mli)
	conserv_rand[m] = list(mli)


def re_move(s):
	'''
	returns a list with the frozenset-combinations of
	all different elements of a single regular expression.
	'''
	s = s.replace('{','')
	s = s.replace('}','')
	sli = s.split('|')
	sli_pairs = [frozenset(list(pair)) for pair in itertools.combinations(sli,2)]
	return sli_pairs


def aacomp_det(aa1,aa2,li):
	'''
	Function for the detailed comparison of two aminoacid sequences
	(called within aacomp).
	'''
	#for each element in aa1:
	for n,el1 in enumerate(aa1):
		el2 = aa2[n]
		if el1 != el2: #compare with the corresponding letter in aa2
			li[n] += 1 #adjust the list with +1 if it is different, and not if it is not different
	li[-1] += 1 #update the final element in the list with +1
	return li

def aacomp(gene1,m1,gene2,m2):
	'''
	Function that compares the aminoacid sequences of homologous
	motifs to identify conserved aminoacids.
	'''
	#compile gene names
	g1 = gene1 + '|' + m1 
	g2 = gene2 + '|' + m2
	m = m1[0] #the motif identifier
	#get the random sequence
	mrandomlist = motrandomdx[m] #returns a list with all sequences of this motif in the database
	rand_n = random.randint(0,len(mrandomlist)-1) #pick a random number to serve as index
	s_r = mrandomlist[rand_n]

	#collect sequences to be compared
	s1 = aadict[g1]
	s2 = aadict[g2]

	#motif name, and get the global dictionary for motif lists
	motname = config.translationdict_inv2[m]
	global conservation
	global conserv_rand
	
	#now run comparisons
	### PART 1: compare gene1 + gene 2
	li_mot = conservation[motname]
	limot_updated = aacomp_det(s1,s2,li_mot)
	conservation[motname] = limot_updated #update global dictionary with limot_updated

	### PART 2: compare either gene 1 or gene 2, and random
	#pick gene 1 or 2
	if random.randint(1,2) == 1:
		s_g = s1
	else:
		s_g = s2
	li_mot_ran = conserv_rand[motname]
	limot_updated_ran = aacomp_det(s_g,s_r,li_mot_ran)
	conserv_rand[motname] = limot_updated_ran #update global dictionary with limot_updated

def mordercheck(m,odict):
	'''
	Check the index of a motif in the protein (start with 0 and every subsequent
	hit is +1) with the existing odict, and update odict.
	'''
	m = m.replace('{','').replace('}','').replace('|','')
	for k in m:
		if k in odict:
			odict[k] += 1
		else:
			odict[k] = 0
	return odict

def makeheatmap(table,name,xlab,ylab):
	'''
	Use 2D data (list of lists) to generate a heatmap of the data.
	'''
	data = pl.array(table)
	colourformap = "YlOrBr"
	fig,ax = pl.subplots()
	heatmap = pl.pcolor(data, cmap=colourformap,vmin=0,vmax=1)
	cbar = pl.colorbar(heatmap)
	
	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
	pl.axis('tight') #remove the white bar
	ax.invert_yaxis() #make sure it starts counting from the top
		
	#make the labels
	ax.set_xticklabels(xlab, minor=False, rotation=90)
	ax.set_yticklabels(ylab, minor=False)
		
	# save the figure
	pl.savefig("%s-%s.png" %(seqconshm,name), dpi = 300)
	pl.savefig("%s-%s.svg" %(seqconshm,name), dpi = 300)
	pl.clf()
	pl.close("all")
	



for combo in orthin:
	o1 = config.re2li(combo[1])
	o2 = config.re2li(combo[3])
	o1dict,o2dict = {},{}
	for e,f in zip(o1,o2):
		if e == 'Z': #this is not a motif. Moving on...
			continue

		#determine the order of motifs: update o1dict/o2dict and return the index of the current motifs
		o1dict = mordercheck(e,o1dict)
		o2dict = mordercheck(f,o2dict)

		# Determine orthology and frequency of ambiguous items
		if e.count('{') > 0 or f.count('{') > 0: # ambiguous motif identified.
			# turn re into lists
			eli = re_move(e)
			fli = re_move(f)
			# count each instance (in total)
			for item in eli + fli:
				try:
					ambiguous[item] += 1
				except KeyError:
					continue
			# determine if it is orthologous
			for item in eli:
				if item in fli:
					orthambi[item] += 2
			continue

		# with no ambiguous items: either score a substitution, or compare sequences...
		# score a substitution here:
		if e != f:
			consdict[config.translationdict_inv2[e]][1] += 1 # score a non-conserved in the conservation dictionary
			consdict[config.translationdict_inv2[f]][1] += 1 # for both motifs.
			individual[e] += 1
			individual[f] += 1
			substitutions[frozenset([e,f])] += 1
		# motifs are identical, so compare sequences:
		else:
			consdict[config.translationdict_inv2[e]][0] += 2 # score a 'conserved' in the conservation dictionary for both motifs (which are the same).
			individual[e] += 1
			individual[f] += 1
			# CLAUSE TO ADD: only run this comparison for distant species#
			# fetch species and protein ID
			sp1,prot1 = combo[0].split('|')
			sp2,prot2 = combo[2].split('|')
			if sp1 in group1:
				if sp2 not in group2:
					continue
					# this is passed only if sp1 in group1 and sp2 in group2
			elif sp1 in group2:
				if sp2 not in group1:
					continue
					# this is passed only if sp1 in group2 and sp2 in group1
			else:
				continue
			# what motif number is this?
			m1 = e + '-' + str(o1dict[e])
			m2 = f + '-' + str(o2dict[f])
			aacomp(prot1,m1,prot2,m2)

for s in substitutions:
	if substitutions[s] != 0:
		a,b = list(s)

'''
#MAKE THE HEATMAP FOR AMBIGUOUS APPEARANCE + SUBSTITUTIONS#
doublematrix3,motifcounts = [],[] #for substitution, motif count bar graphs; respectively
for m in config.motiflist2:
	line1,line2,line3 = [],[],[]
	sumsubs = 0
	a = config.translationdict2[m]
	for n in config.motiflist2:
		b = config.translationdict2[n]
		combo = frozenset([a,b])
		# substitution of motifs in orthologs
		if individual[a] > 0:
			line3.append(substitutions[combo]/float(individual[a]))
		else:
			line3.append(0)
		sumsubs += substitutions[combo]
	if individual[a] > 0:
		totalsub = sumsubs/float(individual[a])
	else:
		totalsub = 0
	line3.append(totalsub)
	doublematrix3.append(line3)
	# list for individual motif counts in bar graph
	motifcounts.append(individual[a])
#makeevolviewheatmap(doublematrix3,"substitutions")

#makeevolviewbars(config.motiflist,motifcounts,"counts")
'''

totalambiguous,totalconserved = 0,0
for key in ambiguous:
	totalambiguous += ambiguous[key]
	totalconserved += orthambi[key]

print "motifs as part of an overlapping set: %s, of which conserved in the orthologous site: %s (%s" %(totalambiguous,totalconserved,int(float(totalconserved)/totalambiguous*100)) + "%)"

### MAKE HEATMAP FOR MOTIF CONSERVATION ON SEQUENCE LEVEL ###
##NB this only happens with comparisons distant enough to be specified by 'group1' and 'group2'

countdx = {} #add a dictionary that keeps the counter separate
for m in config.motiflist2:
	consli = list(conservation[m])
	consli_r = list(conserv_rand[m])
	counter = consli[-1]
	#ressum.write("%s,%s\n" %(m,counter))
	countdx[m] = counter
	if counter <10:
		continue
	print "Making heatmap for %s with %s comparisons" %(m,counter)
	consli_new,consli_new_r = [],[]
	for n in range(len(consli)-1):
		c = consli[n]/float(counter)
		consli_new.append(c)
		c_r = consli_r[n]/float(counter)
		consli_new_r.append(c_r)
	conservation[m] = list(consli_new)
	conserv_rand[m] = list(consli_new_r)
	table = [consli_new] + [consli_new_r]
	xlab = []
	ylab = ["orthologous","random"]
	makeheatmap(table,m,xlab,ylab)

