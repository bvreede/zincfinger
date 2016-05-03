#!/usr/bin/python

'''
This script takes the motif sequences (written for zinc fingers, but
can be any kind of motif) of ortholog pairs, and compares them. The output
is the numbers of pairs in each of the following categories:
(1) ortholog pairs with identical motif sequences;
(2) those with motif substitutions (but identical structures);
(3) those with only differences in structure;
(4) those with only additions/deletions explaining the difference between the pair;
(5) all others (i.e. combinations of additions/deletions, substutions and structural differences).
For each pair, a random sequence is generated and subsequently compared
with one of the proteins of the pair, then categorized in the same way.

Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 10 August 2015
'''

import config,csv,itertools,re,random,os
from jellyfish import levenshtein_distance as jld
import pylab as pl
from matplotlib import cm
from numpy import arange


#### WHAT SPECIES TO USE: CHANGE IT HERE!! ###
#spp = config.chor
spp = config.arth
#spp = config.sppall
#spp = config.spp700

### DON'T FORGET TO CHANGE THE OUTPUT NAME ACCORDINGLY!! ###
#outname = "chor"
#outname = "d700"
outname = "arth"
#outname = "sppall"

###########################################


orthfolder = "%s/%s" %(config.mainfolder,config.orthfolder)
seqmfolder = "%s/%s" %(config.mainfolder,config.seqmfolder)
seqpfolder = "%s/%s" %(config.mainfolder,config.seqpfolder)
dbfolder = "%s/%s" %(config.mainfolder,config.dbfolder)
resfolder = "%s/%s" %(config.mainfolder,config.resfolder)

idr = config.idr

#if necessary: run findmotif.py with the orthology option on all input species fasta files
for sp in spp:
	if os.path.exists("%s/%s-%s_allmotifs.txt" %(dbfolder,idr,sp)):
		continue
	else:
		for f in os.listdir("%s/%s" %(config.mainfolder,config.seqpfolder)):
			print "running 'findmotif' on", f
			os.system("findmotif.py %s/%s orthinfo" %(seqpfolder,f))
		break


def updatedx(tempdx,sp):
	'''
	Creates a new dictionary from an existing one
	where the keys are slightly modified to accommodate
	the rest of the script.
	(Original keys are geneID|genename|protID; new keys
	are species|protID.)
	'''
	newdx = {}
	for key in tempdx:
		prot = key.split('|')[2]
		newkey = "%s|%s" %(sp,prot)
		newdx[newkey] = tempdx[key]
	return newdx

def simple_wordcomp(i,j):
	'''
	Calculate pairwise distances for all the strings collected. As strings
	are regular expressions, first expand the re. and then calculate all distances
	pairwise. Return the minimal distance (with and without spaces).
	'''
	# translate strings[i] and strings[j] to all possible expressions
	sisplit = [part.split('|') for part in re.split(r'\{(.*?)\}',i)]
	sjsplit = [part.split('|') for part in re.split(r'\{(.*?)\}',j)]
	si,sj = [],[]
	for x in itertools.product(*sisplit):
		si.append(''.join(x))
	for x in itertools.product(*sjsplit):
		sj.append(''.join(x))
	# check jld of all strings[i] options against all strings[j] options
	distance = []
	for sin in si:
		for sjn in sj:
			distance.append(jld(sin,sjn))
	# return the minimum distance
	return min(distance)

def longregex_wordcomp(i,j):
	'''
	Less precise but faster option for distance calculation between strings
	with a lot of regular expressions (expansion of many regular expression
	can generate thousands if not millions of individual comparisons, slowing
	down the script significantly).
	'''
	li = config.re2li(i)
	lj = config.re2li(j)
	lboth = li+lj
	used = [e for e in lboth if len(e) == 1] #which individual elements are used? These should not be duplicated in the alphabet translation of re elements.
	# make a dictionary to translate regular expression elements
	tempdict = {}
	for e in used:
		tempdict[e] = e
	k = 0
	for e in lboth:
		if len(e) > 1:
			if e not in tempdict:
				tdk = '' #temporary key; otherwise the dictionary changes during use
				for key in tempdict:
					d = simple_wordcomp(e,key)
					if d == 0: #for regular expressions that are closely related, the same translation can be used
						tdk = tempdict[key]
						pass
				if tdk != '':
					tempdict[e] = tdk
					continue
				# if it gets here, no similar key has been found. So make one!
				while config.alphabet[k] in used:
					k += 1
				tempdict[e] = config.alphabet[k]
	# now translate each into a new string
	ni,nj = "",""
	for e in li:
		if len(e) == 1:
			ni += e
		else:
			ni += tempdict[e]
	for e in lj:
		if len(e) == 1:
			nj += e
		else:
			nj += tempdict[e]
	# now run the comparison between strings
	d = jld(ni,nj)
	return d

def randomorth(gene):
	'''
	Generate a random motif sequence based on identical ZF structure
	(including conservation of spacing) but with random frequency-dependent
	sampling of motifs.
	'''
	regex = ['{', '|'] #components of regular expressions, to be excluded
	# NB } is not used because it indicates the re itself, and thus needs to be substituted by an element, too
	rehit = 0
	seq_g = msequencedx[gene]
	sp = gene.split('|')[0]
	motli = randommotifs[sp]
	ranseq = ""
	for l in seq_g:
		if l == 'Z':
			ranseq += 'Z'
		elif l in regex: #regex component followed by a motif
			rehit = 1
		elif rehit == 1: #catches the motif following the regex indicator
			rehit = 0
		else:
			ranseq += random.choice(motli)
	return ranseq

def lengthre(s):
	'''
	determine how many *actual* elements are in a string
	that contains regular expressions
	'''
	n = len(s) - s.count('{') - s.count('}') - s.count('|') * 2
	return n

def zindex(li,z):
	'''
	returns a list of the indices of string 'z' in list 'li'
	'''
	zi = [n for n,i in enumerate(li) if i == z]
	return zi


def piechart(data,labels,name):
	'''
	makes a piechart with the data
	'''
	#define the color scheme
	cs=cm.RdBu(arange(4)/3.)
	cs = [cs[0],cs[1],cs[3],cs[2]]
	#make the pie
	pl.figure(1,figsize=(10,9.8))
	pl.pie(data,startangle=90, colors=cs)
	#make the legend
	pl.legend(labels,loc=(-0.05,0.9))
	pl.savefig("%s/%s/%s-conservationpie-%s.png" %(config.mainfolder,config.imgfolder,idr,name))


# GET INPUT AND GENERATE (1) list of orth combinations and (2) motif sequence dictionary
if __name__ == "__main__":
	orthcombos = []
	msequencedx = {}
	randommotifs = {}
	for sp in spp:
		# make a list with random motif elements for this species and save it in the dictionary
		motifdb = open("%s/%s-%s_allmotifs.txt" %(dbfolder,idr,sp))
		motli = [line.strip() for line in motifdb]
		randommotifs[sp] = motli
		# open appropriate motif sequence files and import them into a dictionary
		msequences = open("%s/%s-%s_protstring.fa" %(seqmfolder,idr,sp))
		tempdx = config.fastadicter(msequences) #dictionary 1: original dictionary from fasta file
		newdx = updatedx(tempdx,sp) #dictionary 2: updated keys
		msequencedx.update(newdx) #import dictionary 2 to the main sequencedx.
		# open all ortholog files and import combinations as frozensets
		orthologs = csv.reader(open("%s/%s-allorth.csv" %(orthfolder,sp)))
		for o in orthologs:
			if o[2] in spp: #only save if species combination is in the list you want to check
				o1 = "%s|%s" %(o[0],o[1])
				o2 = "%s|%s" %(o[2],o[3])
				orthcombos.append(frozenset([o1,o2]))
	# make it into a set, so that duplicates are removed
	orthcombos = list(set(orthcombos))

	# generate 10 categories:
	# Ortholog-identical
	# Ortholog-substitution
	# Ortholog-structure
	# Ortholog-addition/subtraction
	# Ortholog-other
	orthid,orthsub,orthstruc,orthadd,orthother = 0,0,0,0,0
	# Random-identical
	# Random-substitution
	# Random-structure
	# Random-addition/subtraction
	# Random-other
	ranid,ransub,ranstruc,ranadd,ranother = 0,0,0,0,0
	notcounted = 0

	# generate list for further detailed comparisons
	detcomp,detcompran = [],[]
	
	# for each ortholog combination:
	for y,x in enumerate(orthcombos):
		if y%100 == 0: #only when number is divisible by 100
			print "Ortholog combination %s of %s..." %(y+1,len(orthcombos))
		# identify the ortholog combination
		xli = list(x)
		m,n = xli
		if m not in msequencedx or n not in msequencedx:
			notcounted+=1
			continue
		# retrieve motif sequence for orthologs
		mseq = msequencedx[m]
		nseq = msequencedx[n]
		# shuffle to identify the random combination
		random.shuffle(xli)
		o,p = xli
		# retrieve motif sequence for random
		oseq = msequencedx[o]
		pseq = randomorth(p)
		# for both sequence combinations:
		for n,combo in enumerate([[mseq,nseq],[oseq,pseq]]): #enumerate to identify ortholog (0) v random (1) comparison
			s1,s2 = combo
			# calculate lowest levenshtein distance: first determine what program to use!
			if s1.count('|') + s2.count('|') > 16:
				d=longregex_wordcomp(s1,s2)
			else:
				d=simple_wordcomp(s1,s2)
			if d == 0: # if distance is 0, the motif sequences are identical
				if n == 0: #ortholog combo
					orthid += 1 # add to 'identical'
					detcomp.append(x) # if this was ortholog combo, save their names for further processing.
					continue
				else:
					detcompran.append(combo) # if this was a random combo, save their sequences
					ranid += 1
					continue
			# if lengths of sequences and Z indices are the same, substitution explains the difference
			l1 = config.re2li(s1)
			l2 = config.re2li(s2)
			z1 = zindex(l1,'Z')
			z2 = zindex(l2,'Z')
			if len(l1) == len(l2) and z1 == z2:
				if n == 0: #ortholog combo
					orthsub += 1 #add to 'substitution'
					detcomp.append(x) # if this was ortholog combo, save their names for further processing.
					continue
				else:
					detcompran.append(combo) # if this was a random combo, save their sequences
					ransub += 1
					continue
			# remove Z and calculate levenshtein distance again
			s1_ = s1.replace('Z','')
			s2_ = s2.replace('Z','')
			if s1_.count('|') + s2_.count('|') > 16:
				d_=longregex_wordcomp(s1_,s2_)
			else:
				d_=simple_wordcomp(s1_,s2_)
			if d_ == 0: # if distance is 0, structure explains the difference
				if n == 0: #ortholog combo
					orthstruc += 1
					continue
				else:
					ranstruc += 1
					continue
			# if distance is equal to difference in length, addition of motifs explains the difference
			if d_ == abs(lengthre(s1_) - lengthre(s2_)):
				if n == 0:
					orthadd += 1
				else:
					ranadd += 1
			else:
				if n == 0:
					orthother += 1
				else:
					ranother += 1

# OUTPUT OF RESULTS #
if __name__ == "__main__":
	outcsv = open("%s/%s-%s_orthcomp.csv" %(resfolder,idr,outname), "w")
	outcsv.write("Orthologs:\nidentical,%s\nsubstitution,%s\nstructure,%s\naddition,%s\nother,%s\n\n" %(orthid,orthsub,orthstruc,orthadd,orthother))
	outcsv.write("Random:\nidentical,%s\nsubstitution,%s\nstructure,%s\naddition,%s\nother,%s\n\n" %(ranid,ransub,ranstruc,ranadd,ranother))
	outcsv.write("Total:\northologs,%s\nrandom,%s\nnot_assessed,%s" %((orthid+orthsub+orthstruc+orthadd+orthother),(ranid+ransub+ranstruc+ranadd+ranother),notcounted))
	outcsv.close()

	print "ORTHOLOGS:\n-identical %s\n-substitution %s\n-structure %s\n-addition %s\n-other %s\n" %(orthid,orthsub,orthstruc,orthadd,orthother)
	print "RANDOM:\n-identical %s\n-substitution %s\n-structure %s\n-addition %s\n-other %s\n" %(ranid,ransub,ranstruc,ranadd,ranother)
	print "TOTAL: %s/%s (%s not counted)" %((orthid+orthsub+orthstruc+orthadd+orthother),(ranid+ransub+ranstruc+ranadd+ranother),notcounted)
	labels = ["identical","substitution",'addition','addition+substitution']
	data4pie = [orthid,orthsub,orthadd,orthother+orthstruc]
	data4pie_ran = [ranid,ransub,ranadd,ranother+ranstruc]
	piechart(data4pie,labels,outname+"-orth")
	piechart(data4pie_ran,labels,outname+"-ran")

# PROCESSING DETAILED COMPARISONS TO A NEW FILE TO INPUT ELSEWHERE #
if __name__ == "__main__":
	detcompcsv = open("%s/%s-%s_orthcomp-detail.csv" %(resfolder,idr,outname), "w")
	detcomprancsv = open("%s/%s-%s_orthcomp-detail-random.csv" %(resfolder,idr,outname), "w")
	for d in detcomp:
		xli = list(d)
		m,n = xli
		mseq = msequencedx[m]
		nseq = msequencedx[n]
		detcompcsv.write("%s,%s,%s,%s\n" %(m,mseq,n,nseq))
	for d in detcompran:
		mseq,nseq = d
		txt = "random_sequence"
		detcomprancsv.write("%s,%s,%s,%s\n" %(txt,mseq,txt,nseq))
	detcompcsv.close()
	detcomprancsv.close()
