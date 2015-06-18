import config,csv,itertools,re,random
from jellyfish import levenshtein_distance as jld

idr = "150602-SM00355" #the identifier for all input files (motif sequences)


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



# GET INPUT AND GENERATE (1) list of orth combinations and (2) motif sequence dictionary
if __name__ == "__main__":
	# add option to use only limited species here:
	spp = ['dmel','nvec']
	orthfolder = "%s/%s" %(config.mainfolder,config.orthfolder)
	seqfolder = "%s/%s" %(config.mainfolder,config.seqfolder)
	dbfolder = "%s/%s" %(config.mainfolder,config.dbfolder)
	orthcombos = []
	msequencedx = {}
	randommotifs = {}
	for sp in spp:
		# make a list with random motif elements for this species and save it in the dictionary
		motifdb = open("%s/%s-%s_hmmallmotifs.txt" %(dbfolder,idr,sp))
		motli = [line.strip() for line in motifdb]
		randommotifs[sp] = motli
		# open appropriate motif sequence files and import them into a dictionary
		msequences = open("%s/%s-%s_hmmprotstring.fa" %(seqfolder,idr,sp))
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

	# generate list for further detailed comparisons
	detcomp = []
	
	# for each ortholog combination:
	for y,x in enumerate(orthcombos):
		print "Ortholog combination %s of %s..." %(y,len(orthcombos))
		# identify the ortholog combination
		xli = list(x)
		
		m,n = xli
		if m not in msequencedx or n not in msequencedx:
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
			# calculate lowest levenshtein distance
			d=simple_wordcomp(s1,s2)
			if d == 0: # if distance is 0, the motif sequences are identical
				if n == 0: #ortholog combo
					orthid += 1 # add to 'identical'
					detcomp.append(x) # if this was ortholog combo, save their names for further processing.
					continue
				else:
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
					ransub += 1
					continue
			# remove Z and calculate levenshtein distance again
			s1_ = s1.replace('Z','')
			s2_ = s2.replace('Z','')
			d_ = simple_wordcomp(s1_,s2_)
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

# For ortholog and for random: make a pie chart with the results (or just output them as numbers and manually make a pie chart, whatever)
print "ORTHOLOGS:\n-identical %s\n-substitution %s\n-structure %s\n-addition %s\n-other %s\n" %(orthid,orthsub,orthstruc,orthadd,orthother)
print "RANDOM:\n-identical %s\n-substitution %s\n-structure %s\n-addition %s\n-other %s\n" %(ranid,ransub,ranstruc,ranadd,ranother)
print "TOTAL: %s/%s" %((orthid+orthsub+orthstruc+orthadd+orthother),(ranid+ransub+ranstruc+ranadd+ranother))

