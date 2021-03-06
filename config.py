import re,math,datetime,os.path
import matplotlib.pyplot as plt
import numpy as np
from random import shuffle



# path info
mainfolder = "/home/barbara/Dropbox/zftest"			# Adjust this with initial setup
hmmerbin = "hmmer-3.1b2-linux-intel-x86_64/binaries"		# Location of hmmer binaries, adjust this with initial setup
pfamc2h2 = "hmmer-3.1b2-linux-intel-x86_64/zf-C2H2.hmm"		# Location of pfam hmm model, adjust this with initial setup
ensemblsource = "original_ensembl_dbs"				# Location of zipped ensembl data, adjust this with initial setup
seqfolder = "sequences"
seqpfolder = "sequences/protein"
seqmfolder = "sequences/motifs"
resfolder = "results"
imgfolder = "images"
dbfolder = "data"
orthfolder = "orthologs"
compfolder = "data/compara"
hmmfolder = "data/hmm"

# identifier name
idrpath = "%s/idr.txt" %mainfolder
if os.path.exists(idrpath):
	idrfile = open(idrpath)
	for line in idrfile:
		if len(line.strip()) != 0:
			idr = line.strip()
	idrfile.close()
else:
	idrfile = open(idrpath, "w")
	day = str(datetime.date.today()).replace('-','')
	idr = "%s-ZF" %day #the identifier for all input files (motif sequences)
	idrfile.write(idr)
	idrfile.close()
	

# define motifs
# a complete dataset
moCC = range(1,7) #distances between CC
moCH = range(7,18) #distances between CH
moHH = range(1,7) #distances between HH

# all possible combinations of the above:
motiflist1 = ['%s_%s_%s' %(m,n,o) for m in moCC for n in moCH for o in moHH]


# custom motif list, taken from folder_action result file
motiflist2=[]
semot ="%s/%s/selected_motifs.txt" %(mainfolder,resfolder)
if os.path.isfile(semot):
	mlf = open(semot)
	for line in mlf:
		motiflist2.append(line.strip())

	
standard = ['2_12_3','2_12_4','2_12_5','4_12_3','4_12_4','4_12_5']
alternative = [m for m in motiflist2 if m not in standard]

#distances before and after each C/H
plink,alink = 0,0


#all species							# Adjust this with initial setup
#sppall = ['pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam','dmel','tcas','isca','dpul','smar','turt','atri','atha','crei','cmer','osat','ppat','smoe','slyc','aque','bmal','cele','hrob','lgig','mlei','nvec','sman','spur','tadh','drer','ggal','hsap','mmus','xtro']
sppall=['tthe', 'pfal', 'pinf', 'gthe', 'ehux', 'bnat', 'glam', 'lmaj', 'cmer', 'crei', 'ppat', 'smoe', 'atri', 'osat', 'atha','slyc', 'mlei', 'aque', 'tadh', 'nvec', 'spur', 'drer', 'xtro', 'ggal', 'mmus', 'hsap', 'sman', 'lgig', 'hrob', 'bmal', 'cele', 'isca','turt', 'smar', 'tcas', 'dmel', 'dpul']



#individual groups
chor = ['drer','ggal','hsap','mmus','xtro'] #chordates/vertebrates
arth = ['dmel','tcas','isca','dpul','smar','turt'] #arthropods
eani = ['nvec','mlei','tadh','sman','aque'] #other animals
plan = ['atri','ppat','crei','atha','slyc','osat','smoe','cmer'] #plants
prot = ['pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam'] #protists
anim = ['hrob','spur','lgig','bmal','cele','drer','ggal','hsap','mmus','xtro','dmel','tcas','isca','dpul','smar','turt','nvec','mlei','tadh','sman','aque'] #all animals

#species with >700mya distance
spp700 = ['xtro','spur','tcas','sman','cele','nvec','mlei','tadh','aque','atha','crei','cmer','pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam',]


# Define zf-domains (by the C-H distances), and save in various dictionaries:
	# motiflength: dictionary of motifs and their lengths
	# motifdict: dictionary of motifs and their regular expression
	# translationdict: dictionary to translate ZF-domains to a single letter
def make_motif_dict(motiflist):
	motiflength,motifdict = {},{}
	for m in motiflist:
		cc,ch,hh = m.split('_')
		motif = '[A-Z]{%s}C[A-Z]{%s}C[A-Z]{%s}H[A-Z]{%s}H[A-Z]{%s}' %(plink,cc,ch,hh,alink) # construct the regular expression
		remotif = re.compile(motif)
		motifdict[m] = remotif
		d = [int(x) for x in cc,ch,hh]
		dl = sum(d) + 4
		motiflength[m] = dl
	return motiflength,motifdict

motiflength1,motifdict1 = make_motif_dict(motiflist1)
motiflength2,motifdict2 = make_motif_dict(motiflist2)

#short alphabet
alphabet = '''ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_+=[];\<,>.?/'''

#long alphabet
alphabet = '''~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`'''
translationdict1 = {motif: alphabet[a] for a,motif in enumerate(motiflist1)} # dictionary of motifs and the corresponding string element
translationdict2 = {motif: alphabet[a] for a,motif in enumerate(motiflist2)} # dictionary of motifs and the corresponding string element


translationdict_inv1 = {a: motif for motif,a in translationdict1.items()}
translationdict_inv2 = {a: motif for motif,a in translationdict2.items()}


def fastadicter(fastadb):
	'''
	From an open fasta file, generate a dictionary containing the header as key, and the
	sequence as value.
	'''
	fastadict = {}
	sequence,header = "",""
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


def makebarchart(values,labels,name):
	'''
	Makes a simple bar chart with values on y and labels on x.
	'''
	fig = plt.figure()
	ind = np.arange(len(values))
	plt.bar(ind,values,color="c")
	plt.xticks(ind + 0.5, labels, rotation=90)
	plt.savefig("%s.svg" %name)

def getColour(maxcol):
	'''
	Function to translate numbers into a hex colour.
	Input required: the total number of colours needed. Returns
	a list of colours as long as (or longer) than the number.
	'''
	a = 1/3.
	n = int(math.pow(maxcol,a)) # the number of elements there have to be from 00-FF (minus one, because int is rounded down)
	k = 255/n # the space in integers from 0-255 between each element
	CC,CR,CG,CB,colours = [],[],[],[],[]
	# construct the list of elements from 00-FF
	for i in range(n+1):
		hn = hex(i*k)[2:]
		if len(hn) < 2:
			hn = hn+hn
		CC.append(hn)
	#red: pick each element (n+1)^2 times before moving on to the next
	for c in CC:
		for r in range(pow((n+1),2)):
			CR.append(c)
	#green: pick each element (n+1) times before moving on to the next; repeat (n+1) times
	for g in range(n+1):
		for c in CC:
			for h in range(n+1):
				CG.append(c)
	#blue, pick each element once before moving on to the next, repeat (n+1)^2 times
	for b in range(pow((n+1),2)):
		for c in CC:
			CB.append(c)
	for X,red in enumerate(CR):
		colour = '#' + red + CG[X] + CB[X]
		colours.append(colour)
	shuffle(colours)
	return colours


def sortdata(thresh,nclust,cldict,outfolder):
	'''
	interpret the clustering and apply it to a file with data (e.g. visualization,
	GO terms, etc), so that these are clustered similarly.
	Turned into a function so it can be repeated with different thresholds.
	Defunct function, only kept here for archival purposes.
	'''
	global pID,infile,ftchead #fix this before applying the function anywhere!
	for n in range(1,nclust+1):
		clusterfile = open("%s/%s_cluster%s.csv" %(outfolder,infile.split('/')[-1][:-4],n), "w")
		#write header
		lcollect = ""
		for item in ftchead:
			lcollect += item + ','
		clusterfile.write("%s\n" %lcollect[:-1])
		# write content
		for line in ftc:
			if cldict[line[pID]] == n:
				lcollect = ""
				for item in line:
					lcollect += item + ','
				clusterfile.write("%s\n" %lcollect[:-1])
		clusterfile.close()

def re2li(s):
	'''
	Translate a string containing regular expressions to a list
	where each element of the re occupies a single item.
	'''
	reli = []
	flag = 0
	re_ele = ['{','|','}']
	for i in s:
		if i not in re_ele:
			if flag == 0:
				reli.append(i)
			if flag == 1:
				i_ += i
		else:
			flag = 1
			if i == '{':
				i_ = i
			elif i == '}':
				i_ += i
				reli.append(i_)
				flag = 0
			else:
				i_ += i
	return reli

def re2str(s):
	'''
	Translate a string containing regular expressions to a string
	where all regular expression indicators ({,},|) are removed.
	'''
	s = s.replace('{','')
	s = s.replace('}','')
	s = s.replace('|','')
	return s

