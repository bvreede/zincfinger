import re,math
import matplotlib.pyplot as plt
import numpy as np
from random import shuffle

# path info
mainfolder = "/home/barbara/Dropbox/zinc_finger_data/newdata"
seqfolder = "sequences"
resfolder = "results"
imgfolder = "images"
dbfolder = "databases"
evfolder = "evolview"
orthfolder = "orthologs"


# define motifs
# a complete dataset
moCC = [1,2,3,4,5,6] #distances between CC
moCH = [7,8,9,10,11,12,13,14,15,16,17] #distances between CH
moHH = [1,2,3,4,5,6] #distances between HH

# turn this on if you want to search for all possible combinations of the above: (and don't forget to turn the custom list off!)
#motiflist = ['%s_%s_%s' %(m,n,o) for m in moCC for n in moCH for o in moHH]
# turn this on for custom motif list
#motiflist = ['2_7_4','2_8_3','2_9_3','2_10_5','2_11_3','2_11_4','2_12_2','2_12_3','2_12_4','2_12_5','2_12_6','2_13_3','2_13_4','2_14_3','2_14_4','2_15_4','3_8_3','4_12_3','4_12_4','4_15_3']

motiflist = ['1_12_3', '1_12_6', '1_7_3', '2_10_1', '2_11_3', '2_11_4', '2_12_2', '2_12_3', '2_12_4', '2_12_5', '2_12_6', '2_13_2', '2_13_3', '2_13_4', '2_14_3', '2_14_4', '2_15_4', '2_17_4', '2_7_4', '2_8_3', '2_9_3', '3_12_3', '3_12_4', '3_8_3', '4_12_3', '4_12_4', '4_12_6', '4_15_3', '5_14_3', '5_15_3', '5_15_4', '5_16_2', '5_7_6', '6_12_3', '6_14_6', '6_15_3', '6_15_4', '6_15_5', '6_17_4']


#distances before and after each C/H
plink,alink = 0,0



#individual groups
chor = ['drer','ggal','hsap','mmus','xtro'] #chor
arth = ['dmel','tcas','isca','dpul','smar','turt'] #arth
eani = ['nvec','mlei','tadh','sman','aque'] #eani
plan = ['atri','ppat','crei','atha','slyc','osat','smoe','cmer'] #plan
prot = ['pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam'] #prot
anim = ['hrob','spur','lgig','bmal','cele','drer','ggal','hsap','mmus','xtro','dmel','tcas','isca','dpul','smar','turt','nvec','mlei','tadh','sman','aque']

#all species
sppall = ['pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam','dmel','tcas','isca','dpul','smar','turt','atri','atha','crei','cmer','osat','ppat','smoe','slyc','aque','bmal','cele','hrob','lgig','mlei','nvec','sman','spur','tadh','drer','ggal','hsap','mmus','xtro']



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

motiflength,motifdict = make_motif_dict(motiflist)

#short alphabet
alphabet = '''ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_+=[];\<,>.?/'''

#long alphabet
alphabet = '''~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`'''
translationdict = {motif: alphabet[a] for a,motif in enumerate(motiflist)} # dictionary of motifs and the corresponding string element


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

plt.rcParams['xtick.labelsize']=28

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
