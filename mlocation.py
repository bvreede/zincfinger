#!/usr/bin/python

'''
This script can be used to determine where in a connected series of motifs
certain motifs are found (i.e.: are they found at the start, middle, or
end of a series). The input for this needs to be the result of a findmotif.py
run on a protein sequence.

Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 10 August 2015
'''

import config,sys,os
import numpy as np
import matplotlib.pyplot as plt

#settings for images
plt.rcParams['xtick.labelsize']=11

#dictionary of species sets to translate user input to a list in the config file
setdict = {'sppall': config.sppall,'spp700':config.spp700, 'chor': config.chor,'arth': config.arth, 'eani': config.eani, 'plan':config.plan, 'prot': config.prot,'anim':config.anim}

#check input given: input is needed, and it needs to exist in the dictionary of species sets
exit_message = "USAGE: mlocation.py species-set (optional e.g. arth for arthropodes, or sppall for all species in the database)\nNB: check that the species set exists in the mlocation.py dictionary."
if len(sys.argv) <= 1:
	sys.exit(exit_message)

try:
	sppset = setdict[sys.argv[1]]
except KeyError:
	sys.exit(exit_message)



##GENERATE INPUT FILE##
#concatenate the protstring motif sequences for the selected species
dbdir = config.mainfolder + "/" + config.seqmfolder
infile = dbdir + '/' + config.idr + '-' + sys.argv[1][:4] + '_protstring.fa'

infile_write = open(infile,"w")

for f in os.listdir(dbdir):
	#isolate the species identifier
	spp_idr = f.split(config.idr + '-')[1][:4]
	#confirm that this identifier exists in the list of species specified by the user
	if spp_idr not in sppset:
		continue
	#open the file and write its contents to the new infile
	readfasta = open(dbdir + '/' + f)
	for line in readfasta:
		infile_write.write(line)

infile_write.close()

###NAME OUTPUT FILES###
infilebrev = config.idr + '-' + sys.argv[1] #specifier used to name outputfiles
outtotalsname = "%s/%s/%s_NAMEtotals.csv" %(config.mainfolder,config.resfolder,infilebrev)
outfigure = "%s/%s/%s_NAME.svg" %(config.mainfolder,config.imgfolder,infilebrev)



### OPTION: only treat non-ambiguous zfs! ###
# set this to 0 if you want to treat all.
na = 1

# define the output range and labels for each analysis
CCdata = ['CC','C-C distance',range(1,8)] 
CHdata = ['CH','C-H distance',range(7,18)]
HHdata = ['HH','H-H distance',range(1,8)]


#output = range(7,18)


lengths_of_series = []
lengths_of_proteins = []

# start the dictionary containing the results
# results are collected per motif and for the total, and consist of a list of four elements:
# [total count] [found in first] [found in middle] [found in last]
resultdict = {m:[0,0,0,0] for m in config.motiflist2}
resultdict['total'] = [0,0,0,0]


# make dictionary out of the input file
infi = open(infile)
fadict = config.fastadicter(infi)


def re2str(s):
	'''
	Takes a regular expression and turns it into a string of letters only.
	'''
	s = s.replace('{','')
	s = s.replace('}','')
	s = s.replace('|','')
	return s

def removere(s):
	p = s.count('{')
	for k in range(p):
		start = s.index('{')
		end = s.index('}') + 1
		s = s[:start] + s[end:]
	return s

def re2list(s):
	'''
	Takes a string that contains regular expression elements, and returns
	a list with the first, middle, and last elements separated, and all
	regular expression characters removed.
	'''
	#identify first
	if s[0] == '{':
		ends1=s.index('}') + 1
		if na == 1:
			s1 = ''
		else:
			s1 = re2str(s[:ends1])
		s = s[ends1:] #remove first element from s
	else:
		s1= s[0]
		s = s[1:] #remove first element from s
	#identify last
	if s[-1] == '}':
		starts2 = s.rfind('{')
		if na == 1:
			s2 = ''
		else:
			s2 = re2str(s[starts2:])
		s = s[:starts2] #remove last element from s
	else:
		s2 = s[-1]
		s = s[:-1]  #remove last element from s
	# now remove re elements from the leftover middle
	if na == 1:
		s = removere(s)
	else:
		s = re2str(s)
	reli = [s1,s,s2]
	return reli

def first(aa):
	global resultdict
	for a in aa:
		m = config.translationdict_inv2[a]
		resultdict[m][0] += 1
		resultdict[m][1] += 1
		resultdict['total'][0] += 1
		resultdict['total'][1] += 1

def last(aa):
	global resultdict
	for a in aa:
		m = config.translationdict_inv2[a]
		resultdict[m][0] += 1
		resultdict[m][3] += 1
		resultdict['total'][0] += 1
		resultdict['total'][3] += 1

def middle(aa):
	global resultdict
	for a in aa:
		m = config.translationdict_inv2[a]
		resultdict[m][0] += 1
		resultdict[m][2] += 1
		resultdict['total'][0] += 1
		resultdict['total'][2] += 1






# start the analysis!
for protein in fadict:
	ssplit = fadict[protein].strip().split('Z')
	prot = protein.replace('Z','')
	prl = prot.count('{') + prot.count('}') + prot.count('|') * 2
	protl = len(prot) - prl
	lengths_of_proteins.append(protl)
	for s in ssplit:
		recount = s.count('{') + s.count('}') + s.count('|') * 2
		if len(s) - recount < 2:
			continue
		else:
			lengths_of_series.append(len(s) - s.count('{') - s.count('}') - s.count('|') * 2)
			sli = re2list(s)
			first(sli[0])
			last(sli[2])
			middle(sli[1])

percdict = {}
for r in resultdict:
	resli = resultdict[r]
	try:
		percli = [100*rl/float(resli[0]) for rl in resli[1:]]
		percli.append(resli[0])
	except ZeroDivisionError:
		continue
	percdict[r] = percli



def barchart_loc(aadata):
	'''
	Make bar charts with pairs of bars grouped for easy comparison.
	'''
	global outtotalsname
	global outfigure
	# now assemble the dataset:
	bardata = {}
	for n in range(18):
		cfi,cmi,cen,total = 0,0,0,0
		for p in percdict:
			if p == 'total':
				continue
			CC,CH,HH = p.split('_')
			if aadata[0] == 'CC':
				key = CC
			elif aadata[0] == 'CH':
				key = CH
			elif aadata[0] == 'HH':
				key = HH
			else:
				sys.exit("Calling barchart with the wrong info.")
			if int(key) == n:
				c1 = percdict[p][0]
				cm = percdict[p][1]
				c2 = percdict[p][2]
				ct = float(percdict[p][3])
				cfi += c1*ct
				cmi += cm*ct
				cen += c2*ct
				total += ct
		if total == 0:
			bardata[n] = [0,0,0,0]
		else:
			bardata[n] = [cfi/total,cmi/total,cen/total,total]

	# get data for plots:
	set1,set2,set3,totals = [],[],[],[]
	for n in aadata[2]:
		set1.append(bardata[n][0])
		set2.append(bardata[n][1])
		set3.append(bardata[n][2])
		totals.append(bardata[n][3])

	n_groups = len(set1)
	fig, ax = plt.subplots()
	index = np.arange(n_groups)
	bar_width = 0.25

	rects1 = plt.bar(index, set1, bar_width,alpha=0.7, color='c',label='First')
	rects2 = plt.bar(index + bar_width+0.04, set2, bar_width,alpha=0.8, color='g',label='Middle')
	rects3 = plt.bar(index + 2*bar_width+0.08, set3, bar_width,alpha=0.5,color='r',label='Last')

	plt.xlabel(aadata[1])
	plt.ylabel('Frequency')
	plt.title('Appearance of motif in series')
	plt.xticks(index + bar_width, (list(aadata[2])))
	plt.legend()

	plt.tight_layout()

	def autolabel(rects):
	    # attach some text labels
	    for n,rect in enumerate(rects):
	        height = rect.get_height()
	        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%totals[n], ha='center', va='bottom')

	#autolabel(rects1)
	#autolabel(rects2)

	name = aadata[0]
	outtotalsname_loc = outtotalsname.replace('NAME',name)
	outtotals = open(outtotalsname_loc, "w")
	
	for n,t in enumerate(totals):
		outtotals.write("%s,%s,%s,%s,%s\n" %(aadata[2][n],t,set1[n],set2[n],set3[n]))
	outtotals.close()

	outfigure_loc = outfigure.replace('NAME',name)
	plt.savefig(outfigure_loc)
	plt.clf()

barchart_loc(CCdata)
barchart_loc(CHdata)
barchart_loc(HHdata)

