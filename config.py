import re
import matplotlib.pyplot as plt
import numpy as np

# path info
mainfolder = "/home/barbara/Dropbox/zinc_finger_data/newdata"
seqfolder = "sequences"
resfolder = "results"
imgfolder = "images"
dbfolder = "databases"
evfolder = "evolview"


# define motifs
# a complete dataset
moCC = [1,2,3,4,5,6] #distances between CC
moCH = [7,8,9,10,11,12,13,14,15,16,17] #distances between CH
moHH = [1,2,3,4,5,6] #distances between HH

# turn this on if you want to search for all possible combinations of the above: (and don't forget to turn the custom list off!)
motiflist = ['%s_%s_%s' %(m,n,o) for m in moCC for n in moCH for o in moHH]
# turn this on for custom motif list
#motiflist = ['2_7_4','2_8_3','2_9_3','2_10_5','2_11_3','2_11_4','2_12_2','2_12_3','2_12_4','2_12_5','2_12_6','2_13_3','2_13_4','2_14_3','2_14_4','2_15_4','3_8_3','4_12_3','4_12_4','4_15_3']

motiflist = ['1_13_3', '1_13_4', '1_15_1', '1_7_4', '1_8_3', '1_8_6', '2_10_1', '2_10_4', '2_11_3', '2_11_4', '2_11_5', '2_12_1', '2_12_2', '2_12_3', '2_12_4', '2_12_5', '2_12_6', '2_13_3', '2_13_4', '2_14_1', '2_14_2', '2_14_3', '2_14_4', '2_14_5', '2_15_3', '2_15_4', '2_15_5', '2_16_1', '2_16_2', '2_16_3', '2_16_4', '2_16_5', '2_16_6', '2_17_4', '2_7_4', '2_8_3', '2_8_4', '2_8_6', '2_9_2', '2_9_3', '2_9_4', '2_9_5', '2_9_6', '3_11_1', '3_11_3', '3_11_4', '3_12_4', '4_10_1', '4_10_3', '4_10_4', '4_10_5', '4_12_1', '4_12_2', '4_12_3', '4_12_4', '4_12_5', '4_12_6', '4_15_3', '5_10_5', '5_11_1', '5_12_3', '5_12_4', '5_14_5', '5_15_3', '5_15_4', '5_9_1', '5_9_3', '5_9_4', '5_9_5', '6_10_3', '6_12_3', '6_12_4', '6_12_5', '6_12_6', '6_14_6', '6_15_3', '6_15_4', '6_8_3', '6_8_4', '6_8_5']

#distances before and after each C/H
plink,alink = 0,0

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

alphabet = '''ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_+=[];\<,>.?/'''

#alphabet = '''~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`'''
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



