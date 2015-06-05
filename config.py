import re

# path info
mainfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/newdata"
seqfolder = "sequences"
resfolder = "results"
imgfolder = "images"
dbfolder = "databases"
evfolder = "evolview"


# define motifs
# a complete dataset
moCC = [2,3,4] #distances between CC
moCH = [7,8,9,10,11,12,13,14,15] #distances between CH
moHH = [2,3,4,5,6] #distances between HH

# turn this on if you want to search for all possible combinations of the above: (and don't forget to turn the custom list off!)
motiflist = ['%s_%s_%s' %(m,n,o) for m in moCC for n in moCH for o in moHH]
# turn this on for custom motif list
#motiflist = ['2_7_4','2_8_3','2_9_3','2_10_5','2_11_3','2_11_4','2_12_2','2_12_3','2_12_4','2_12_5','2_12_6','2_13_3','2_13_4','2_14_3','2_14_4','2_15_4','3_8_3','4_12_3','4_12_4','4_15_3']


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

alphabet = '''ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`'''
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


