import csv,sys

if len(sys.argv) <= 1:
	sys.exit("USAGE: python evolview_protannot.py path/to/inputfile")

infile = sys.argv[1]
#resfolder = ('/').join(infile.split('/')[:-1])
#outname = infile.split('/')[-1][:-4] + '.txt'
outname = infile[:-4] + '.txt'


# load source file, open results file
source = csv.reader(open(infile))
result = open(outname, "w")

# define colours of motifs
'''
Makes a hex colour based on the motif sequence.
NOT ROBUST AGAINST OTHER COLOUR COMBINATIONS!
'''
def findcolour(CC,CH,HH):
	CCdict = {'2':'20','4':'D0'}
	CHdict = {'8':'40','12':'A0','15':'F0'}
	HHdict = {'3':'00','4':'A0','5':'E0'}
	try:
		mclr = '#' + CHdict[CH] + HHdict[HH] + CCdict[CC]
	except KeyError:
		sys.exit("I question your motifs: they don't exist in the colour scheme of this visualization module.")
	return mclr


'''
writes the first lines of the results file
'''
def headerprint(colourdict):
	domains, colours = '',''
	for d in colourdict:
		domains += d + ','
		colours += colourdict[d] + ','
	result.write('##\n!groups\t%s\n!colors\t%s\n!showLegends\t0\n##\n' %(domains[:-1],colours[:-1]))

# read the first line and determine the sequence of all motifs, and the colour of the motifs
colourdict,msequence,mlen = {},{},{}
for i,line in enumerate(source):
	if i == 0: #only on the first (header) line
		for n,m in enumerate(line[4:]):
			if len(m) > 0:
				CC,CH,HH = m.split('_')
				mcolour = findcolour(CC,CH,HH)
				# add colour, sequence, motif length to dictionaries
				colourdict[m] = mcolour
				msequence[m] = n + 4
				mlen[m] = int(CC)+int(CH)+int(HH)+4
		# print header
		headerprint(colourdict)
		continue
	# for every line following that, read the information and generate a line in the output file
	#name tab length tab begin,end,motif tab begin,end,motif tab
	annot = line[1] + '|' + line[2] + '\t' + line[3] + '\t' #name is line[2], length is line[3]
	for m in msequence:
		hits = line[msequence[m]]
		if len(hits) > 0:
			hits = hits.split('|')
			for b in hits:
				e = int(b)+mlen[m]
				annot += b + ',' + str(e) + ',' + m + '\t'
	result.write('%s\n' %annot)

# close the output file
result.close()
