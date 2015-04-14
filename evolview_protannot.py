#!/usr/bin/python

'''
This script can be used to interpret the results from 'findmotif.py', and
put these into a format that can be read by evolview for annotation
as 'protein domains'.
Optional is the customization of the colours; this option is not robust
to different (higher) numbers of domains, but can be easily adjusted.
However, a robust alternative is available by not opting in to the customized
colour palet.

Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 14 April 2015
'''

import csv,sys,math
from random import shuffle

if len(sys.argv) <= 1:
	sys.exit("USAGE: evolview_protannot.py path/to/inputfile [C]\nInputfile is csv file with motifhits per protein.\nC is not necessary, but if if is indicated the script will use a customized colour palette for the individual motifs. Make sure it is encoded properly in the script!")

infile = sys.argv[1]
#resfolder = ('/').join(infile.split('/')[:-1])
#outname = infile.split('/')[-1][:-4] + '.txt'
outname = infile[:-4] + '.txt'


# load source file, open results file
source = csv.reader(open(infile))
result = open(outname, "w")

if sys.argv[2] == 'C':
	custom = 1
else:
	custom = 0

custompalet = ['#b62020','#ff2020','#fc913a','#f9d62e','#fff797','#c6c386',
'#dff79e','#8ae429','#4d7f17','#204c39','#6b4423',
'#272d70','#0392cf','#83adb5','#66bbae',
'#bfb5b2','#a200ff','#6c2a7d','#ff99cc','#ff5588']



# define colours of motifs
'''
Function to translate numbers into a hex colour.
Input required: the total number of colours needed. Returns
a list of colours as long as (or longer) than the number.
'''
def getColour(maxcol):
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
		if custom == 1:
			colours = custompalet
		else:
			colours = getColour(len(line)-4)
		for n,m in enumerate(line[4:]):
			if len(m) > 0:
				CC,CH,HH = m.split('_')
				mcolour = colours[n] #findcolour(CC,CH,HH)
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
