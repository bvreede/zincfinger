'''
This script can be used to visualize the motif hits in the resulting database
from a motif search performed by 'findmotif_RE.py', or similar looking databases.
The database has the following columns:
Gene_stable_ID, Gene_name, Protein_stable_ID, Sequence_length,domain1,domain2,...,domainn,
	(where domain1 etc. are motif names according to the format P_Q_R, where P, Q, R, are
	integers describing the number of residues between C and C, C and H, H and H, respectively
	in a C2H2 zinc finger motif.)
The column headers are present in the database (and required for this script to function).
The product from this script is an .svg image that can be read in a browser or opened with e.g.
inkscape (we love inkscape!).
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 11 December 2014
'''

import random,csv

### DEFINE PARAMETERS AND INPUT FILES ###
# infile: dbfolder/nameinfile_files.csv
# outfile: dbfolder/nameoutfile_files.svg
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results"
nameinfile = "motifhits"
nameoutfile = "motifviz"
files = ['dmel','tcas','smar','dpul','isca']

### COLOURS AND LENGTS OF MOTIFS ###
motlen = {}
motclr = {}
# define all possible colours (combinations of 'indvhex', except for greyscale)
clralphabet = []
indvhex = ['0','4','8','c','f']
for p in range(len(indvhex)):
	for q in range(len(indvhex)):
		for r in range(len(indvhex)):
			clr = '#' + indvhex[p] + indvhex[q] + indvhex[r]
			if p == q == r: # we don't like no greyscale
				continue
			else:
				clralphabet.append(clr)
			

'''
determine the order at which the columns should be read. Also
fills up information on the colour and length dictionaries.
'''
def order(columns):
	length,callseq = [],[]
	for c in columns:
		if c not in motlen:		# first, some household stuff: this is a new motif, so assign length + colour
			CC,CH,HH = c.split('_')
			clen = int(CC)+int(CH)+int(HH)+4		# calculate the length for this motif
			motlen[c] = clen				# put the length for the motif in the dictionary
			ri = random.randint(0,len(clralphabet)-1)	# pick random colour from clralphabet
			motclr[c] = clralphabet[ri]			# assign it to that motif in the dictionary
			clralphabet.pop(ri)				# remove the colour from the alphabet
		length.append(motlen[c])
	for l in length:			# define the order in which columns need to be assessed, by long to short motifs
		i = length.index(max(length))		# get the index of the max value in length
		callseq.append(i)
		length[i] = 0				# change the value in length to 0
	return callseq

for f in files:
	db = csv.reader(open("%s/%s_%s.csv" %(dbfolder,nameinfile,f)))
	for i,line in enumerate(db):
		if i == 0:
			# the following is because of the possibly empty last element of header:
			if line[-1] == '':
				columns = line[4:-1]
			else:
				columns = line[4:]
			callseq = order(columns)



'''
actheader = ['Gene_stable_ID', 'Gene_name', 'Protein_stable_ID', 'Pfam_ID', 'SMART_ID', 'GO_domain', 'GO_term_accession', 'GO_term_name']
headseq = []
headline = ""
for name in actheader:
	try:
		i = curheader.index(name) #finds the ACTUAL header element in the CURRENT header and returns the index
	except ValueError:
		print "%s column not found in database. Continuing anyway..." %name
		i = 'N'
	headseq.append(i)
	headline += name + ',' # write the actual header for the final database
o.write("%s\n" %headline[:-1]) # and write it to file

'''


# determine the order in which to draw: large domains first
# so read first line, get indices from the largest one up to smaller ones


"""


for i in range (1,38):
	f = 'dmel-G-cluster' + str(i)
	files.append(f)
#files = ['dmel1','tcas','smar','dpul','isca']
image_width = 4000
image_height = 20000
label_align = 200 #x axis of label
linedist = 40 #distance between different lines
motifs_for_legend = ['2_8_3','2_12_3','2_12_4','2_12_5','4_12_3','4_12_4','4_15_3']



### DEFINE VISUALIZATION OF MOTIFS ###
colordict = {'2_8_3': '#c55', '2_12_3': '#5c5', '2_12_4': '#55c', '2_12_5': '#cc5','4_12_3': '#c5c','4_12_4': '#5cc','4_15_3': '#a60'}

def draw_arrow(motif,L,drx,X,Y):
	block = '<rect style="fill:%s;fill-opacity:1;stroke:none" width="%s" height="16" x="%s" y="%s" />' %(colordict[motif],L,X,Y)
	if drx == 'fwd':
		arrowhead = '8,-8 -8,-8'
		Xa = X+L
	elif drx == 'rev':
		arrowhead = '-8,-8 8,-8'
		Xa = X
	arrow = '<path style="fill:%s;fill-opacity:1;stroke:none" d="m %s,%s 0,16 %s z" />' %(colordict[motif],Xa,Y,arrowhead)
	return block,arrow

def draw_gene(name,seqlen,Y):
	Yn = Y+14
	Yg = Y+8
	labelX = label_align - 10
	label = '<text><tspan x="%s" y="%s" style="font-size:14px;fill:#000;fill-opacity:1;font-family:Helvetica;-inkscape-font-specification:Sans;text-align:end;text-anchor:end">%s</tspan></text>' %(labelX,Yn,name)
	gene = '<path style = "fill:none;stroke:#000;stroke-width:2;stroke-linecap:butt;stroke-opacity:1" d="m %s,%s %s,0" />' %(label_align,Yg,seqlen)
	return gene,label


### PREPARING DRAW, CALCULATE PARAMETERS. DRAW MOTIF ###
def draw(motif,p,Y,outputimg):
	Lm = motif.split('_')
	L = int(Lm[0])+int(Lm[1])+int(Lm[2])+4 #size of the motif: four residues plus nucleotides between
	X = label_align + p
	drx = 'fwd' #default direction
	if motif[-4:] == "_inv":
		motif = motif[:-4]
		drx = 'rev'
	block,arrow = draw_arrow(motif,L,drx,X,Y)
	outputimg.write("%s\n%s\n" %(block,arrow))

### MAKING THE LEGEND ###
def draw_legend(leg_mot,i,outputimg):
	Y = 20
	Ylabel = Y + 30
	X = i*70+label_align
	Xlabel = X - 5
	block,arrow = draw_arrow(leg_mot,20,'fwd',X,Y)
	label = '<text><tspan x="%s" y="%s" style="font-size:12px;fill:#000;fill-opacity:1;font-family:Helvetica;-inkscape-font-specification:Sans">%s</tspan></text>' %(Xlabel,Ylabel,leg_mot)
	outputimg.write("%s\n%s\n%s\n" %(block,arrow,label))

### GO THROUGH ALL POSSIBLE CATEGORIES (SEE LIST) ###
for s in files:
	inputdb = open("%s/results/motifhits_%s.csv" %(dbfolder,s))
	outputimg = open("%s/results/motifimg_%s.svg" %(dbfolder,s), "w")
	outputimg.write('<svg width="%s" height="%s" id="svg2" xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink">\n' %(image_width,image_height))
	for i in range(7): # making the legend
		leg_mot = motifs_for_legend[i]
		draw_legend(leg_mot,i,outputimg)
	linecount = 0
	maxlen = [] # to calculate the max width for the image
	for line in inputdb:
		linecount += 1
		l = line.strip().split(',')[:-1] #remove last empty item (due to trailing comma in csv)
		if l[0] == "Gene_stable_ID": #this is the header
			motiflist = l[4:] #collects the motifs in the header line
		else:
			Y = linecount * linedist
			name,seqlen = l[1],l[3]
			if name == "":
				name=l[0]
			maxlen.append(int(seqlen))
			gene,label = draw_gene(name,seqlen,Y)
			outputimg.write("%s\n%s\n" %(gene,label))
			motifpos = l[4:] #list of positions for a certain motif.
			for k in range(len(motifpos)): #for each motif:
				if len(motifpos[k]) > 0: #if the motif actually has hits in this gene
					positions = motifpos[k].split('|')
					for p in positions:
						p = int(p)
						draw(motiflist[k],p,Y,outputimg)
	outputimg.write('</svg>')
	if linecount * linedist > image_height:
		print "Your image is incomplete: sequences for " + s + " extend to height = " + str(linecount * linedist)
	if max(maxlen) + label_align > image_width:
		print "Your image is incomplete: sequences for " + s + " extend to width = " + str(max(maxlen) + label_align)
"""
