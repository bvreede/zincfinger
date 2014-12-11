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

image_width = 4200
image_height = 28000
label_align = 400 #x axis of gene labels
linedist = 40 #distance between different lines

### COLOURS AND LENGTS OF MOTIFS ###
motlen = {}
motclr = {}
motcount = {}
# define all possible colours (combinations of 'indvhex', except for greyscale)
clralphabet = []
indvhex = ['0','8','f']#['0','4','8','c','f']
for p in range(len(indvhex)):
	for q in range(len(indvhex)):
		for r in range(len(indvhex)):
			clr = '#' + indvhex[p] + indvhex[q] + indvhex[r]
			#if sum([p,q,r]) > len(indvhex)+1:
			if [p,q,r].count(len(indvhex)-1) > 1: # we don't like no pastels
				continue
			else:
				clralphabet.append(clr)


### DEFINE VISUALIZATION OF MOTIFS ###
def draw_arrow(m,X,Y,out):
	block = '<rect style="fill:%s;fill-opacity:1;stroke:none" width="%s" height="16" x="%s" y="%s" />' %(motclr[m],motlen[m],X,Y)
	arrowhead = '7,-8 -7,-8'
	Xa = X + motlen[m]
	arrow = '<path style="fill:%s;fill-opacity:1;stroke:none" d="m %s,%s 0,16 %s z" />' %(motclr[m],Xa,Y,arrowhead)
	motcount[m] += 1
	out.write("%s\n%s\n" %(block,arrow))

def draw_gene(name,seqlen,Y,out):
	Yn = Y+14
	Yg = Y+8
	labelX = label_align - 10
	label = '<text><tspan x="%s" y="%s" style="font-size:14px;fill:#000;fill-opacity:1;font-family:Helvetica;-inkscape-font-specification:Sans;text-align:end;text-anchor:end">%s</tspan></text>' %(labelX,Yn,name)
	gene = '<path style = "fill:none;stroke:#000;stroke-width:2;stroke-linecap:butt;stroke-opacity:1" d="m %s,%s %s,0" />' %(label_align,Yg,seqlen)
	out.write("%s\n%s\n" %(gene,label))	

def draw_legend(m,i,out):
	Y = 20
	Ylabel = Y + 30
	X = i*70+label_align
	Xlabel = X - 5
	draw_arrow(m,X,Y,out)
	label = '<text><tspan x="%s" y="%s" style="font-size:12px;fill:#000;fill-opacity:1;font-family:Helvetica;-inkscape-font-specification:Sans">%s</tspan></text>' %(Xlabel,Ylabel,m)
	out.write("%s\n" %(label))

'''
Makes a hex colour based on the motif sequence.
AWFUL AWFUL AWFUL SYSTEM NEEDS TO BE FIXED.
'''
def findcolour(CC,CH,HH):
	CCdict = {'2':'2','4':'e'}
	CHdict = {'8':'4','12':'a','15':'e'}
	HHdict = {'3':'0','4':'a','5':'e'}
	mclr = CCdict[CC] + CHdict[CH] + HHdict[HH]
	return mclr

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
			mclr = findcolour(CC,CH,HH)
			#ri = random.randint(0,len(clralphabet)-1)	# pick random colour from clralphabet
			motclr[c] = mclr# clralphabet[ri]			# assign it to that motif in the dictionary
			#clralphabet.pop(ri)				# remove the colour from the alphabet
			motcount[c] = -1
		length.append(motlen[c])
	for l in length:			# define the order in which columns need to be assessed, by long to short motifs
		i = length.index(max(length))		# get the index of the max value in length
		callseq.append(i)
		length[i] = 0				# change the value in length to 0
	return callseq


'''
Go through the files and make a visualization for each.
'''
for f in files:
	db = csv.reader(open("%s/%s_%s.csv" %(dbfolder,nameinfile,f)))
	out = open("%s/%s_%s.svg" %(dbfolder,nameoutfile,f), "w")
	out.write('<svg width="%s" height="%s" id="svg2" xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink">\n' %(image_width,image_height))
	linecount = 2
	maxlen = [] # to calculate the max width for the image
	start = 4 # the index of the first motif info
	for i,line in enumerate(db):
		if i == 0: # this is the header
			# the following is because of the possibly empty last element of header:
			if line[-1] == '':
				columns = line[start:-1]
			else:
				columns = line[start:]
			callseq = order(columns)
			for p,motif in enumerate(columns):
				draw_legend(motif,p,out)
		else:
			linecount+=1
			name = line[0] + '|' + line[1] + '|' + line[2]
			seqlen = int(line[3])
			maxlen.append(seqlen) # collecting the gene lengths to ensure they are smaller than image
			mdata = line[start:]
			# draw the gene
			Y = linecount * linedist
			draw_gene(name,seqlen,Y,out)
			# draw motifs based on the order in callseq
			for cs in callseq:
				motifname = columns[cs]
				if len(mdata[cs]) == 0:
					continue
				hits = mdata[cs].split('|')
				for h in hits:
					X = int(h) + label_align # the border where gene starts
					draw_arrow(motifname,X,Y,out)
	if linecount * linedist > image_height:
		print "Your image is incomplete: sequences for " + f + " extend to height = " + str(linecount * linedist)
	if max(maxlen) + label_align > image_width:
		print "Your image is incomplete: sequences for " + f + " extend to width = " + str(max(maxlen) + label_align)
	out.write('</svg>')
	out.close()

print motcount,motclr
