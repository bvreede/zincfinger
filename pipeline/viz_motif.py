#!/usr/bin/python
'''
This script can be used to visualize the motif hits in the resulting database
from a motif search performed by 'findmotif.py', or similar looking databases.
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

import random,csv,sys,os,math

print "USAGE: viz_motif.py [C] (C is not necessary, but if it is indicated the script will use a customized colour palette for the individual motifs. Make sure it is encoded properly in the script!)\n\nNB: this is not an error message; the script is running, don't worry."

if sys.argv[1] == 'C':
	custom = 1
else:
	custom = 0

custompalet = ['#b62020','#ff2020','#fc913a','#f9d62e','#fff797','#c6c386',
'#dff79e','#8ae429','#4d7f17','#204c39','#6b4423',
'#272d70','#0392cf','#83adb5','#66bbae',
'#bfb5b2','#a200ff','#6c2a7d','#ff99cc','#ff5588']


'''
Defining input/output folders
'''
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/images"
infolder = "tovisualize"
outfolder = "visualized"

'''
Verifying/creating input/output files/folders
'''
if not os.path.exists(dbfolder):
	sys.exit("Could not find database directory. Verify path in the code.")
if not os.path.exists("%s/%s" %(dbfolder,infolder)):
	sys.exit("Could not find input directory. Verify path in the code.")
if not os.path.exists("%s/%s" %(dbfolder,outfolder)):
	os.system("mkdir %s/%s" %(dbfolder,outfolder))

files = []
for filename in os.listdir("%s/%s" %(dbfolder,infolder)):
	if filename[-4:] == ".csv":
		files.append(filename)

'''
Defining image parameters
'''
image_width = 4200
image_height = 28000
label_align = 400 #x axis of gene labels
linedist = 40 #distance between different lines

motlen = {} #dictionary of motif lengths
motclr = {} #dictionary of motif colours
motcount = {} #keeps track of how often a single motif is used


'''
Function to translate a number into a hex colour.
Input required: the total number of colours needed.
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
	#combine the r, g, b lists to generate hex colours
	for X,red in enumerate(CR):
		colour = '#' + red + CG[X] + CB[X]
		colours.append(colour)
	random.shuffle(colours)
	return colours



'''
define visualization of motifs/genes/labels
'''
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
determine the order at which the columns should be read. Also
fills up information on the colour and length dictionaries.
'''
def order(columns):
	length,callseq = [],[]
	if custom == 0:
		colours = getColour(len(columns))
	else:
		colours = custompalet
	for k,c in enumerate(columns):
		if c not in motlen:		# first, some household stuff: this is a new motif, so assign length + colour
			if len(c) == 0:
				continue
			CC,CH,HH = c.split('_')
			clen = int(CC)+int(CH)+int(HH)+4		# calculate the length for this motif
			motlen[c] = clen				# put the length for the motif in the dictionary
			mclr = colours[k]
			motclr[c] = mclr
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
	db = csv.reader(open("%s/%s/%s" %(dbfolder,infolder,f)))
	out = open("%s/%s/%s_viz.svg" %(dbfolder,outfolder,f[:-4]), "w")
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
	if linecount > 2:
		if linecount * linedist > image_height:
			print "Your image is incomplete: sequences for " + f + " extend to height = " + str(linecount * linedist)
		if max(maxlen) + label_align > image_width:
			print "Your image is incomplete: sequences for " + f + " extend to width = " + str(max(maxlen) + label_align)
	out.write('</svg>')
	out.close()
