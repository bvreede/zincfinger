

### DEFINE PARAMETERS AND INPUT FILES ###
dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"
files = []
for i in range (1,26):
	f = 'dmel-cluster' + str(i)
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

