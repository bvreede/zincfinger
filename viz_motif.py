

### DEFINE PARAMETERS AND INPUT FILES ###
dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"
files = ['dmel','tcas','smar','dpul','isca']
image_width = 500
image_height = 500



### DEFINE VISUALIZATION OF MOTIFS ###
colordict = {'2_8_3': '#c55', '2_12_3': '#5c5', '2_12_4': '#55c', '2_12_5': '#cc5','4_12_3': '#c5c','4_12_4': '#5cc','4_15_3': '#a60'}

def draw_arrow(motif,L,drx,X,Y):
	block = '<rect style="fill:%s;fill-opacity:1;stroke:none" width="%s" height="20" x="%s" y="%s" />' %(colordict[motif],L,X,Y)
	if drx == 'fwd':
		arrowhead = '10,-10 -10,-10'
		Xa = X+L
	elif drx == 'rev':
		arrowhead = '-10,-10 10,-10'
		Xa = X
	arrow = '<path style="fill:%s;fill-opacity:1;stroke:none" d="m %s,%s 0,20 %s z" />' %(colordict[motif],Xa,Y,arrowhead)
	return block,arrow

def draw_gene:
	

### PREPARING DRAW, CALCULATE PARAMETERS. DRAW MOTIF ###
def draw(motif,p,Y,outputimg):
	Lm = motif.split('_')
	L = int(Lm[0])+int(Lm[1])+int(Lm[2])+4 #size of the motif: four residues plus nucleotides between
	X = 
	Y = linecount * 50

	drx = 'fwd' #default direction
	if motif[-4:] == "_inv":
		motif = motif[:-4]
		drx = 'rev'
	block,arrow = draw_arrow(motif,L,drx,X,Y)
	outputimg.write("%s\n%s\n" %(block,arrow))


for s in files:
	inputdb = open("%s/results/motifhits_%s.csv" %(dbfolder,s))
	outputimg = open("%s/results/motifimg_%s.svg" %(dbfolder,s), "w")
	outputimg.write('<svg width="%s" height="%s" id="svg2" xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink">\n' %(image_width,image_height))
	linecount = 0
	for line in inputdb:
		linecount += 1
		l = line.strip().split(',')[:-1] #remove last empty item (due to trailing comma in csv)
		if l[0] == "Gene_stable_ID":
			motiflist = l[4:]
		else:
			Y = linecount * 50
			name,seqlen = l[1],l[3]
			gene,label = draw_gene(name,seqlen,linecount)
			outputimg.write("%s\n%s\n" %(gene,label))
			motifpos = l[4:] #list of positions for a certain motif.
			for k in range(len(motifpos)): #for each motif:
				if len(motifpos[k]) > 0: #if the motif actually has hits in this gene
					positions = motifpost[k].split('|')
					for p in positions:
						draw(motiflist[k],p,Y,outputimg)
	outputimg.write('</svg>')

