import os,math
from random import shuffle

indb = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/data_from_seetharam/fasta_files_seetharam"
dmelfasta = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-dmel_seq.fasta"

prot2class = {}
classes = []
for filename in os.listdir(indb):
	if filename[0:23] == "Drosophila_melanogaster":
		nameli = filename.split('-')
		prot = nameli[2].split('.')[0]
		prot2class[prot] = nameli[1]
		classes.append(nameli[1])



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



#give classes a colour
classes = list(set(classes))
colours  = getColour(len(classes))
ccdict = {}
for n,cl in enumerate(classes):
	ccdict[cl] = colours[n]

#prot to evolview name
fasta = open(dmelfasta)
namedict = {}
for line in fasta:
	if line[0] == '>':
		line = line.strip()
		a,b = line.split('|')[1:]
		name = a + '|' + b
		namedict[b] = name


# generate the evolview-readable doc
evolview = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/seetharam_evolview_labels.txt", "w")
evolview.write(" ## leaf background color\n\n")
for key in prot2class:
	evolview.write("%s\t%s\tprefix\n" %(namedict[key],ccdict[prot2class[key]]))
evolview.close()


