import os,math
from random import shuffle

indb = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/data_from_seetharam/fasta_files_seetharam"
dmelfasta = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-dmel_seq.fasta"
iscafasta = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-isca_seq-trans.fasta"
dpulfasta = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-dpul_seq-trans.fasta"
tcasfasta = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-tcas_seq-trans.fasta"



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
Function to blast the protein to a database and return the first hit.
'''
def blast(filename,spp):
	blastdb = "/home/barbara/data/zincfingers/%szfsTRANS.fa" %spp
	query = "%s/%s" %(indb,filename)
	out = "%s/tempout.txt" %(indb)
	blast = "blastp -query %s -db %s -out %s" %(query,blastdb,out)
	os.system(blast)
	readout = open(out)
	for line in readout:
		if line[0] == '>':
			line = line.strip()
			hit  = line[2:]
			break
	os.remove(out)
	hit = hit.replace('-111-','(')
	hit = hit.replace('-222-',')')
	hit = hit.replace('dm:', 'dm-')
	a,b,c = hit.split('-333-')
	hit = b + '|' + c
	return hit


#proteins to class
DMprot2class = {}
ISprot2class = {}
TCprot2class = {}
DPprot2class = {}
classes = []
countdm,countis,counttc,countdp = 0,0,0,0
for filename in os.listdir(indb):
	nameli = filename.split('-')
	if len(nameli) > 1:
		classes.append(nameli[1])
	if filename[0:23] == "Drosophila_melanogaster":
		prot = nameli[2].split('.')[0]
		DMprot2class[prot] = nameli[1]
		countdm +=1
	elif filename[0:3] == "Ixo":
		prot = nameli[2].split('.')[0]
		ISprot2class[prot] = nameli[1]
		countis +=1
	elif filename[0:4] == "Trib":
		prot = blast(filename,'tcas')
		TCprot2class[prot] = nameli[1]
		counttc +=1
	elif filename[0:3] == "Dap":
		prot = blast(filename,'dpul')
		DPprot2class[prot] = nameli[1]
		countdp +=1


#class to color
classes = list(set(classes))
colours  = getColour(len(classes))
ccdict = {}
for n,cl in enumerate(classes):
	ccdict[cl] = colours[n]


#Drosophila only: genes to class
dmfasta = open(dmelfasta)
DMgene2class = {}
genenames = []
for line in dmfasta:
	if line[0] == '>':
		line = line.strip()
		a,b = line.split('|')[1:] #gene name in a; prot name in b
		#if prot name in prot2class, save gene name with same class
		if b in DMprot2class:
			DMgene2class[a] = DMprot2class[b]
		#also, save all names in a list, as they will appear in the evolview visualization
		name = a + '|' + b
		genenames.append(name)

#Ixodes: 
isfasta = open(iscafasta)
ISname2class = {}
for line in isfasta:
	if line[0] == '>':
		line = line.strip()
		a,b = line.split('|')[1:] #gene name in a; prot name in b
		#save all names in a list, as they will appear in the evolview visualization
		name = a + '|' + b
		name = name.replace('dm:', 'dm-')
		genenames.append(name)
		#if prot name in prot2class, save name with same class
		c = b.replace('-PA','')
		if b in ISprot2class:
			ISname2class[name] = ISprot2class[b]
		elif c in ISprot2class:
			ISname2class[name] = ISprot2class[c]

#Tribolium and Daphnia: 
tcfasta = open(tcasfasta)
for line in tcfasta:
	if line[0] == '>':
		line = line.strip()
		a,b = line.split('|')[1:] #gene name in a; prot name in b
		#save all names in a list, as they will appear in the evolview visualization
		name = a + '|' + b
		name = name.replace('dm:', 'dm-')
		genenames.append(name)

dpfasta = open(dpulfasta)
for line in dpfasta:
	if line[0] == '>':
		line = line.strip()
		a,b = line.split('|')[1:] #gene name in a; prot name in b
		#save all names in a list, as they will appear in the evolview visualization
		name = a + '|' + b
		name = name.replace('dm:', 'dm-')
		genenames.append(name)


print countis,counttc,countdp

# generate the evolview-readable doc
evolview = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/seetharam_evolview_labels.txt", "w")
evolview.write(" ## leaf background color\n\n")
for gene in genenames:
	gn = gene.split('|')[0]
	if gn in DMgene2class:
		evolview.write("%s\t%s\tprefix\n" %(gene,ccdict[DMgene2class[gn]]))
	elif gene in ISname2class:
		evolview.write("%s\t%s\tprefix\n" %(gene,ccdict[ISname2class[gene]]))
		countis -= 1
	elif gene in TCprot2class:
		evolview.write("%s\t%s\tprefix\n" %(gene,ccdict[TCprot2class[gene]]))
		counttc -= 1
	elif gene in DPprot2class:
		evolview.write("%s\t%s\tprefix\n" %(gene,ccdict[DPprot2class[gene]]))
		countdp -= 1
evolview.close()

for key in DPprot2class:
	if key not in genenames:
		print key, DPprot2class[key]
	
for key in TCprot2class:
	if key not in genenames:
		print key, TCprot2class[key]

print countis,counttc,countdp

print len(ISname2class),len(TCprot2class),len(DPprot2class)
