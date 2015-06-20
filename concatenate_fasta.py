#!/bin/python
import config,csv


#give a name to the new file
newcat = "arth"
namechange = 1 # set to 0 if no name change is required.
orthsource = 'dmel' # the comparison species for ortholog determination; if namechange is 0, then it won't be used.

if namechange == 1:
	nc = "NC"
else:	
	nc = ""

#which files need to be included
#spp = ['drer','ggal','hsap','mmus','xtro'] #chor
spp = ['dmel','tcas','isca','dpul','smar','turt'] #arth
#spp = ['nvec','mlei','tadh','sman','aque'] #eani
#spp = ['atri','ppat','crei','atha','slyc','osat','smoe','cmer'] #plan
#spp = ['pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam'] #prot
#spp = ['hrob','spur','lgig','bmal','cele','drer','ggal','hsap','mmus','xtro','dmel','tcas','isca','dpul','smar','turt','nvec','mlei','tadh','sman','aque'] #anim
#spp = ['hrob','spur','lgig','bmal','cele','drer','ggal','hsap','mmus','xtro','dmel','tcas','isca','dpul','smar','turt','nvec','mlei','tadh','sman','aque', 'atri','ppat','crei','atha','slyc','osat','smoe','cmer','pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam'] #all species / alls


newfile = "%s/%s/150602-SM00355-%s%s_seq.fasta" %(config.mainfolder,config.seqfolder,newcat,nc)
nfile = open(newfile, "w")

filelist = ["%s/%s/150602-SM00355-%s_seq.fasta" %(config.mainfolder,config.seqfolder,sp) for sp in spp]


if namechange == 1:
	# first, make a dictionary to convert protein ID to name of the reference species
	prot2name = {}
	infilep2n = open("%s/%s/150602-SM00355-%s_seq.fasta" %(config.mainfolder,config.seqfolder,orthsource))
	for line in infilep2n:
		if line[0] == '>':
			# save protein: name in the prot2name dictionary
			header = line.strip().split('|')
			prot2name[header[2]] = header[1]
	# then, for each species (minus the reference species), extract ortholog information
	ortholog = {}
	for sp in spp:
		orthfile = csv.reader(open("%s/%s/%s-allorth.csv" %(config.mainfolder,config.orthfolder,sp)))
		for line in orthfile:
			if line[2] == orthsource: #only collect ortholog info for the reference species
				try:
					ortholog[line[1]] = prot2name[line[3]] #assign the name corresponding to the orthologous protein to the protein of the species currently in review.
				except KeyError:
					continue

for n,f in enumerate(filelist):
	fo = open(f)
	sp = spp[n]
	for line in fo:
		if namechange == 1 and sp != orthsource:
			if line[0] == ">":
				header = line.strip().split('|')
				try:
					orth = ortholog[header[2]]
					orthhead = "/%s-%s" %(orthsource[:2],orth)
				except KeyError:
					orthhead = ''
				if len(header[1]) > 0:
					oldname = "/%s" %header[1]
				else:
					oldname = header[1]
				newname = "%s%s%s" %(sp,oldname,orthhead)
				line = "%s|%s|%s\n" %(header[0],newname,header[2])
		nfile.write(line)
	nfile.write('\n')

nfile.close()
