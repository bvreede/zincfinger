import os,csv

zffolder = "/home/barbara/data/zincfingers"
seqfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences"

"""
# replace pipes and other characters in the fasta files (do this before making dbs)
for fasta in os.listdir(zffolder):
	f = open("%s/%s" %(zffolder,fasta))
	f2 = open("%s/%s2" %(zffolder,fasta), "w")
	for line in f:
		line = line.replace('(','-111-')
		line = line.replace(')','-222-')
		line = line.replace('|','-333-')
		f2.write(line)
	f2.close()
	os.remove("%s/%s" %(zffolder,fasta))
	os.rename("%s/%s2" %(zffolder,fasta), "%s/%s" %(zffolder,fasta))


# prep databases by individually saving sequences
def prepdb(spp):
	infasta = open("%s/%szfs.fa" %(zffolder,spp))
	outfolder = "%s/%s" %(zffolder,spp)
	if os.path.exists(outfolder):
		print "Folder %s already exists." %spp
	else:
		os.system("mkdir %s" %outfolder)

	for line in infasta:
		line = line.strip()
		if len(line) == 0:
			ofile.write('\n')
			ofile.close()
		elif line[0] == '>':
			name = line[1:]
			ofile = open("%s/%s" %(outfolder,name), "w")
		else:
			ofile.write(line)

# prep blast db
def prepblastdb(spp):
	os.system("makeblastdb -in %s/%szfs.fa -dbtype prot" %(zffolder,spp))

# only necessary to do this once per species
for spp in ['dpul','dmel','tcas','turt','smar','isca']:
	prepdb(spp)
	prepblastdb(spp)


# blast each file of dmel against the five other species, and save the name of its top hit
# output: csv file with dmel in column 1, spp in column 2; named: dmel-spp
for spp in ['dpul','tcas','turt','smar','isca']:
	outdb = open("%s/dmel-%s.csv" %(zffolder,spp), "w")
	for f in os.listdir("%s/dmel" %(zffolder)):
		if f[0] == '.':
			print "skipping hidden file ", f
			continue
		blastcommand = "blastp -query %s/dmel/%s -db %s/%szfs.fa -out %s/tempout.txt" %(zffolder,f,zffolder,spp,zffolder)
		os.system(blastcommand)
		blastres = open("%s/tempout.txt" %zffolder)
		for line in blastres:
			if line[0] == '>':
				line = line.strip()
				outdb.write("%s,%s\n" %(f,line[1:]))
				break
	print "Blast from dmel to %s complete." %spp
	outdb.close()

# blast each file from the five species against dmel and save the name of its top hit
# output: csv file with dmel in column 1, spp in column 2; named: spp-dmel
for spp in ['dpul','tcas','turt','smar','isca']:
	outdb = open("%s/%s-dmel.csv" %(zffolder,spp), "w")
	for f in os.listdir("%s/%s" %(zffolder,spp)):
		if f[0] == '.':
			print "skipping hidden file ", f
			continue
		blastcommand = "blastp -query %s/%s/%s -db %s/dmelzfs.fa -out %s/tempout.txt" %(zffolder,spp,f,zffolder,zffolder)
		os.system(blastcommand)
		blastres = open("%s/tempout.txt" %zffolder)
		for line in blastres:
			if line[0] == '>':
				line = line.strip()
				outdb.write("%s,%s\n" %(line[1:],f))
				break
	print "Blast from %s to dmel complete." %spp
	outdb.close()
"""

# make a reciprocal check of each pair of csv files
gndict = open("%s/genenamedict.csv" %(seqfolder), "w")
for spp in ['dpul','tcas','turt','smar','isca']:
	todmel = csv.reader(open("%s/%s-dmel.csv" %(zffolder,spp)))
	todmeldx = {}
	for g in todmel:
		spg = g[1].strip()
		dmg = g[0].split('-333-')[1] #the middle argument is the gene name
		spg = spg.replace('-333-','|')
		spg = spg.replace('-111-','(')
		spg = spg.replace('-222-',')')
		dmg = dmg.replace('-111-','(')
		dmg = dmg.replace('-222-',')')
		todmeldx[spg] = dmg
	fromdmel = csv.reader(open("%s/dmel-%s.csv" %(zffolder,spp)))
	for g in fromdmel:
		spg = g[1].strip()
		dmg = g[0].split('-333-')[1] #the middle argument is the gene name
		spg = spg.replace('-333-','|')
		spg = spg.replace('-111-','(')
		spg = spg.replace('-222-',')')
		dmg = dmg.replace('-111-','(')
		dmg = dmg.replace('-222-',')')
		if todmeldx[spg] == dmg:
			gndict.write('%s,%s\n' %(spg,dmg))
	print "Made dictionary for %s." %spp
gndict.close()

# read the gene translation dictionary and rewrite the species fasta files
gtransdx = {}
gtrans = csv.reader(open("%s/genenamedict.csv" %(seqfolder)))
for line in gtrans:
	gtransdx[line[0]] = line[1]
 
for spp in ['dpul','tcas','turt','smar','isca']:
	fasta = open('%s/150111-SM00355-%s_seq.fasta' %(seqfolder,spp))
	fasta2 = open('%s/150111-SM00355-%s_seq-trans.fasta' %(seqfolder,spp), "w")
	for line in fasta:
		if line[0] == '>':
			totalID = line.strip()
			totalID = totalID[1:]
			if totalID in gtransdx:
				gID,gna,pID = line.split('|')
				line = '%s|%s(dm:%s)|%s' %(gID,gna,gtransdx[totalID],pID)
		fasta2.write(line)
	print "Rewrote fasta for %s." %spp
	fasta2.close()
