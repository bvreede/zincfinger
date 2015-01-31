import os

zffolder = "/home/barbara/data/zincfingers"

# replace pipe with underscore in the fasta files (do this before making dbs)
'''
for fasta in os.listdir(zffolder):
	f = open("%s/%s" %(zffolder,fasta))
	f2 = open("%s/%s2" %(zffolder,fasta), "w")
	for line in f:
		line = line.replace('|','_')
		f2.write(line)
	f2.close()
	os.remove("%s/%s" %(zffolder,fasta))
	os.rename("%s/%s2" %(zffolder,fasta), "%s/%s" %(zffolder,fasta))
'''


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
#for spp in ['dpul','dmel','tcas','turt','smar','isca']:
#	prepdb(spp)
#	prepblastdb(spp)

'''
# blast each file of dmel against the five other species, and save the name of its top hit
# output: csv file with dmel in column 1, spp in column 2; named: dmel-spp
for spp in ['dpul','tcas','turt','smar','isca']:
	outdb = open("%s/dmel-%s.csv" %(zffolder,spp), "w")
	for f in os.listdir("%s/%s" %(zffolder,spp)):
		if f[0] == '.':
			print "skipping hidden file ", f
			continue
		blastcommand = "blastp -query %s/%s/%s -db %s/%szfs.fa -out %s/tempout.txt" %(zffolder,spp,f,zffolder,spp,zffolder)
		os.system(blastcommand)
		blastres = open("%s/tempout.txt" %zffolder)
		for line in blastres:
			if line[0] == '>':
				line = line.strip()
				outdb.write("%s,%s\n" %(f,line[1:]))
				break
'''			

# blast each file from the five species against dmel and save the name of its top hit
# output: csv file with spp in column 1, dmel in column 2; named: spp-dmel
for spp in ['dpul','tcas','turt','smar','isca']:
	outdb = open("%s/%s-dmel.csv" %(zffolder,spp), "w")
	for f in os.listdir("%s/dmel" %(zffolder)):
		if f[0] == '.':
			print "skipping hidden file ", f
			continue
		blastcommand = "blastp -query %s/dmel/%s -db %s/dmelzfs.fa -out %s/tempout.txt" %(zffolder,f,zffolder,zffolder)
		os.system(blastcommand)
		blastres = open("%s/tempout.txt" %zffolder)
		for line in blastres:
			if line[0] == '>':
				line = line.strip()
				outdb.write("%s,%s\n" %(line[1:],f))
				break


