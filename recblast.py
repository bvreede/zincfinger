#!/usr/bin/python

'''
This script can be used, either entirely or in parts, to
(1) generate blast databases and query files from any fasta file;
(Fasta headers should look like: >[geneID]|[genename]|[protID])
(2) blast sequences from one fasta file against the other (and back). 
NB!!! Any fasta file that is the source for ortholog
labeling needs to have gene names specified in all genes.
(3) look at the resulting blast hits, extract the reciprocal best
blast hits, and rewrite the initial fasta inputfile to include the
ortholog.
(4) write a csv file which has, per Drosophila (or other source species)
gene, the orthologs of other species in columns.
Author: Barbara Vreede
Date: 4 February 2015
Contact: b.vreede@gmail.com
'''

import os,csv,sys

#CUSTOMIZE: zffolder is for the databases and individual query files; seqfolder is for the original fasta files
zffolder = "/home/barbara/data/zincfingers"
seqfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences"

#CUSTOMIZE: the names given to databases, and the names the fasta files currently have. This requires
#all fasta files used for the reciprocal blast to have the same name structure (+ individualized specifier);
#eg: 150111_dmel.fasta, 150111_tcas.fasta, 150111_isca.fasta would have '150111_' as prefix and '.fasta' as suffix.
dbsuffix = "zfs.fa" #the name that will be given to databases as they are indexed for blast
dbcurprefix = "150111-SM00355-" #the current fasta file name (prefix)
dbcursuffix = "_seq.fasta" #the current fasta file name (suffix)
dbnewsuffix = "_seqRW.fasta" #the rewritten fasta file name (suffix)

#CUSTOMIZE: the individual species specifiers for each fasta file ('ilistA') and the comparative species ('ilistB').
ilistA = ['isca','smar','turt','tcas','dpul']
comp = ['dmel']

#CUSTOMIZE: what do you want to do?
newdbs = 0 #set to 1 if you want to make new databases for blasting and searching; else set to 0
redoblast = 0 #set to 1 if you want to re-do the blasting; set to 0 if you want to use existing result files
rewritefa = 0 #set to 1 if you want to rewrite your original fasta files to include the dmel/comp ortholog name; else set to 0
tableout = 0 #set to 1 if you want a csv file with orthologs; else set to 0

#END CUSTOMIZATION.
ilistC = ilistA + comp
comp = comp[0] #comp should be a string, but was briefly a list for concatenation


'''
PART 0:
VERIFY EXISTENCE OF OUTPUTFOLDER AND COPY FILES THERE.
'''
if not os.path.exists(zffolder):
	os.system("mkdir %s" %zffolder)
else:
	answer = raw_input("The folder %s already exists. Continuing may overwrite existing files. Do you want to continue? (y/n)" %zffolder)
	if answer == 'n':
		sys.exit("Quitting...")
	elif answer != 'y':
		sys.exit("Answer invalid. Quitting anyway.")

for filename in os.listdir(seqfolder):
	if newdbs == 0:
		break
	k = len(dbcurprefix)
	if filename[:k] == dbcurprefix:
		# make sure only to take those files from the database folder that follow the name pattern indicated in customization
		if dbcursuffix not in filename:
			continue
		spp = filename.replace(dbcursuffix,'')[k:]
		if spp not in ilistC:
			continue
		os.system("cp %s/%s %s/%s%s" %(seqfolder,filename,zffolder,spp,dbsuffix))


'''
PART 1:
PREPARE DATABASES AND QUERY FILES FROM FASTA FILES DOWNLOADED FROM BIOMART
'''
# replace pipes and other characters in the fasta files (do this before making dbs and sequences!)
for fasta in os.listdir(zffolder):
	# escape if you don't want to do this
	if newdbs == 0:
		break
	# check if the fasta file is one of the original transfers
	spp = fasta.replace(dbsuffix,'')
	if spp not in ilistC:
		continue
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

# prep databases for future reciprocal searches by individually saving all sequences in the fasta
def prepdb(spp):
	infasta = open("%s/%s%s" %(zffolder,spp,dbsuffix))
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
	os.system("makeblastdb -in %s/%s%s -dbtype prot" %(zffolder,spp,dbsuffix))

# only necessary to do this once per species
for spp in ilistC:
	if newdbs == 0: #escape if you don't want to make a new db
		break
	prepdb(spp)
	prepblastdb(spp)

'''
PART 2:
PERFORM BLASTS FROM ALL SPECIES IN 'ILISTA' TO THE SPECIES IN 'COMP' AND BACK
'''
# blast each file of dmel against the five other species, and save the name of its top hit
# output: csv file with dmel (or other comparison) in column 1, spp in column 2; named: dmel-spp
for spp in ilistA:
	blastdb =  "%s/%s%s" %(zffolder,spp,dbsuffix)
	if redoblast == 0:
		break
	if not os.path.exists(blastdb):
		sys.exit("No blast db found. Create the databases, then run the script again.")
	outdb = open("%s/%s-%s.csv" %(zffolder,comp,spp), "w")
	for f in os.listdir("%s/%s" %(zffolder,comp)):
		if f[0] == '.':
			print "skipping hidden file ", f
			continue
		blastcommand = "blastp -query %s/%s/%s -db %s -out %s/tempout.txt" %(zffolder,comp,f,blastdb,zffolder)
		os.system(blastcommand)
		blastres = open("%s/tempout.txt" %zffolder)
		for line in blastres:
			if line[0] == '>':
				line = line.strip()
				outdb.write("%s,%s\n" %(f,line[1:]))
				break
	print "Blast from %s to %s complete." %(comp,spp)
	os.remove("%s/tempout.txt" %zffolder)
	outdb.close()

# blast each file from the five species against dmel and save the name of its top hit
# output: csv file with dmel in column 1, spp in column 2; named: spp-dmel
for spp in ilistA:
	if redoblast == 0:
		break
	outdb = open("%s/%s-%s.csv" %(zffolder,spp,comp), "w")
	for f in os.listdir("%s/%s" %(zffolder,spp)):
		if f[0] == '.':
			print "skipping hidden file ", f
			continue
		blastcommand = "blastp -query %s/%s/%s -db %s/%s%s -out %s/tempout.txt" %(zffolder,spp,f,zffolder,comp,dbsuffix,zffolder)
		os.system(blastcommand)
		blastres = open("%s/tempout.txt" %zffolder)
		for line in blastres:
			if line[0] == '>':
				line = line.strip()
				outdb.write("%s,%s\n" %(line[1:],f))
				break
	print "Blast from %s to %s complete." %(spp,comp)
	os.remove("%s/tempout.txt" %zffolder)
	outdb.close()

'''
PART 3:
CHECK WHICH RECIPROCAL HITS EXIST AND REWRITE THE FASTA FILES
'''
# make a reciprocal check of each pair of csv files
if redoblast == 1:
	gndict = open("%s/genenamedict.csv" %(seqfolder), "w")
	hitcollect = [] #list of lists containing hits to dmel genes per species
for spp in ilistA:
	if redoblast == 0:
		break
	hcspp = [spp] #the list containing all dmel hits for this species (and a header with the species name)
	todmel = csv.reader(open("%s/%s-%s.csv" %(zffolder,spp,comp)))
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
	fromdmel = csv.reader(open("%s/%s-%s.csv" %(zffolder,comp,spp)))
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
			hcspp.append(spg)
		else:
			hcspp.append("")
	hitcollect.append(hcspp)
	print "Made dictionary for %s." %spp
if redoblast == 1:
	gndict.close()

# read the gene translation dictionary and rewrite the species fasta files
gtransdx = {}
gndictcsv = "%s/genenamedict.csv" %(seqfolder)
if not os.path.exists(gndictcsv):
	sys.exit("No reciprocal best blast hit dictionary was found. Run the script again, and include the blast step.")
gtrans = csv.reader(open(gndictcsv))
for line in gtrans:
	gtransdx[line[0]] = line[1]
 
for spp in ilistA:
	if rewritefa == 0:
		break
	fasta = open('%s/%s%s%s' %(seqfolder,dbcurprefix,spp,dbcursuffix))
	fasta2 = open('%s/%s%s%s' %(seqfolder,dbcurprefix,spp,dbnewsuffix), "w")
	for line in fasta:
		if line[0] == '>':
			totalID = line.strip()
			totalID = totalID[1:]
			if totalID in gtransdx:
				gID,gna,pID = line.split('|')
				line = '%s|%s(%s-%s)|%s' %(gID,gna,comp[:2],gtransdx[totalID],pID)
		fasta2.write(line)
	print "Rewrote fasta for %s." %spp
	fasta2.close()
