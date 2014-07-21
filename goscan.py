'''
This script scans all proteins, scores the GO-terms, and saves the results
in a separate csv file, where each GO term has a column.
'''

#definition of terms and files:
gofile = "/home/barbara/Dropbox/zinc_finger_data/data/goterms_curated3.csv" #file that contains the GO-terms
pg = "G" #"P" = protein; "G" = gene
dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data/databases/" #folder for the in and output dbs
species = ['smar','dmel','tcas','isca','dpul'] # species abbreviations used in file names


# gets all the GO-term accession numbers in a list
goterms = open(gofile)
golist = []
for line in goterms:
	golist.append(line.split(',')[1])
golist.sort()
print len(golist)

# goes through GO list and scores whether the term occurs in the protein
def gocheck(querylist,outdb):
	if querylist[0] == "GO_term_accession":
		outdb.write("Gene_stable_ID,Gene_name,Protein_stable_ID,")
		for i in golist:
			outdb.write("%s," %i)
		outdb.write("\n")
	else:
		for i in golist:
			#outdb.write("%s," %(querylist.count(i)))
			if i in querylist:
				outdb.write("1,")
			else:
				outdb.write("0,")
		outdb.write("\n")
	

# goes through the database that needs to be checked (per organism)
# Make sure the database is curated and GO terms are in the right columns!
def dbscan(org):
	protdb = open("%s140720-SM00355-%s2.csv" %(dbfolder,org))
	outdb = open("%s140720-SM00355-%s-GO-%s.csv" %(dbfolder,org,pg), "w")
	if pg == "P":
		testID = "Protein_stable_ID"
		idno = 2
	elif pg == "G":
		testID = "Gene_stable_ID"
		idno = 0
	else:
		print "ERROR! Use P or G to indicate Protein or Gene info collection"
	querylist = [] # collects all GO-accessions from a single gene/protein
	for line in protdb:
		geneID = line.split(',')[0:3]
		if geneID[idno] == testID:
			querylist.append(line.split(',')[6]) # adds the GO-accession to the querylist
		else:
			gocheck(querylist,outdb)
			querylist = [] #empty querylist for each new gene/protein
			outdb.write("%s,%s,%s," %(geneID[0],geneID[1],geneID[2]))
			testID = geneID[idno] #replaces the testID with the gene/proteinID to ensure that all GO terms for this protein will be collected.
			querylist.append(line.split(',')[6]) # adds the GO-accession to the querylist
	gocheck(querylist,outdb) #for last protein
	
for org in species: #all species used in this project
	dbscan(org)
