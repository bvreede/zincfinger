'''
This script scans all proteins, scores the GO-terms, and saves the results
in a separate csv file, where each GO term has a column.
'''

# gets all the GO-term accession numbers in a list
goterms = open("/home/barbara/Dropbox/zinc_finger_data/data/goterms.csv")
golist = []
for line in goterms:
	golist.append(line.split(',')[0])
golist.sort()

def gocheck(querylist,outdb):
	for i in golist:
		if querylist[0] == "Go_term_accession":
				continue
		else:
			outdb.write("%s," %(querylist.count(i)))
		outdb.write("\n")
	

# goes through the database that needs to be checked (per organism)
# Make sure the database is curated and GO terms are in the right columns!
def dbscan(org):
	protdb = open("/home/barbara/Dropbox/zinc_finger_data/data/databases/140720-SM00355-%s2.csv" %org)
	outdb = open("/home/barbara/Dropbox/zinc_finger_data/data/databases/140720-SM00355-%s-GO.csv" %org, "w")
	testID = "Protein_stable_ID" #the first line in the csv doc
	querylist = [] # collects all GO-accessions from a single protein
	for line in protdb:
		geneID = line.split(',')[0:3]
		if geneID[2] == testID:
			querylist.append(line.split(',')[6]) # adds the GO-accession to the querylist
		else:
			gocheck(querylist,outdb)
			querylist = [] #empty querylist for each new protein
			outdb.write("%s,%s,%s," %(geneID[0],geneID[1],geneID[2]))
			testID = geneID[2] #replaces the testID with the proteinID to ensure that all GO terms for this protein will be collected.
			querylist.append(line.split(',')[6]) # adds the GO-accession to the querylist
	gocheck(querylist,outdb) #for last protein
	
for org in ['smar','dmel','tcas','isca','dpul']: #all species used in this project
	dbscan(org)
