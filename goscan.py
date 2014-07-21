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

def gocheck(query):
	for i in golist:
		if query == i:
			outdb.write("1,")
		else:
			outdb.write(",")
#WORKING ON THIS! Need to make sure that it scans once per protein, right now it scans for every GO entry.

# goes through the database that needs to be checked (per organism)
# Make sure the database is curated and GO terms are in the right columns!
def dbscan(org):
	protdb = open("/home/barbara/Dropbox/zinc_finger_data/data/databases/140720-SM00355-%s2.csv" %org)
	outdb = open("/home/barbara/Dropbox/zinc_finger_data/data/databases/140720-SM00355-%s-GO.csv" %org, "w")
	testID = "Protein_stable_ID"
	for line in protdb:
		geneID = line.split(',')[0:3]
		if geneID[2] == testID:
			continue
		outdb.write("%s,%s,%s," %(geneID[0],geneID[1],geneID[2]))
		testID = geneID[2] # to test whether the proteinID has already been used

for org in ['smar','dmel','tcas','isca','dpul']: #all species used in this project
	dbscan(org)
