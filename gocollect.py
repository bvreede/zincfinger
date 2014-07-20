'''
Script takes the curated csv files containing GO terms (downloaded from BioMart; GO domain, accession, and name)
and collects all GO terms in one single file, removing the duplicates.
It assumes that there is a single domain/name combination for each accession number.
'''


goterms = open("/home/barbara/Dropbox/zinc_finger_data/data/goterms.csv", "w")
godic = {}

def gofind(org):
	db = open("/home/barbara/Dropbox/zinc_finger_data/data/databases/140720-SM00355-%s2.csv" %org) #curated databases with GO-terms in columns 5-7.
	for line in db:
		godom,goacc,goname = line.strip().replace('"','').split(",")[5:8] #'"' replace is due to bug in the database! Pay attention to this in the results
		godic[goacc] = godom + ":" + goname

for org in ['smar','dmel','tcas','isca','dpul']: #all species used in this project
	gofind(org)

for key in godic:
	if key == "" or key == "GO_term_accession":
		continue
	goterms.write("%s,%s\n" %(key,godic[key])) #output: 'Go Accession', 'Go Domain:Go name'

