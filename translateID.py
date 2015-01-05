import sys,csv

to_translate = open(sys.argv[1])
IDlibrary = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/databases/150105-SM00355-dmel_IDs.csv"))
transl_name = sys.argv[1][:-4] + '_translID.txt'
translated = open(transl_name, "w")

IDlist = []
for l in IDlibrary:
	IDlist.append(l)

TRlist = []
symbols = []
for prot in to_translate:
	prot = prot.strip()
	for ID in IDlist:
		if ID[2] == prot:
			if ID[3] != '':
				TRlist.append(ID[3])
			elif ID[4] != '':
				TRlist.append(ID[4])
			elif ID[5] != '':
				TRlist.append(ID[5])
			else:
				symbols.append(ID[1])

TRlist = list(set(TRlist))
symbols = list(set(symbols))

for t in TRlist:
	translated.write("%s\n" %t)

if len(symbols) > 0:
	translated.write("\nSYMBOLS:\n")
	for s in symbols:
		translated.write("%s\n" %s)
translated.close()

	
