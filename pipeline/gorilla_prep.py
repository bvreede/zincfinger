import os,csv

infolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/1.1547-clusters"

# prep library
IDlibrary = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/databases/150105-SM00355-dmel_IDs.csv"))
IDlist = []
for l in IDlibrary:
	IDlist.append(l)

for filename in os.listdir(infolder):
	# open file
	f = csv.reader(open("%s/%s" %(infolder,filename)))
	# get column 3 translated into list
	li = []
	for p in f:
		li.append(p[2])
	# translate IDs for GOrilla use
	TRlist,symbols = [],[]
	for i in li:
		for ID in IDlist:
			if ID[2] == i:
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
	o = open("%s/%s_transID.txt" %(infolder,filename[:-4]), "w")
	for t in TRlist:
		o.write("%s\n" %t)
	if len(symbols) > 0:
		o.write("\nSYMBOLS:\n")
		for s in symbols:
			o.write("%s\n" %s)
	o.close()

