import csv

current = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/nameclusters6.csv"))
newgroups = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/averageclusters1.csv"))
out = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/nameclusters7A.csv", "w")


prot2clus = {}
for line in current:
	prot2clus[line[2]] = line[3]

for line in newgroups:
	out.write("%s,%s,%s,%s,%s," %(line[0],line[1],line[2],line[3],line[4]))
	if line[2] in prot2clus:
		out.write(prot2clus[line[2]])
	out.write("\n")

out.close()
