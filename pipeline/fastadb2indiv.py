#original db
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences"
faname = "150111-SM00355-"
faappend = "_seq.fasta"
db = open("%s/%sall2%s" %(dbfolder,faname,faappend))

#read original db
#assemble species names
dbdix = {}
head, seq = "", ""
spp = []
for line in db:
	if line[0] == ">":
		if len(head) > 0:
			dbdix[head] = seq
			seq = ""
		head = line.strip()
		sp = head.split('|')[1][:4].lower() #head is eg. >ISCW007582|ISCA-(dm-Zn72D)|ISCW007582-PA so this selects ISCA and makes it lower case.
		if sp not in spp:
			spp.append(sp)
	else:
		seq += line.strip()
dbdix[head] = seq

#open result dbs for all possible species
for sp in spp:
	spdb = open("%s/%s%s2%s" %(dbfolder,faname,sp,faappend), "w")
	spdb.close()

#for each fasta entry in original db
#write it in the new species file (append)
for head in dbdix:
	sp = head.split('|')[1][:4].lower()
	spdb = open("%s/%s%s2%s" %(dbfolder,faname,sp,faappend), "a")
	spdb.write("%s\n%s\n\n" %(head,dbdix[head]))
	spdb.close()



