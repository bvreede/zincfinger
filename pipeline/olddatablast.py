import os,csv

indb = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/Additional_file_2"
groups = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/originalgroups.csv"))
outfinal = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/nameclusters1.csv", "w")

name2group = {}
for name in groups:
	name2group[name[0]] = name[1]


'''
Function to blast the protein to a database and return the first hit.
'''
def blast(filename,spp):
	if spp == 'dmel':
		blastdb = "/home/barbara/data/zincfingers/%szfs.fa" %spp
	else:
		blastdb = "/home/barbara/data/zincfingers/%szfsTRANS.fa" %spp
	query = "%s/%s" %(indb,filename)
	out = "%s/tempout.txt" %(indb)
	blast = "blastp -query %s -db %s -out %s" %(query,blastdb,out)
	os.system(blast)
	readout = open(out)
	result = 0
	for line in readout:
		if line[0] == '>':
			line = line.strip()
			hit  = line[2:]
			result = 1
			break
	os.remove(out)
	if result == 0:
		hit = "no hit"
	else:
		hit = hit.replace('-111-','(')
		hit = hit.replace('-222-',')')
		hit = hit.replace('dm:', 'dm-')
		a,b,c = hit.split('-333-')
		hit = b + '|' + c
	return hit


for fasta in os.listdir(indb):
	spp = fasta.split('_')[0].lower()
	hit = blast(fasta,spp)
	outfinal.write("%s,%s,%s\n" %(hit,name2group[fasta],fasta))

outfinal.close()
