import csv, sys

maindb = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters"
infile = "nameclusters4.csv"

#fasta files with sequences
fasta1 = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-" #prefix
fasta2 = "_seq-trans.fasta" #postfix
spp = ['turt','tcas','dmel','isca','smar','dpul']
allfasta = ["%s%s%s" %(fasta1,s,fasta2) for s in spp]



'''
Read the input from the cluster file (= gene names and their assigned cluster)
and save to local memory.
'''
infileo = csv.reader(open("%s/%s" %(maindb,infile)))

clusternames,clusterinfo,singlecluster = [],[],[]
curcluster = ''
#singlecluster will collect information per cluster
#clusternames will save all cluster names, in order of appearance
#clusterinfo is the total 3d matrix with clusters in D1 (0), genes in D2 (1), and ID,name,protID in D3 (2)
for line in infileo:
	if line[3] in clusternames:
		if curcluster != line[3]:
			sys.exit("Script is not robust to clusters that are not in order. Sort and try again!")
		else:
			singlecluster.append(line[0:3])
	else:
		clusternames.append(line[3])
		if len(singlecluster) >= 1:
			clusterinfo.append(singlecluster)
		singlecluster = []
		singlecluster.append(line[0:3])
		curcluster = line[3]
clusternames.append(line[3])
if len(singlecluster) >= 1:
	clusterinfo.append(singlecluster)

'''
Save all fasta files to memory so sequences can be recalled per cluster.
'''
sequences = {}
for f in allfasta:
	fo = open(f)


'''

input file: nameclusters4

output:
- per cluster: fasta file with sequences (send to Ariel and have him make sequence trees)
- per cluster: motif scores in separate file (and run the viz script on this)
- Drosophila genes only; gorilla input list
- Drosophila genes only; with gorilla output: make visualization (heatmap)

'''

for i,c in enumerate(clusterinfo):
	clusterid = clusternames[i]
	
