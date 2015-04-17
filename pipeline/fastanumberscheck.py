infile = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-smar_seqTR.fasta")


counthead = 0
genes = []
for line in infile:
	if line[0] == '>':
		counthead += 1
		gene = line.split('|')[0]
		genes.append(gene)

genes = set(genes)

print "total isoforms: ", counthead
print "total genes: ", len(genes)
