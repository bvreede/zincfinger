allseqs = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/motifseq_allz.fasta")

lenseqdix = {}
seqdix = {}
seqdix2 = {}
alllens = []
switch = 0
for line in allseqs:
	if line[0] == '>':
		switch = 1
		head = line.strip()[1:]
	elif switch == 1:
		seq1 = line.strip()
		n = seq1.count('|')
		seq = seq1.replace('{','')
		seq = seq.replace('}','')
		seq = seq.replace('|','')
		seq = seq.replace('Z','')
		lenseq = len(seq) - n
		switch = 2
		alllens.append(lenseq)
	if switch == 2:
		lenseqdix[head] = lenseq
		seqdix[head] = seq1
		seqdix2[head] = seq

print max(alllens)

for key in lenseqdix:
	if lenseqdix[key] == max(alllens):
		print key
		print seqdix[key]
		print seqdix2[key]
