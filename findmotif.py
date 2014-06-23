from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.motifs import Instances
#m=motifs.Motif(alphabet=IUPAC.unambiguous_dna)

memeout = open("/home/barbara/Dropbox/zinc_finger_data/meme.txt")
fastadb = open("/home/barbara/Dropbox/zinc_finger_data/databases/zfonly-dmel-aaseq.fa")

### DEFINE MOTIFS ###
# 2-8-3
# 2-12-3
# 2-12-4
# 2-12-5
# 4-12-3
# 4-12-4
# 4-15-3
# C2HC
# P-DLS

motifsM = list(motifs.parse(memeout, "MEME"))
mainmotif = motifsM[0].consensus#.instances[1].pvalue
allmotifs = motifsM[0].instances

### SEARCH SEQUENCES ###
# for each motif: go through fasta file (or go through fasta file and search for multiple?)
# at each iteration of this motif: note start and end site on .csv file 
# also note: name, spp, sequence, motif type, group
# iteration of motif (order in the protein; for alignment purposes?)

#read fasta file and put in dictionary
fastadict = {}
sequence = ""
header = ""
for line in fastadb:
	if line[0] == ">":
		if header != "":
			fastadict[header] = sequence
		header = line[1:].strip()
		sequence = ""
	else:
		sequence += line.strip()
fastadict[header] = sequence

for key in fastadict:
	test_seq = Seq(fastadict[key],mainmotif.alphabet)
	for pos,seq in Instances(allmotifs).search(test_seq):
		print pos,seq#.tostring()

#Instances.search(Instances,test_seq)

#for pos,seq in mainmotif.search.Instances(test_seq):
#	print pos,seq.tostring()


### DRAW IMAGE ###
# for each group:
# draw line for length of sequence
# draw coloured box for each motif
# 
