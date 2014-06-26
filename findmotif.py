from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.motifs import Instances
<<<<<<< HEAD
=======
import matplotlib
>>>>>>> 8c9b40ec43cc16009a8b928d26a95d3b6a5f678e

fastadb = open("/home/barbara/Dropbox/zinc_finger_data/databases/zfonly-dmel-aaseq.fa")
outputdb = open("motifhits.csv", "w")

### DEFINE MOTIFS ###
# 2-8-3
meme_2_8_3 = open("/home/barbara/Dropbox/zinc_finger_data/meme/2_8_3.txt")
# 2-12-3
meme_2_12_3 = open("/home/barbara/Dropbox/zinc_finger_data/2_12_3.txt")
# 2-12-4
meme_2_12_4 = open("/home/barbara/Dropbox/zinc_finger_data/2_12_4.txt")
# 2-12-5
meme_2_12_5 = open("/home/barbara/Dropbox/zinc_finger_data/2_12_5txt")
# 4-12-3
meme_4_12_3 = open("/home/barbara/Dropbox/zinc_finger_data/4_12_3.txt")
# 4-12-4
meme_4_12_4 = open("/home/barbara/Dropbox/zinc_finger_data/4_12_4.txt")
# 4-15-3
meme_4_15_3 = open("/home/barbara/Dropbox/zinc_finger_data/4_15_3.txt")
# C2HC
meme_C2HC = open("/home/barbara/Dropbox/zinc_finger_data/C2HC.txt")
# P-DLS
<<<<<<< HEAD
=======
meme_PDLS = open("/home/barbara/Dropbox/zinc_finger_data/PDLS.txt")

>>>>>>> 8c9b40ec43cc16009a8b928d26a95d3b6a5f678e
motifsM = list(motifs.parse(memeout, "MEME"))
motifname = "2-12-3"
mainmotif = motifsM[0].consensus
allmotifs = motifsM[0].instances

### SEARCH SEQUENCES ###
# for each motif: go through fasta file (or go through fasta file and search for multiple?)
# at each iteration of this motif: note start site on .csv file 
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
	hits = [key]
	seqs = [fastadict[key]]
	outputdb.write("%s,%s" %(key,motifname))
	for pos,seq in Instances(allmotifs).search(test_seq):
		hits.append(pos)
		seqs.append(seq)
<<<<<<< HEAD
		#print pos,seq#.tostring()
=======
		print pos,seq#.tostring()
>>>>>>> 8c9b40ec43cc16009a8b928d26a95d3b6a5f678e
	outputdb.write(",%s" %(len(hits)-1))
	if len(hits) > 1:
		print hits
		for i in range(len(hits)-1):
			pos = hits[i+1]
			seq = seqs[i+1]
			outputdb.write(",%s,%s" %(pos,seq))
	outputdb.write("\n")
<<<<<<< HEAD
=======


>>>>>>> 8c9b40ec43cc16009a8b928d26a95d3b6a5f678e



### DRAW IMAGE ###
# for each group:
# draw line for length of sequence
# draw coloured box for each motif
# 
