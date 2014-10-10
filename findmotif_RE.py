'''
This script can be used to detect distinct C2H2 zinc finger motifs
in a protein sequence.
The output is (1) a csv file that scores the total number of domains
per type; (2) a fasta file with domains in order (+ a readme for the
fasta file); (3) a csv file that can be read by an image interpreter
(+ a readme for the file). 
'''
import re

### SPECIFY INFORMATION: DATA TO USE###
### input consists of: dbfolder/seqfolder/prefix-species_suffix
### output consists of: dbfolder/resfolder/
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
seqfolder = "sequences"
resfolder = "results"
prefix = "140720-SM00355-"
suffix = "_seq.fasta"
species = ["dmel","tcas","dpul","isca","smar"]
motiflist = ['2_8_3','2_12_3','2_12_4','2_12_5','4_12_3','4_12_4','4_15_3'] # use only numbers: C-H distances separated by _

### END OF CUSTOM INFO. DON'T FUCK WITH THE FOLLOWING ###

'''
Define regular expressions for all zf-domains (by the C-H distances), and save in a
dictionary. Each domain is only represented as a forward domain..
'''
motifdict = {}
for m in motiflist:
	cc,ch,hh = m.split('_')
	motif = 'C[A-Z]{%s}C[A-Z]{%s}H[A-Z]{%s}H' %(cc,ch,hh)
	remotif = re.compile(motif)
	motifdict[m] = remotif

'''
From an open fasta file, generate a dictionary containing the header as key, and the
sequence as value.
'''
def fastadicter(fastadb):
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
	return fastadict

'''
Turn the csv file into a matrix with single genes in the rows, and a motif count
in the columns.
Only takes into account the first protein per gene
'''
def matrix(sp):
	prev_output = open("%s/%s/motifhits_%s.csv" %(dbfolder,resfolder,sp)) # outputfile that saves motif locations and sequence length
	new_output = open("%s/%s/motifhits_%s-count.csv" %(dbfolder,resfolder,sp), "w") #outputfile that saves n motifs
	testID = "Gene_stable_ID"
	for line in prev_output:
		l = line.strip().split(',')
		if l[0] == testID:
			if l[0] == "Gene_stable_ID":
				lr = line.replace('Sequence_length,','')
				new_output.write(lr)
			continue
		else:
			new_output.write("%s,%s,%s," %(l[0],l[1],l[2]))
			for i in range(4,len(l)-1):
				if len(l[i]) != 0:
					new_output.write(str(len(l[i].split('|'))))
					new_output.write(",")
				else:
					new_output.write("0,")
			new_output.write("\n")

'''
Per species, read the fasta file, open an output file, and scan all sequences for the presence of
motifs.
'''
for sp in species:
	fastadb = open("%s/%s/%s-%s_%s" %(dbfolder,seqfolder,prefix,sp,suffix))
	outputdb = open("%s/%s/motifhits_%s.csv" %(dbfolder,resfolder,sp), "w")
	outputdb.write("Gene_stable_ID,Gene_name,Protein_stable_ID,Sequence_length,")
	for m in motifdict:
		outputdb.write("%s," %m)
	outputdb.write("\n")
	fastadict = fastadicter(fastadb)
	for key in fastadict:
		ids = key.split('|')
		seqlen = len(fastadict[key]) #length of the sequence
		outputdb.write("%s,%s,%s,%s," %(ids[0],ids[1],ids[2],seqlen)) #turn the header name into gene ID/name/prot ID
		for m in motifdict: #go through each motif and find all instances in the sequence
			domain = motifdict[m]
			positions = ''
			for i in domain.finditer(fastadict[key]):
				positions += str(i.start())
				positions += '|'
			outputdb.write('%s,' %(positions[:-1])) #remove final pipe from total positions
		outputdb.write("\n")
	outputdb.close()
	matrix(sp)
	


