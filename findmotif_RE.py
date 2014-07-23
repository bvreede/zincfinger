import re

dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data/"
species = ["dmel","tcas","dpul","isca","smar"]
motiflist = ['2_8_3','2_12_3','2_12_4','2_12_5','4_12_3','4_12_4','4_15_3']#,'C2HC','PDLS'] --- use only numbers C-H distances separated by _

'''
Define regular expressions for all zf-domains (by the C-H distances), and save in a
dictionary. Each domain is represented both as a forward and an inverse.
'''
motifdict = {}
for m in motiflist:
	cc,ch,hh = m.split('_')
	motif = 'C[A-Z]{%s}C[A-Z]{%s}H[A-Z]{%s}H' %(cc,ch,hh)
	motif_inv = 'H[A-Z]{%s}H[A-Z]{%s}C[A-Z]{%s}C' %(hh,ch,cc)
	remotif = re.compile(motif)
	remotif_inv = re.compile(motif_inv)
	motifdict[m] = remotif	
	minv = m + "_inv"
	motifdict[minv] = remotif_inv

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
	prev_output = open("%sresults/motifhits_%s.csv" %(dbfolder,sp)) # outputfile that saves motif locations and sequence length
	new_output = open("%sresults/motifhits_%s-count.csv" %(dbfolder,sp), "w") #outputfile that saves n motifs
	testID = "Gene_stable_ID"
	for line in prev_output:
		l = line.strip().split(',')
		if l[0] == testID:
			if l[0] == "Gene_stable_ID":
				new_output.write(line)
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
	fastadb = open("%ssequences/140720-SM00355-%s_seq.fasta" %(dbfolder,sp))
	outputdb = open("%sresults/motifhits_%s.csv" %(dbfolder,sp), "w")
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
	


