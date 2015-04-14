import csv,re

# read fasta file into memory
fasta = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-allz_seq.fasta")
nametrans = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/identified_sequences_ariel.csv"))
newfasta = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/150111-SM00355-all2_seq.fasta", "w")


names = {}
for line in nametrans:
	protein = line[0].split('|')[1]
	names[protein] = line[1]



# turn this on if you want to search for all possible combinations of the above: (and don't forget to turn the custom list off!)
#motiflist = ['%s_%s_%s' %(m,n,o) for m in moCC for n in moCH for o in moHH] 

# turn this on for custom motif list
motiflist = ['2_7_4','2_8_3','2_9_3','2_10_5','2_11_3','2_11_4','2_12_2','2_12_3','2_12_4','2_12_5','2_12_6','2_13_3',
'2_13_4','2_14_3','2_14_4','2_15_4','3_8_3','4_12_3','4_12_4','4_15_3']

'''
Define regular expressions for all zf-domains (by the C-H distances), and save in a
dictionary. Each domain is only represented as a forward domain.
'''
motifdict = {} # dictionary of motifs and their regular expression
plink = 4
alink = 4
for m in motiflist:
	cc,ch,hh = m.split('_')
	motif = '[A-Z]{%s}C[A-Z]{%s}C[A-Z]{%s}H[A-Z]{%s}H[A-Z]{%s}' %(plink,cc,ch,hh,alink) # construct the regular expression
	remotif = re.compile(motif)
	motifdict[m] = remotif

def zfcheck(sequence):
	result = 0
	for m in motiflist:
		#search if motif appears
		if re.search(motifdict[m],sequence) != None:
			result = 1
			break
	return result	

seqdict,transdict = {},{}
# for each header:
header,sequence = '',''
for line in fasta:
	if line[0] == '>':
		#save previous header/sequence combo
		if header != '':
			if zfcheck(sequence) == 1:
				seqdict[header] = sequence
		header = line.strip()
		header = header.replace(':','-') #fix for evolview
		nameli = header[1:].split('|')
		if nameli[2] in names:
			testnew = '(' + names[nameli[2]][5:] + ')'
			nameli[1] = nameli[1] + testnew
		#spacer only necessary when there is a name in the genename field
		if len(nameli[1]) > 0:
			spacer = '-'
		else:
			spacer = ''
		#add species identifier to the name
		if nameli[2][:2] == 'SM':
			nameli[1] = 'SMAR' + spacer + nameli[1]
		elif nameli[2][:2] == 'TC':
			nameli[1] = 'TCAS' + spacer + nameli[1]
		elif nameli[2][:2] == 'IS':
			nameli[1] = 'ISCA' + spacer + nameli[1]
		elif nameli[2][:2] == 'EF':
			nameli[1] = 'DPUL' + spacer + nameli[1]
		elif nameli[2][:2] == 'te':
			nameli[1] = 'TURT' + spacer + nameli[1]
		elif nameli[2][:2] == 'FB':
			nameli[1] = 'DMEL' + spacer + nameli[1]
		#save in dictionary
		transdict[header] = ">%s|%s|%s" %(nameli[0],nameli[1],nameli[2])
		#empty sequence
		sequence = ''
	else:
		sequence += line.strip()
if zfcheck(sequence) == 1:
	seqdict[header] = sequence

# what is left: write to a new fasta file
for key in seqdict:
	header = transdict[key]
	sequence = seqdict[key]
	newfasta.write("%s\n%s\n\n" %(header,sequence))

newfasta.close()
