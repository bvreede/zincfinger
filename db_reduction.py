dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
inputfile = open("%s/results/motifhits_dmel-count.csv" %dbfolder)
outputfile = open("%s/databases/140720-SM00355-dmel-motifhits-G.csv" %dbfolder, "w")

'''
import the file into local memory
'''
M = []
for line in inputfile:
	l = line.strip().split(',')
	if l[0] == "Gene_stable_ID":
		head = l
		continue
	M.append(l)
M.sort() # sort the genes so that all proteins from the same gene are grouped together
motiflist = ['2_8_3','2_12_3','2_12_4','2_12_5','4_12_3','4_12_4','4_15_3']

'''
write first line of the outputfile: header
'''
for n in head[:3]:
	outputfile.write("%s," %n)
for m in motiflist:
	outputfile.write("%s," %m)
outputfile.write("\n")

'''
takes line that has scores for motifs both inverse and forward,
and collapses that to find scores for the motifs in total
'''
def collapseinv(head,line):
	cdict = {}
	cl = line[:3]
	for i in range(3,len(head)-1):
		cdict[head[i]] = line[i]
	for m in motiflist:
		m_inv = m + "_inv"
		total = int(cdict[m]) #+ int(cdict[m_inv])
		#print cdict[m], cdict[m_inv], total
		cl.append(total)
	return cl

'''
takes gene-specific matrix with all observations
related to that matrix, and reduce them to a single line.
Also add up the _inv with the fwd motifs.
'''
def matrix_to_list(gm):
	a= len(gm)
	b= len(gm[0])
	mtl = gm[0][0:3] # writes the first three terms (gene ID, name, protID) to the output list
	gmt = zip(*gm) # transposes the matrix
	for m in range(3,b): #from first to last column of results (which is now in rows)
		mtl.append(max(gmt[m]))
	return mtl

'''
go through local memory to generate individual matrices
of all genes, and reduce them to single lines.
Write single lines to final document.
'''
lcheck = ""
gm = []
for line in M:
	lc = collapseinv(head,line)
	if lc[0] == lcheck:
		gm.append(lc)
		continue
	lcollect = ""
	if gm != []:
		mtl = matrix_to_list(gm)
		for q in mtl:
			outputfile.write("%s," %q)
		outputfile.write("\n")
	gm = []
	gm.append(lc)
	lcheck = lc[0]
