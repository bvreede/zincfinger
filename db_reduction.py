dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"
inputfile = open("%s/results/motifhits_dmel-count.csv" %dbfolder)
outputfile = open("%s/databases/140720-SM00355-dmel-motifhits-G.csv" %dbfolder, "w")

'''
import the file to be reduced into local memory
'''
M = []
for line in inputfile:
	l = line.strip().split(',')
	if l[0] == "Gene_stable_ID":
		outputfile.write(line)
		continue
	M.append(l)
M.sort() # sort the genes so that similar genes are grouped together

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
	for m in range(3,b-1): #from first to last column of results (which is now in rows); -1 because last trailing comma causes empty field
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
	if line[0] == lcheck:
		gm.append(line)
		continue
	lcollect = ""
	if gm != []:
		mtl = matrix_to_list(gm)
		for q in mtl:
			outputfile.write("%s," %q)
		outputfile.write("\n")
	gm = []
	gm.append(line)
	for i in line:
		lcollect += i
		lcollect += ','
	outputfile.write("%s\n" %lcollect[:-1])
	lcheck = line[0]
