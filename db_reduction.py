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

M.sort()

lcheck = ""
for line in M:
	if line[0] == lcheck:
		continue
	lcollect = ""
	for i in line:
		lcollect += i
		lcollect += ','
	outputfile.write("%s\n" %lcollect[:-1])
	lcheck = line[0]
