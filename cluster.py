import scipy, pylab, csv
import scipy.cluster.hierarchy as sch
from collections import Counter

species = "dmel"
dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"
motiffile = "%s/databases/140720-SM00355-%s-motifhits-G.csv" %(dbfolder,species) #the file used for the clustering
infile = "%s/results/motifhits_%s.csv" %(dbfolder,species) #the file to apply the clustering to, and split into new files
inputfile = csv.reader(open(motiffile)) 
clustermeth = "weighted"
clusterorder = "%s/results/clustering_%s-%s.csv" %(dbfolder,species,clustermeth)
orderfile = open(clusterorder, "w")

'''
Read the input file and extract:
- matrix of observations (GOmatrix)
- labels of GOterms or motifs
- labels of genes (both numbers and names)
provides options to transpose the matrix
'''
observations = [] #will collect the actual matrix of observations
genenumbers = [] #will collect gene numbers (column 1)
genenames = [] #will collect gene names
for line in inputfile:
        if line[0] == "Gene_stable_ID": # gather GOterm/motif labels
            termlabels = line[3:] # skipping the first 3 columns
            termlabels = termlabels[:-1] # ignoring the last column (which is '' due to the trailing comma)
            continue
	genenumbers.append(line[0])
	genenames.append(line[1])
        GOline = [] #collects the observations per line
        for i in range(3,len(line)-1): #line ends in comma, so this removes the last empty item, and removes the first identifiers
            GOline.append(float(line[i]))
        observations.append(GOline)
GOmatrix = scipy.array(observations) # matrixify the observations
#GOmatrix = scipy.transpose(GOmatrix) # flips the matrix --- comment out if necessary!

'''
Calculate the clusters. Uses a cutoff value to define the number
of clustered categories that will be made.
'''
cutoff = 25 # determine cutoff: number of categories to be formed(1.1547)
L = sch.fclusterdata(GOmatrix, cutoff, criterion='maxclust', method=clustermeth) 
S = set(L) #turns the clustering into a set so as to remove duplicates
Llist = list(L) #turns the clustering into a list, so it may be indexed

cldict = {}
for i in range(len(genenumbers)):
	orderfile.write("%s,%s\n" %(genenumbers[i],Llist[i]))
	cldict[genenumbers[i]] = Llist[i]
orderfile.close()


print len(S) #returns the total number of categories (should be equal to cutoff)
print Counter(Llist) #counts instances per category

#for each category, go into GO and motif files, and save the genes seperately
for c in range(1,cutoff+1):
	outfile = infile[:-4] + "-cluster" + str(c) + '.csv'
	selection = open(outfile, "w") #resultsfile (file that is cluster-specific)
	to_select = open(infile) #file to read
	for line in to_select:
		l = line.strip().split(',')
		if l[0] == "Gene_stable_ID":
			selection.write(line)
		elif l[0] in cldict:
			if cldict[l[0]] == c:
				selection.write(line)
	selection.close()
	to_select.close()
