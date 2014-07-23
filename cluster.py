import scipy, pylab, csv
import scipy.cluster.hierarchy as sch
from collections import Counter

dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"

#generate array of csv file
# The motifs file path
GOfile = csv.reader(open("%s/results/motifhits_dmel-count.csv" %(dbfolder))) 
# the Go terms path
#GOfile = csv.reader(open("%s/databases/140720-SM00355-dmel-GO-G.csv" %(dbfolder)))

GOmatrix = []
genes = []
feature_labelsU = []
for line in GOfile:
        if line[0] == "Gene_stable_ID":
            # gather gene labels, skipping the first 3 columns, and 
            # ignoring the last column (which is '' due to the trailing comma"
            feature_labelsT = line[3:] #for labels of inverse/transposed 
            feature_labelsT = feature_labelsT[:-1]
            continue
	feature_labelsU.append(line[1]) #for labels of untransposed matrix: gene names
        GOline = []

        for i in range(3,len(line)-1): #line ends in comma, so this removes the last empty item, and removes the first line
            GOline.append(float(line[i]))
        GOmatrix.append(GOline)

        # add the gene identifier to genes labels
        genes.append(line[0])


D = scipy.array(GOmatrix)
Y = sch.linkage(D, method='weighted')
cutoff = 25#1.1547
#Z = sch.dendrogram(Y, labels=feature_labelsU)#, orientation='right')#labels=feature_labels, orientation='right')

#T=sch.fcluster(Y, cutoff,'distance')
#print T
#Z=sch.dendrogram(Y, color_threshold=cutoff)

L = sch.fclusterdata(D, cutoff, criterion='maxclust', method='weighted') #criterion='inconsistent', metric='euclidean', depth=2, , R=None) #
S = set(L)
Llist = list(L)
print len(S)

print Counter(Llist)
'''
for i in range(len(S)):
	if Llist.count(i+1) > 10:
		print i+1, Llist.count(i+1)

'''
