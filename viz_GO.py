import scipy
import pylab
import csv
import scipy.cluster.hierarchy as sch


#generate array of csv file
GOfile = csv.reader(open("/home/barbara/Dropbox/zinc_finger_data/data/results/motifhits_dmel-count.csv")) 
#GOfile = csv.reader(open("/home/barbara/Dropbox/zinc_finger_data/data/databases/140720-SM00355-dmel-GO-G.csv"))
GOmatrix = []
for line in GOfile:
	if line[0] == "Gene_stable_ID":
		continue
	GOline = []
	for i in range(3,len(line)-1): #line ends in comma, so this removes the last empty item, and removes the first line
		GOline.append(float(line[i]))
	GOmatrix.append(GOline)
D = scipy.array(GOmatrix)

# Compute and plot dendrogram.
fig = pylab.figure()
axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
Y = sch.linkage(D, method='weighted')
Z = sch.dendrogram(Y, orientation='right')
axdendro.set_xticks([])
axdendro.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
index = Z['leaves']
D = D[index,:]
D = D[:,index]
im = axmatrix.matshow(D, aspect='auto', origin='lower')
axmatrix.set_xticks([])
axmatrix.set_yticks([])

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
pylab.colorbar(im, cax=axcolor)

# Display and save figure.
fig.show()
fig.savefig("dendrogram.png")
