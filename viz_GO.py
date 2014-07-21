# marc needs to do this for shit to work
#import matplotlib
#matplotlib.use('Agg')

import scipy
import pylab
import csv
import scipy.cluster.hierarchy as sch

dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"

#generate array of csv file
# The motifs file path
GOfile = csv.reader(open("%s/results/motifhits_dmel-count.csv" %(dbfolder))) 
# the Go terms path
#GOfile = csv.reader(open("%s/databases/140720-SM00355-dmel-GO-G.csv" %(dbfolder)))

GOmatrix = []
genes = []
for line in GOfile:
        if line[0] == "Gene_stable_ID":
            # gather gene labels, skipping the first 3 columns, and 
            # ignoring the last column (which is '' due to the trailing comma"
            feature_labels = line[3:]
            feature_labels = feature_labels[:-1]
            continue
        GOline = []

        for i in range(3,len(line)-1): #line ends in comma, so this removes the last empty item, and removes the first line
            GOline.append(float(line[i]))
        GOmatrix.append(GOline)

        # add the gene identifier to genes labels
        genes.append(line[0])


D = scipy.array(GOmatrix)

# flippin' da matrix!!!!
#D = scipy.transpose(D)

# Compute and plot dendrogram.

# create a wide figure, sizes are in "inches" WTF???
fig = pylab.figure(figsize=(25, 10))

# add axes as: left bottom width height
axdendro = fig.add_axes([0.09,0.1,0.2,0.8])

Y = sch.linkage(D, method='weighted')
Z = sch.dendrogram(Y, orientation='left') #labels=feature_labels, orientation='right')

# get correct row label orders
index = Z['leaves']
D = D[index,:]
#D = D[:,index]

axdendro.set_xticks([])

# Plot distance matrix.
x_start = 0.45
axmatrix = fig.add_axes([x_start, 0.1,1-x_start-0.15,0.8])

im = axmatrix.matshow(D, aspect='auto', origin='right', cmap='PuRd')

# hide all t ticks
#axmatrix.set_yticks([])
#axmatrix.set_xticks([])

# draw super tiny gene names here
#axmatrix.set_xticks(range(len(genes)))
#axmatrix.set_xticklabels(genes, rotation=270, fontsize=0.2)

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
pylab.colorbar(im, cax=axcolor)

# Display and save figure.
# fig.show()
fig.savefig("dendrogram.png", dpi=1200)

