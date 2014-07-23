'''
Marc-specific specifications
'''
#import matplotlib
#matplotlib.use('Agg')

import scipy, pylab, csv
import scipy.cluster.hierarchy as sch

'''
Defining which files to use: use comment-out to switch between visualizing GO-term clustering
and motif clustering.
'''
dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"
# The motifs file path
GOfile = csv.reader(open("%s/results/motifhits_dmel-count.csv" %(dbfolder))) 
# the Go terms path
#GOfile = csv.reader(open("%s/databases/140720-SM00355-dmel-GO-G.csv" %(dbfolder)))


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
for line in GOfile:
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
Create the actual dendrogram
'''
fig = pylab.figure(figsize=(10, 20)) # create the figure, sizes are in inches
axdendro = fig.add_axes([0.09,0.1,0.2,0.8]) # add axes as: left bottom width height
Y = sch.linkage(GOmatrix, method='weighted') # calculate the clustering; methods are: weighted, single, average, complete, centroid, median, ward
cutoff = 4
Z = sch.dendrogram(Y, color_threshold=cutoff, orientation='right') #create the dendrogram. Optional: add 'labels=termlabels,' (but especially with larger matrices it makes sense to add them later!)

# get correct row label orders
index = Z['leaves']
GOmatrix = GOmatrix[index,:]
#GOmatrix = GOmatrix[:,index] # this makes it crash, but keeping it because it was in the original download...

axdendro.set_xticks([])
#axdendro.set_yticks([])
#axdendro.set_yticks(range(len(genenames)))
axdendro.set_yticklabels(genenames, fontsize=2)

'''
Create the heatmap
'''
x_start = 0.45
axmatrix = fig.add_axes([x_start, 0.1,1-x_start-0.15,0.8])
im = axmatrix.matshow(GOmatrix, aspect='auto', origin='right', cmap='PuRd') #cmap = color pattern. Play with this :)

axmatrix.set_yticks([]) # hides all t ticks
axmatrix.set_xticks([])
axmatrix.set_xticks(range(len(termlabels))) 
axmatrix.set_xticklabels(termlabels, rotation=90, fontsize=10) # writes the labels on the heatmap. Make sure to pick the right labels when the matrix is transposed!

axcolor = fig.add_axes([0.91,0.1,0.02,0.8]) # Plots the colorbar
pylab.colorbar(im, cax=axcolor)

'''
That's it! Save the figure...
'''
fig.savefig("dendrogram.png", dpi=1200)

