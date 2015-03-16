'''
the point of this is to calculate: of each known ortholog (per spp), how great are the odds that they are in the same cluster as the Dmel?
'''


import csv


dbfolder = "/home/barbara/Dropbox/shared_work/zinc-finger-data/data"
clusters = "results/clustering-string_allz-average.csv"
orthologs = "sequences/dmel-orthologs.csv"


def orthos(genename):
	'''
	Look up per gene what their Dmel ortholog is. NB this may be multiple, as it could be multiple isoforms.
	'''
	return isoforms #list of all isoforms
 

def clusters(genes): #genes is a list
	'''
	Look up the cluster they are assigned (both the Dmel and the Spp gene)
	'''
	return clusters #list of clusters, same order as genes



#Collect all genes that are orthologs per species

#If one of the isoforms has the same cluster hit as the Spp ortholog, score as 1



