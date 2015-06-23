import config


infilebrev = infile.split('/')[-1].split('.')[0]
heatmaptxt = "%s/%s/%s_conservation-heatmap" %(config.mainfolder,config.evfolder,infilebrev)


#Make a dictionary for the species (abbreviation to full name)
species=open("%s/species.txt" %config.mainfolder)
spdict = {sp.strip().split()[0][0].lower()+sp.strip().split()[1][:3]: sp.strip() for sp in species}

#Make dictionary for 


def makeevolviewheatmap(doublematrix,name):
	'''
	Takes a list of lists (heatmap data) and uses motifs as labels
	and returns a text file that can be read by evolview.
	'''
	# Open the heatmap text file:
	evhm = open("%s-%s.txt" %(heatmaptxt,name), "w")
	evhm.write(" #heatmap\n !legendTitle\tFrequency of overlap\n !showLegends\t1\n !colorgradient\tfloralwhite,orange,red,purple,navy\
\n !colorgradientMarkLabel\t0,0.2,0.4,0.6,0.8,1\n # -- heatmap column labels --\n !showHeatMapColumnLabel\t1\n !heatmapColumnLabels\t")
	# heatmap requires motiflist in sequence, with commas between them
	motifscomma = ','.join(config.motiflist) 
	evhm.write("%s\n" %motifscomma)
	evhm.write(" # -- heatmap --\n !heatmap\tmargin=1,colwidth=18,roundedcorner=1\n # -- show data value\n !showdataValue\tshow=0,fontsize=12,fontitalic=0,textalign=start\n\n")
	for i,line in enumerate(doublematrix):
		evhm.write("%s\t" %config.motiflist[i])
		linew = ""
		for l in line:
			linew += str(l)
			linew += ','
		evhm.write("%s\n" %linew[:-1])
	evhm.close()





