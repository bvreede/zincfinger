import config,csv,sys

if len(sys.argv) <= 1:
	sys.exit("USAGE: python evolview_motifcount.py path/to/inputfile (input is the summary file where total motifcounts are stored after 'findmotif.py'. Usually in 'results' folder. Ensure you know when it was created and check it before use!)")


#INPUT FILES
infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('.')[0]
infile_nonamb = infile.split('.')[0] + "-nonambg.csv"

#OUTPUT FILES
heatmaptxt = "%s/%s/%s_heatmap" %(config.mainfolder,config.evfolder,infilebrev)
bargraphtxt = "%s/%s/%s_bargraph" %(config.mainfolder,config.evfolder,infilebrev)

def makeevolviewheatmap(resultsdict,name):
	'''
	Takes a dictionary with hit results and uses corresponding species as labels.
	Returns a text file that can be read by evolview.
	'''
	# Open the heatmap text file:
	evhm = open("%s-%s.txt" %(heatmaptxt,name), "w")

	# Write the evolview headers
	evhm.write(" #heatmap\n !legendTitle\tFrequency of occurrence\n !showLegends\t1\n !colorgradient\tfloralwhite,orange,red,purple,navy\
\n !colorgradientMarkLabel\t0,0.2,0.4,0.6,0.8,1\n # -- heatmap column labels --\n !showHeatMapColumnLabel\t1\n !heatmapColumnLabels\t")
	motifscomma = ','.join(motifs) # heatmap requires motiflist in sequence, with commas between them
	evhm.write("%s\n" %motifscomma)
	evhm.write(" # -- heatmap --\n !heatmap\tmargin=1,colwidth=18,roundedcorner=1\n # -- show data value\n !showdataValue\tshow=0,fontsize=12,fontitalic=0,textalign=start\n\n")

	# Write the results
	for key in resultsdict:
		evhm.write("%s\t" %spdict[key]) #the label of the species: this is how it will be recognized in the tree by evolview
		evhm.write("%s\n" %resultsdict[key]) #the results of percentage calculations		
	evhm.close()


def makeevolviewbargraph(d_total,d_na,name):
	'''
	Takes a dictionary with hit results and uses corresponding species as labels.
	Returns a text file that can be read by evolview.
	'''
	# Open the heatmap text file:
	evbg = open("%s-%s.txt" %(bargraphtxt,name), "w")

	# Write the evolview headers
	evbg.write(" !groups\tnon-ambiguous motifs,ambiguous motifs\n !colors\tnavy,darkkhaki\n !itemheight\t15px\n !showdataValue\tshow=0,fontsize=12,fontitalic=0,textalign=end\n\n")

	# Write the results
	for key in d_total:
		t = int(d_total[key])
		na = int(d_na[key])
		a = t-na
		sp = spdict[key]
		evbg.write("%s\t%s,%s\n" %(sp,na,a))
	evbg.close()
		



#Make a dictionary for the species (abbreviation to full name)
species=open("%s/species.txt" %config.mainfolder)
spdict = {sp.strip().split()[0][0].lower()+sp.strip().split()[1][:3]: sp.strip() for sp in species}

#Make dictionary for the motif hit results (frequency and bargraph)
results,barstotal = {},{}
f = csv.reader(open(infile))
for n,line in enumerate(f):
	if n == 0:
		#these are the headers; extract the motifs in sequence
		motifs = line[1:]
	else:
		key = line[0].split('-')[-1]
		hits = line[1:]
		#calculate the total
		total = 0
		for h in hits:
			total += int(h)
		#calculate the percentages
		perc = []
		for h in hits:
			p = float(h)/total
			perc.append(str(p))
		value = ','.join(perc)
		#save in the motif hit results dictionary
		results[key] = value
		barstotal[key] = total


# For a double-category bar, get also the information from the nonambiguous counts
barsna = {}
f2 = csv.reader(open(infile_nonamb))
for n,line in enumerate(f2):
	if n == 0:
		continue #motifs are already extracted
	else:
		key = line[0].split('-')[-1]
		hits = line[1:]
		#calculate the total
		total = 0
		for h in hits:
			total += int(h)
		barsna[key] = total

makeevolviewheatmap(results,"occurrence")
makeevolviewbargraph(barstotal,barsna,"total")




