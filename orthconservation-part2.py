import sys,config,csv,itertools


if len(sys.argv) <= 1:
	sys.exit("USAGE: python orthconservation-part2.py path/to/inputfile (input file is the '-detailed.csv' output of orthconservation.py).\nOutputfolders are indicated in the script; edit the script if you want to alter them.")

infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]

### OUTPUT FILES ###
heatmaptxt = "%s/%s/%s_comservation-heatmap" %(config.mainfolder,config.evfolder,infilebrev)

orthin = [line for line in csv.reader(open(infile))]

# make dictionaries where motif combos can be counted.
ambiguous = {} #dictionary to count appearance of ambiguous combinations
for m in config.motiflist:
	m_ = config.translationdict[m]
	for n in config.motiflist:
		n_ = config.translationdict[n]
		mn = frozenset([m_,n_])
		ambiguous[mn] = 0

# copy combination dictionaries
orthambi = dict(ambiguous) # same as ambiguous but this one to count ambiguous combinations that are orthologous
substitutions = dict(ambiguous) #same as ambiguous, to count substitutions between orthologs

# make dictionary where motifs individually can be counted
individual = {}
for m in config.motiflist:
	n = config.translationdict[m]
	individual[n] = 0

def re_move(s):
	'''
	returns a list with the frozenset-combinations of
	all different elements of a single regular expression.
	'''
	s = s.replace('{','')
	s = s.replace('}','')
	sli = s.split('|')
	sli_pairs = [frozenset(list(pair)) for pair in itertools.combinations(sli,2)]
	return sli_pairs

def makeevolviewheatmap(doublematrix,name):
	'''
	Takes a list of lists (heatmap data) and uses motifs as labels
	and returns a text file that can be read by evolview.
	'''
	# Open the heatmap text file:
	evhm = open("%s-%s.txt" %(heatmaptxt,name), "w")
	evhm.write(" #heatmap\n !legendTitle\tFrequency of overlap\n !showLegends\t1\n !colorgradient\tfloralwhite,orange,red,purple,blue\
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



for combo in orthin:
	o1 = config.re2li(combo[1])
	o2 = config.re2li(combo[3])
	for e,f in zip(o1,o2):
		if e == 'Z':
			continue
		# Determine orthology and frequency of ambiguous items
		elif e.count('{') > 0 or f.count('{') > 0: # ambiguous motif identified.
			# turn 
			eli = re_move(e)
			fli = re_move(f)
			# count each instance (in total)
			for item in eli + fli:
				try:
					ambiguous[item] += 1
				except KeyError:
					continue
			# determine if it is orthologous
			for item in eli:
				if item in fli:
					orthambi[item] += 2
		# with no ambiguous items: either score a substitution, or compare sequences...
		# score a substitution here:
		elif e != f:
			individual[e] += 1
			individual[f] += 1
			substitutions[frozenset([e,f])] += 1

#print individual
for s in substitutions:
	if substitutions[s] != 0:
		a,b = list(s)
		#print s, substitutions[s], individual[a], individual[b]

#MAKE THE HEATMAP FOR AMBIGUOUS APPEARANCE#
doublematrix1,doublematrix2 = [],[]
for m in config.motiflist:
	line1,line2 = [],[]
	for n in config.motiflist:
		a = config.translationdict[m]
		b = config.translationdict[n]
		combo = frozenset([a,b])
		line1.append(ambiguous[combo])
		if ambiguous[combo] > 50:
			freq = float(orthambi[combo])/ambiguous[combo]
		else:
			freq = 0
		line2.append(freq)
	doublematrix1.append(line1)
	doublematrix2.append(line2)
makeevolviewheatmap(doublematrix1,"ambiguous")
makeevolviewheatmap(doublematrix2,"orthambi")

totalambiguous,totalconserved = 0,0
for key in ambiguous:
	totalambiguous += ambiguous[key]
	totalconserved += orthambi[key]

print "ambiguous: %s, of which conserved: %s (%s" %(totalambiguous,totalconserved,int(float(totalconserved)/totalambiguous*100)) + "%)"
