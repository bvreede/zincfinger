import sys,config,csv,itertools


if len(sys.argv) <= 1:
	sys.exit("USAGE: python orthconservation-part2.py path/to/inputfile (input file is the '-detailed.csv' output of orthconservation.py).\nOutputfolders are indicated in the script; edit the script if you want to alter them.")

infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]


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
#print ambiguous
for a in ambiguous:
	if ambiguous[a] != 0:
		print a, ambiguous[a], orthambi[a], str(int(float(orthambi[a])/ambiguous[a] * 100))+'%'
#print orthambi
	
	
# then check if they are conserved (per mot combo: yes or no) -- go over both lists, so yes orthologs will be counted double!
# 
