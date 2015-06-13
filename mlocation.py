import config,sys




if len(sys.argv) <= 1:
	sys.exit("USAGE: mlocation.py path/to/inputfile (input needs to be a motif sequence fasta file).\nPlease mind that the motif list and alphabet used to make this sequence file are unchanged for this analysis!")

infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]

### OPTION: only treat non-ambiguous zfs! ###
# set this to 0 if you want to treat all.
na = 1

# output




# start the dictionary containing the results
# results are collected per motif and for the total, and consist of a list of four elements:
# [total count] [found in first] [found in middle] [found in last]
resultdict = {m:[0,0,0,0] for m in config.motiflist}
resultdict['total'] = [0,0,0,0]

# dictionary to translate a letter back to a motif
mtranslate = {config.alphabet[a]: motif for a,motif in enumerate(config.motiflist)}


# make dictionary out of the input file
infi = open(infile)
fadict = config.fastadicter(infi)


def re2str(s):
	'''
	Takes a regular expression and turns it into a string of letters only.
	'''
	s = s.replace('{','')
	s = s.replace('}','')
	s = s.replace('|','')
	return s

def removere(s):
	p = s.count('{')
	for k in range(p):
		start = s.index('{')
		end = s.index('}') + 1
		s = s[:start] + s[end:]
	return s

def re2list(s):
	'''
	Takes a string that contains regular expression elements, and returns
	a list with the first, middle, and last elements separated, and all
	regular expression characters removed.
	'''
	#identify first
	if s[0] == '{':
		ends1=s.index('}') + 1
		if na == 1:
			s1 = ''
		else:
			s1 = re2str(s[:ends1])
		s = s[ends1:] #remove first element from s
	else:
		s1= s[0]
		s = s[1:] #remove first element from s
	#identify last
	if s[-1] == '}':
		starts2 = s.rfind('{')
		if na == 1:
			s2 = ''
		else:
			s2 = re2str(s[starts2:])
		s = s[:starts2] #remove last element from s
	else:
		s2 = s[-1]
		s = s[:-1]  #remove last element from s
	# now remove re elements from the leftover middle
	if na == 1:
		s = removere(s)
	else:
		s = re2str(s)
	reli = [s1,s,s2]
	return reli

def first(aa):
	global resultdict
	for a in aa:
		m = mtranslate[a]
		resultdict[m][0] += 1
		resultdict[m][1] += 1
		resultdict['total'][0] += 1
		resultdict['total'][1] += 1

def last(aa):
	global resultdict
	for a in aa:
		m = mtranslate[a]
		resultdict[m][0] += 1
		resultdict[m][3] += 1
		resultdict['total'][0] += 1
		resultdict['total'][3] += 1

def middle(aa):
	global resultdict
	for a in aa:
		m = mtranslate[a]
		resultdict[m][0] += 1
		resultdict[m][2] += 1
		resultdict['total'][0] += 1
		resultdict['total'][2] += 1




# start the analysis!
for protein in fadict:
	ssplit = fadict[protein].strip().split('Z')
	for s in ssplit:
		recount = s.count('{') + s.count('}') + s.count('|') * 2
		if len(s) - recount < 3:
			continue
		else:
			sli = re2list(s)
			first(sli[0])
			last(sli[2])
			middle(sli[1])

#print resultdict
#print "\n"

percdict = {}
for r in resultdict:
	resli = resultdict[r]
	try:
		percli = [int(100*rl/float(resli[0])) for rl in resli[1:]]
		percli.append(resli[0])
	except ZeroDivisionError:
		continue
	percdict[r] = percli


for p in percdict:
	if p == 'total':
		print 'total:', percdict[p]
	elif int(p.split('_')[1]) > 12:
		print p, percdict[p]
