import config,sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['xtick.labelsize']=11

if len(sys.argv) <= 1:
	sys.exit("USAGE: mlocation.py path/to/inputfile (input needs to be a motif sequence fasta file).\nPlease mind that the motif list and alphabet used to make this sequence file are unchanged for this analysis!")

infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]

### OPTION: only treat non-ambiguous zfs! ###
# set this to 0 if you want to treat all.
na = 1

# output range
output = range(7,18)

lengths_of_series = []
lengths_of_proteins = []

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
	prot = protein.replace('Z','')
	prl = prot.count('{') + prot.count('}') + prot.count('|') * 2
	protl = len(prot) - prl
	lengths_of_proteins.append(protl)
	for s in ssplit:
		recount = s.count('{') + s.count('}') + s.count('|') * 2
		if len(s) - recount < 2:
			continue
		else:
			lengths_of_series.append(len(s) - s.count('{') - s.count('}') - s.count('|') * 2)
			sli = re2list(s)
			first(sli[0])
			last(sli[2])
			middle(sli[1])

percdict = {}
for r in resultdict:
	resli = resultdict[r]
	try:
		percli = [100*rl/float(resli[0]) for rl in resli[1:]]
		percli.append(resli[0])
	except ZeroDivisionError:
		continue
	percdict[r] = percli



# now assemble the dataset:
bardata = {}
for n in range(18):
	cfi,cmi,cen,total = 0,0,0,0
	for p in percdict:
		if p == 'total':
			continue
		CC,CH,HH = p.split('_')
		if int(CH) == n:
			c1 = percdict[p][0]
			cm = percdict[p][1]
			c2 = percdict[p][2]
			ct = float(percdict[p][3])
			cfi += c1*ct
			cmi += cm*ct
			cen += c2*ct
			total += ct
	if total == 0:
		bardata[n] = [0,0,0,0]
	else:
		bardata[n] = [cfi/total,cmi/total,cen/total,total]

# get data for plots:
set1,set2,set3,totals = [],[],[],[]
for n in output:
	set1.append(bardata[n][0])
	set2.append(bardata[n][1])
	set3.append(bardata[n][2])
	totals.append(bardata[n][3])


#print resultdict

"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""

n_groups = len(set1)

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.35

rects1 = plt.bar(index, set1, bar_width,alpha=0.7, color='c',label='First')
#rects2 = plt.bar(index + bar_width, set2, bar_width,alpha=opacity, color='g',label='Middle')
rects3 = plt.bar(index + bar_width+0.05, set3, bar_width,alpha=0.5,color='r',label='Last')

plt.xlabel('C-H distance')
plt.ylabel('Frequency')
plt.title('Appearance of motif in series')
plt.xticks(index + bar_width, (list(output)))
#plt.legend()

plt.tight_layout()

def autolabel(rects):
    # attach some text labels
    for n,rect in enumerate(rects):
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%totals[n],
                ha='center', va='bottom')

#autolabel(rects1)
#autolabel(rects2)

name = "CH"
outtotals = open("%s/%s/%s_%stotals.csv" %(config.mainfolder,config.imgfolder,infilebrev,name), "w")
print name,set1,set2,set3

for n,t in enumerate(totals):
	outtotals.write("%s,%s,%s,%s,%s\n" %(output[n],t,set1[n],set2[n],set3[n]))
outtotals.close()

plt.savefig("%s/%s/%s_%s.svg" %(config.mainfolder,config.imgfolder,infilebrev,name))

#print max(lengths_of_series), max(lengths_of_proteins)

