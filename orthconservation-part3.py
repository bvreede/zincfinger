import sys,config,csv#,itertools,random
#import pylab as pl
#import numpy as np

if len(sys.argv) <= 1:
	sys.exit("USAGE: python orthconservation-part3.py path/to/inputfile (input file is the '-detail.csv' output of orthconservation.py).")


### INPUT FILES ###
infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]
orthin = [line for line in csv.reader(open(infile))] #this opens the file with ortholog combinations to investigate

### OUTPUT FILE(S) ###
outfile = "%s/%s/%s_orthcomp-motifcons.csv" %(config.mainfolder,config.resfolder,infilebrev)
out = open(outfile, "w")

#scorekeeping dictionaries: score_con for conserved combinations, and score_non for non-conserved combinations
score_con = {m:0 for m in config.motiflist}
score_con['ambi'] = 0
score_non = {m:0 for m in config.motiflist}
score_non['ambi'] = 0


#go over all motif combinations
for combo in orthin:
	c1,c2 = config.re2li(combo[1]),config.re2li(combo[3])
	for m1,m2 in zip(c1,c2):
		#don't count spaces
		if m1 == 'Z':
			continue
		#m1 and m2 at this point are the two orthologous motifs.
		#check first for ambiguous combinations:
		if m1[0] == '{' or m2[0] == '{':
			#remove regular expression elements to score each motif individually
			m1_s, m2_s = config.re2str(m1),config.re2str(m2)

			#determine first if ambiguity is conserved
			ambscore = 0
			for e in m1_s:
				if e in m2_s:
					ambscore += 1
			if ambscore >= 2: #2 or more ambiguous elements are found in both orthologs and thus scored as conserved
				score_con['ambi'] += 1
			else:
				score_non['ambi'] += 1

			#m1 to m2:
			for e in m1_s:
				eT = config.translationdict_inv[e]
				if e in m2_s:
					score_con[eT] +=1
				else:
					score_non[eT] +=1
			#m2 to m1:
			for e in m2_s:
				eT = config.translationdict_inv[e]
				if e in m1_s:
					score_con[eT] +=1
				else:
					score_non[eT] +=1
		else:
			mot1 = config.translationdict_inv[m1]
			mot2 = config.translationdict_inv[m2]
			if mot1 == mot2:
				score_con[mot1] += 2
			else:
				score_non[mot1] += 1
				score_non[mot2] += 1

def percentage(a,b):
	if a+b == 0:
		return 0
	else:
		return float(a)/(a+b)*100


#calculate total numbers of motifs
mottot = 0
for m in config.motiflist:
	mottot += (score_con[m] + score_non[m])

threshold = mottot/200

print "Standard motifs:"
out.write("Standard motifs:\n")
for m in config.standard:
	if score_con[m] + score_non[m] > threshold:
		print m,'\t',score_con[m],'\t',score_non[m],'\t',percentage(score_con[m],score_non[m])
		out.write("%s,%s,%s,%s\n" %(m,score_con[m],score_non[m],percentage(score_con[m],score_non[m])))

print "\nAlternative motifs:"
out.write("\nAlternative motifs:\n")
for m in config.alternative:
	if score_con[m] + score_non[m] > threshold:
		print m,'\t',score_con[m],'\t',score_non[m],'\t',percentage(score_con[m],score_non[m])
		out.write("%s,%s,%s,%s\n" %(m,score_con[m],score_non[m],percentage(score_con[m],score_non[m])))

print "\nAmbiguous:"
out.write("\nAmbiguous:\n")
print 'ambi:',score_con['ambi'],'\t',score_non['ambi'],'\t',percentage(score_con['ambi'],score_non['ambi'])
out.write("Ambi,%s,%s,%s\n" %(score_con["ambi"],score_non["ambi"],percentage(score_con["ambi"],score_non["ambi"])))

out.close()



