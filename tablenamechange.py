'''
The purpose of this script is to change the names of
table headers so that spaces are removed and replaced by
underscores. This makes the table headers usable in SQL.
'''

import sys

db = sys.argv[1]
r=open(db)
o = open(db[:-4] + "_hc.csv","w")

header=r.next().strip().split(',')	# first line contains table headers; converted into list
headline = ""
for h in header:			
	newh = h.replace(' ','_')	# replaces spaces with underscores
	headline += newh + ","		# converts list back to csv, in string
o.write(headline[:-1])
o.write("\n")

for l in r:
	o.write(l)			# includes the rest of the file in the output
o.close()


