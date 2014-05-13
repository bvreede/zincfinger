'''
The purpose of this script is to change the names of
table headers so that spaces are removed and replaced by
underscores. This makes the table headers usable in SQL.
'''

import csv, sys

db = sys.argv[1]

#r = csv.reader(open(db))
r=open(db)
o = open(db[:-4] + "_hc.csv","w")

header=r.next().strip().split(',')

for h in header:
	newh = h.replace(' ','_')
	o.write(newh + ",")
o.write("\n")

for l in r:
	o.write(l)

o.close()


