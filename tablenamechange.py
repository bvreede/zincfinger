'''
The purpose of this script is to change the names of
table headers so that spaces are removed and replaced by
underscores. This makes the tables usable in SQL.
'''

import csv, sys, os

os

db = sys.argv[1]

r = csv.reader(open(db))
header=r.next()
for title in header:
	newh = title.replace(' ','_')
	print newh

