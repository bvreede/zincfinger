'''
The purpose of this script is to change the names of
table headers so that spaces are removed and replaced by
underscores. This makes the table headers usable in SQL.
Furthermore, the order of the data is checked (and altered
if necessary) so that it is suitable for the rest of the
data processing pipeline.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 2 December 2014
'''

import sys, csv

'''
INPUT PARAMETERS
taken from the command line and subsequently modified
'''
if len(sys.argv) <= 1:
	sys.exit("USAGE: python table_prep.py path/to/inputfolder")
db = sys.argv[1]


r = csv.reader(open(db, "U")) # U = universal newline support. Sometimes these databases have \r instead of \n, so this way only the first line can be read.
o = open(db[:-4] + "_corr.csv","w")

'''
Copy the current header into memory and replace spaces with _
'''
header=r.next()#.strip().split(',')		# first line contains table headers; converted into list
curheader = []
for h in header:			
	curheader.append(h.replace(' ','_'))	# replaces spaces with underscores

'''
Header should be:
Gene_stable_ID,Gene_name,Protein_stable_ID,Pfam_ID,SMART_ID,GO_domain,GO_term_accession,GO_term_name
Determine the sequence of these elements in the current database by comparing
newheader to actheader.
ACTUAL = what the list is supposed to look like
CURRENT = what the list now looks like.
Write a list of numbers that represent the order of the ACTUAL elements in
the order of the CURRENT elements.
Report an error in case any element in the CURRENT content is not present in
the ACTUAL list, and vice versa. In case an ACTUAL element is not present
in the CURRENT list, simply fill out empty content in its place.
'''
actheader = ['Gene_stable_ID', 'Gene_name', 'Protein_stable_ID', 'Pfam_ID', 'SMART_ID', 'GO_domain', 'GO_term_accession', 'GO_term_name']
headseq = []
headline = ""
for name in actheader:
	try:
		i = curheader.index(name) #finds the ACTUAL header element in the CURRENT header and returns the index
	except ValueError:
		print "%s column not found in database. Continuing anyway..." %name
		i = 'N'
	headseq.append(i)
	headline += name + ',' # write the actual header for the final database
o.write("%s\n" %headline[:-1]) # and write it to file


for name in curheader: # check if any CURRENT columns are not present in ACTUAL (and therefore not used). Report if need be.
	if name in actheader:
		pass
	else:
		print "Unidentified column (%s) found in database. Continuing anyway..." %name

'''
Use the header sequence to shuffle all columns to the appropriate order
in the outputfile. Place an empty column in case of an N. Also do a final
check on the csv fields (such as: no commas allowed!)
'''
for line in r:
	nl = ''
	for n in headseq:
		try:
			k = line[n].replace(',',';') # sometimes, mysteriously, commas make their way into csv fields.
			nl += k + ','
		except TypeError: # if the desired column does not exist, the headseq line contains a letter
			nl += ',' # and the column in the new file simply stays empty.
	o.write("%s\n" %nl[:-1])
o.close()

