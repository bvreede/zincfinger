import os

infile = "/home/barbara/Dropbox/shared_work/zinc_finger_data/supplementary_files/Additional_file_2.fasta"
openinfile = open(infile)

fname = infile.split('/')[-1].split('.')[0] #the name of the input file without extension
outfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/sequences/%s" %fname

if os.path.exists(outfolder):
	print "Outputfolder already exists."
else:
	os.system("mkdir %s" %outfolder)

for line in openinfile:
	line = line.strip()
	if len(line) == 0:
		ofile.write('\n')
		ofile.close()
	elif line[0] == '>':
		try:
			print "Completed file " + name
			ofile.write('\n')
			ofile.close()
			name = line[1:].replace('/','-')
			name = name.replace('(','')
			name = name.replace(')','')
			ofile = open("%s/%s" %(outfolder,name), "w")
		except NameError:
			name = line[1:]
			ofile = open("%s/%s" %(outfolder,name), "w")
	else:
		ofile.write(line)

