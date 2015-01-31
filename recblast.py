


# prep databases by individually saving files
spp = "turt"
infasta = open("/home/barbara/data/zincfingers/%szfs.fa" %spp)
outfolder = "/home/barbara/data/zincfingers/" + spp

for line in infasta:
	line = line.strip()
	if len(line) == 0:
		ofile.write('\n')
		ofile.close()
	elif line[0] == '>':
		name = line[1:]
		ofile = open("%s/%s" %(outfolder,name), "w")
	else:
		ofile.write(line)
	
	

