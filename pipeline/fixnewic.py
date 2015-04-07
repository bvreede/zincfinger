infile = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/motifhits_allz.txt")
outfile = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/motifhits_allz2.txt", "w")


for line in infile:
	line = line.replace("dm:", "dm-")
	outfile.write(line)

outfile.close()
