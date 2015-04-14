import csv

in1 = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/nameclusters1.csv"))
in2 = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/nameclusters6.csv"))
out = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/evolview_newclusters.txt", "w")

actualnames = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/motifhits_all2.csv"))

namedict = {}
for name in actualnames:
	namedict[name[2]] = "%s|%s" %(name[1],name[2])

colourdict = {'1': '#6fc833', '2': '#ff75cf', '3': '#ff8c00', '4': '#ffd700', '5': '#1e3072', '6': '#8650ac', '7': '#0a6910', '8': '#5092ce', '9': '#0eb0ff', '10': '#bd1f30', '11': '#ffab91', '12': '#addbe1', '13': '#476291', '14': '#830000'}

out.write("## leaf background color\n\n")

"""
for line in in1:
	name = line[0].replace(':','-')
	colour = colourdict[line[1]]
	out.write("%s\t%s\tprefix\n" %(name,colour))
"""

for line in in2:
	#name = "%s|%s" %(line[1],line[2])
	#name = name.replace(':','-')
	name = namedict[line[2]]
	colour = colourdict[line[3]]
	out.write("%s\t%s\tprefix\n" %(name,colour))


out.close()

