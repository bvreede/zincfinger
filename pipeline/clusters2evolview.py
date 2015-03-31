import csv

in1 = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/nameclusters1.csv"))
in2 = csv.reader(open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/nameclusters4.csv"))
out = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/evolview_oldclusters.txt", "w")

colourdict = {'1': '#6fc833', '2': '#ff75cf', '3': '#ff8c00', '4': '#ffab91', '5': '#476291', '6': '#8650ac', '7': '#0a6910', '8': '#5092ce', '9': '#bd1f30', '10': '#0eb0ff', '11': '#ffd700', '12': '#addbe1', '13': '#1e3072', '14': '#830000'}

out.write("## leaf background color\n\n")

"""
for line in in1:
	name = line[0].replace(':','-')
	colour = colourdict[line[1]]
	out.write("%s\t%s\tprefix\n" %(name,colour))
"""

for line in in2:
	name = "%s|%s" %(line[1],line[2])
	name = name.replace(':','-')
	colour = colourdict[line[3]]
	out.write("%s\t%s\tprefix\n" %(name,colour))

out.close()

