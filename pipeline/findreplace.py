group = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/originalgroups.csv")
group2 = open("/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results/interest_clusters/originalgroups3.csv", "w")

for line in group:
	line = line.strip()
	line = line.replace("Completed file ","")
	group2.write(line)
	group2.write(",\n")




