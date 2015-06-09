import urllib2, os, config, sys


# species in pancompara: Arthropods: dmel,smar,turt,dpul,tcas


# spp list (in sequence)
spp = ['isca']#['cele']# ['dmel','smar','turt','dpul','tcas']


spp = ['dmel','tcas','dpul','smar','turt','atri','atha','crei','cmer','osat','ppat','smoe','slyc','aque','bmal','hrob','lgig','mlei','nvec','sman','spur','tadh','drer','ggal','hsap','mmus','xtro','bnat','ehux','glam','gthe','lmaj','pinf','pfal','tthe']


ctype = "pan_homology"
seqfolder = "%s/%s" %(config.mainfolder,config.seqfolder)


def filelist(spp):
	fastas = []
	for sp in spp:
		for f in os.listdir(seqfolder):
			if f.split('_')[-1] == 'seq.fasta':
				if sp in f:
					fastas.append(f)
	if len(fastas) != len(spp):
		sys.exit("ABORT:\nExtracting the filenames from the species list did not go as planned.\nLook \
at the 'filelist' function in the script to solve this problem.")
	return fastas

def fastaheaders(infile):
	'''
	Takes a fasta file with headers >geneID|genename|proteinID
	and returns a list of all geneIDs, all gene names, and all proteinIDs
	as well as a dictionary translating proteinID tot geneID.
	'''
	infile_o = open(infile)
	genes,names,proteins = [],[],[]
	prot2gene = {}
	for line in infile_o:
		if line[0] == ">":
			line = line.strip()
			gnp = line[1:].split('|')
			genes.append(gnp[0])
			names.append(gnp[1])
			proteins.append(gnp[2])
			prot2gene[gnp[2]] = gnp[0]
	return genes,names,proteins,prot2gene
		

def comparahtml(genename):
	'''
	From the geneID, retrieve the compara REST API url.
	See more information at ensemblgenomes.org/info/data/pan_compara
	'''
	url = "http://rest.ensemblgenomes.org/homology/id/%s?compara=%s&content-type=application/json" %(genename,ctype)
	response = urllib2.urlopen(url)
	html = response.read()
	return html

def parsecompara(html,sp,gene):
	#following four lines only to save file instead of parsing right away (so I can run it overnight)
	out = open("%s/compara/%s-%s.txt" %(config.mainfolder,sp,gene),"w")
	for line in html:
		out.write(line)
	out.close()

	orthoscsv = ""
	for sp in spp:
		orthoscsv += ""
	#orthoscsv = orthoscsv[:-1] #remove last comma
	return(orthoscsv)


if __name__ == "__main__":
	fastas = filelist(spp) #retrieve the input file (in this case fasta) to extract gene and protein IDs
	for k,sp in enumerate(spp):
		infile = "%s/%s" %(seqfolder,fastas[k])
		genes,names,proteins,prot2gene = fastaheaders(infile) #read fasta file 
		# open file, write headers
		#out = ""
		for gene in genes: # per gene
			print "Opening compara for %s in %s..." %(gene,sp)
			#out.write(gene) # write gene in output
			html = comparahtml(gene)
			#print html[0:10]
			orthoscsv = parsecompara(html,sp,gene)
			print "Parsing complete."
			# write orthos in output
