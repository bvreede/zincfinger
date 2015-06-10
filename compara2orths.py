import urllib2, os, config, sys


### OPTIONS ###
saving = 1 # Set to 1 if you want to save the compara data locally as a text doc.
parse-online = 0 # Set to 1 if you want to parse ONLINE compara data
parse-local = 0 # Set to 1 if you want to parse LOCAL compara data (as a text file)

# spp list (in sequence)
spp = ['hsap','mmus','xtro'] # ['isca','dmel','tcas','dpul','smar','turt','cele','atri','atha','crei','cmer','osat','ppat','smoe','slyc','aque','bmal','hrob','lgig','mlei','nvec','sman','spur','tadh','bnat','ehux','glam','gthe','lmaj','pinf','pfal','tthe','mmus','drer','ggal','hsap','mmus','xtro']

vertebrates = ['mmus','drer','ggal','hsap','mmus','xtro']

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
		

def comparahtml(gene,sp):
	'''
	From the geneID, retrieve the compara REST API url.
	See more information at ensemblgenomes.org/info/data/pan_compara
	'''
	if sp in vertebrates:
		url = "http://rest.ensembl.org/homology/id/%s?compara=%s&content-type=application/json" %(gene,ctype) #vertebrate url
	else:
		url = "http://rest.ensemblgenomes.org/homology/id/%s?compara=%s&content-type=application/json" %(genename,ctype) #metazoa, plants, protists url
	response = urllib2.urlopen(url)
	html = response.read()
	return html

def savecompara(html,sp,gene):
	'''
	Instead of parsing right away, this function gives the option to save
	the page separately as a text document, so it can be parsed later and locally.
	'''
	out = open("%s/compara/%s-%s.txt" %(config.mainfolder,sp,gene),"w")
	for line in html:
		out.write(line)
	out.close()


def parsecompara(html):
	orthoscsv = ""
	for sp in spp:
		orthoscsv += ""
	orthoscsv = orthoscsv[:-1] #remove last comma
	return(orthoscsv)


if __name__ == "__main__":
	fastas = filelist(spp) #retrieve the input file (in this case fasta) to extract gene and protein IDs
	for k,sp in enumerate(spp):
		infile = "%s/%s" %(seqfolder,fastas[k])
		genes,names,proteins,prot2gene = fastaheaders(infile) #read fasta file 
		# open file, write headers
		#out = ""
		for gene in genes: # per gene
			print "Reading compara for %s in %s..." %(gene,sp)
			#out.write(gene) # write gene in output
			if saving == 1 or parse-online == 1:
				html = comparahtml(gene,sp)
			# save or parse?
			if saving == 1:
				savecompara(html,sp,gene)
				print "Saved file."
			if parse-online == 1:
				orthoscsv = parsecompara(html) # to parse directly from online file
				print "Parsed online data."
			if parse-local == 1:
				orthoscsv = parselocal(sp,gene)
				print "Parsed local data."
			if parse-online == 1 or parse-local == 1:
				# write orthos in output
