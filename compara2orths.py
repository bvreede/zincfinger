#!/bin/python
'''
This script uses the ensembl compara database to find orthologs among the sequences
in the fasta files of all species specified in config.sppall. It outputs a database
of orthologs inside the 'orthologs' folder, one database per species.
NB! Connection to the internet is necessary for this script to run, if 'parseonline'
or 'saving' (see Options) are set to 1!
'''


import urllib2, os, config, sys, re


### OPTIONS ###
saving = 0 # Set to 1 if you want to save the compara data locally as a text doc.
parseonline = 1 # Set to 1 if you want to parse ONLINE compara data
parselocal = 0 # Set to 1 if you want to parse LOCAL compara data (as a text file)

# spp lists
spp = [i for i in config.sppall if i not in config.chor]
vertebrates = config.chor #ENSURE THIS LIST CONTAINS ALL ENSEMBL VERTEBRATE SPECIES IN YOUR DB!

ctype = "pan_homology"
seqfolder = "%s/%s" %(config.mainfolder,config.seqfolder)
comparafolder = "%s/%s" %(config.mainfolder,config.compfolder)
resfolder = "%s/%s" %(config.mainfolder,config.orthfolder)

# regular expressions to search for in json results
proteinre = re.compile('protein_id.*taxon_id')
sppre = re.compile('species.*perc_id')


def filelist(spp):
	'''
	Gets the list of input files that are associated with the species in
	spp. NB! If there are multiple files associated with one species, the
	script aborts. Or if one of the species is missing...
	'''
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
	return genes,names,proteins,prot2gene #probably 'genes' is the only necessary component
		

def comparahtml(gene,sp):
	'''
	From the geneID, retrieve the compara REST API url.
	See more information at ensemblgenomes.org/info/data/pan_compara
	'''
	if sp in vertebrates:
		url = "http://rest.ensembl.org/homology/id/%s?compara=%s&content-type=application/json" %(gene,ctype) #vertebrate url
	else:
		url = "http://rest.ensemblgenomes.org/homology/id/%s?compara=%s&content-type=application/json" %(gene,ctype) #metazoa, plants, protists url
	response = urllib2.urlopen(url)
	html = response.read()
	return html

def savecompara(html,sp,gene):
	'''
	Instead of parsing right away, this function gives the option to save
	the page separately as a text document, so it can be parsed later and locally.
	'''
	out = open("%s/%s-%s.txt" %(comparafolder,sp,gene),"w")
	for line in html:
		out.write(line)
	out.close()

def findorthos(txt):
	'''
	Goes through compara json file to return orthologs in a dictionary
	where keys are proteins and values are the corresponding species.
	'''
	orthosdict = {}
	txtli = txt.split("source")[1:] #splits each ortholog hit (removing the first one, which contains text prior to the source statement)
	for source in txtli:
		# get separate entries for source and ortholog
		sourceentry,orthoentry = source.split("target")
		# in the entry for the ortholog, identify the species
		ospfind = sppre.search(orthoentry)
		ospecies = ospfind.group()[10:-10]
		osp_abbr = ospecies[0] + ospecies.split('_')[1][:3]
		if osp_abbr not in spp: # checks if we're interested in this species at all
			continue
		# in the entry for the source, identify the protein ID
		sprotfind = proteinre.search(sourceentry)
		sprotein = sprotfind.group()[13:-11]
		# in the entry for the ortholog, identify the protein ID	
		oprotfind = proteinre.search(orthoentry)
		oprotein = oprotfind.group()[13:-11]
		orthosdict[osp_abbr + '_' + oprotein] = sprotein
	return orthosdict


if __name__ == "__main__":
	fastas = filelist(spp) #retrieve the input file (in this case fasta) to extract gene and protein IDs
	for k,sp in enumerate(spp):
		print "Reading compara for %s." %sp
		infile = "%s/%s" %(seqfolder,fastas[k])
		genes,names,proteins,prot2gene = fastaheaders(infile) #read fasta file
		if parseonline ==1 or parselocal == 1:
			outfile = "%s/%s-allorth.csv" %(resfolder,sp)
			out = open(outfile, "w")
		for gene in genes: # per gene
			#print "Reading compara for %s in %s..." %(gene,sp)
			if saving == 1 or parseonline == 1:
				html = comparahtml(gene,sp)
			# save or parse?
			if saving == 1:
				savecompara(html,sp,gene)
			if parseonline == 1:
				orthosdict = findorthos(html) # to parse directly from online file
			if parselocal == 1:
				try:
					txtfile = open("%s/%s-%s.txt" %(comparafolder,sp,gene))
					txt = ""
					for line in txtfile:
						txt += line
					orthosdict = findorthos(txt) # to parse from a local file
				except IOError:
					print "No file found for %s in %s." %(gene,sp)
					continue
			# if parsing: write results to file
			if parseonline == 1 or parselocal == 1:
				for ortho in orthosdict:
					if ortho[:4] == sp: #paralog identified; not interesting for our purposes
						continue
					out.write("%s,%s,%s,%s\n" %(sp,orthosdict[ortho],ortho[:4],ortho[5:]))
		if parseonline ==1 or parselocal == 1:
			out.close
					
