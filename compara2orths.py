import urllib2, os, config


# species in pancompara: Arthropods: dmel,smar,turt,dpul,tcas


# spp list (in sequence)
spp = ['dmel','smar','turt','dpul','tcas']
ctype = "pan_homology"
seqfolder = "%s/%s" %(config.mainfolder,config.seqfolder)


def filelist(spp):
	fastas = []
	for sp in spp:
		for f in os.listdir(seqfolder):
			if f.split('_')[-1] == 'seq.fasta':
				if sp in f:
					fastas.append(f)
	return fastas

def fastaheaders(infile):
	infile_o = open(infile)
	genes,names,proteins = [],[],[]
	for line in infile_o:
		if line[0] == ">":
			gnp = line[1:].split('|')
			genes.append(gnp[0])
			names.append(gnp[1])
			proteins.append(gnp[2])
	return genes,names,proteins
		

def comparaurl(genename):
	url = "http://rest.ensemblgenomes.org/homology/id/%s?compara=%s&content-type=application/json" %(genename,ctype)
	return url

def parsecompara(url):
	orthoscsv = ""
	for sp in spp:
		orthoscsv += ""
	orthoscsv = orthoscsv[:-1] #remove last comma
	return(orthoscsv)


if __name__ == "__main__":
	fastas = filelist(spp)
	for k,sp in enumerate(spp):
		print fastas[k]
		# read fasta file and only save headers
		#genes,names,proteins = fastaheaders(infile)
		# open file, write headers
		'''
		out = ""
		for gene in genes: # per gene
			out.write(gene) # write gene in output
			url = comparaurl(gene)
			orthoscsv = parsecompara(url)
			# write orthos in output
		'''
