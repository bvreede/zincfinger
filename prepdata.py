#!/usr/bin/python

'''
This script is used to prepare the system for a series of scripts in analysing zinc finger
domains of protein sequences.
Results folders are created, the protein datasets are unzipped, moved, and 
cleaned up by removing multiple isoforms. It retains the longest isoform.
HMMer is run on the resulting files.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 16 March 2016
'''
import os,config

out = open("%s/out.csv" %config.mainfolder, "w")

def findiso(fadict):
	'''
	A function used to collect the keys of the largest isoforms
	per gene, given a dictionary of sequences.
	'''
	genes = [] # a list of all the genes present in this database
	proteins = {} # a dictionary containing the geneID as key and the header of the longest isoform as value
	for key in fadict:
		try:
			geneID,genename,protein = key.split('|')
		except ValueError: #to prevent crash in case there are only two elements
			continue
		# check if there is already an isoform in the database
		if geneID in genes:
			# this gene already has an isoform in the database; check if the current isoform is longer.
			inDB_ID = proteins[geneID] # retrieve the header of the currently longest isoform for this gene
			inDB_seq = fadict[inDB_ID] # retrieve the sequence for this isoform
			cur_seq = fadict[key] # retrieve the sequence of the isoform currently being assessed
			if len(inDB_seq) < len(cur_seq):
				# the current isoform is longer; replace the entry in the database
				proteins[geneID] = key
		else:
			# this gene is not yet present in the database. Enter it!
			genes += [geneID]
			proteins[geneID] = key
	return proteins


def replace_db(infile):
	'''
	Takes a fasta file with sequences that may contain multiple isoforms
	per gene, and returns a fasta file with only the longest isoform.
	'''
	openfile = open(infile)	

	# retrieve the fasta entries as dictionary
	fadict = config.fastadicter(openfile)

	# remove the original file, open outputfile with same name
	openfile.close()
	os.remove(infile)
	outfile = open(infile,"w")

	# loop through the dictionary to collect protein isoforms (function)
	whichkeys = findiso(fadict)

	# print the new dictionary as fasta to the outputfile
	for key in whichkeys:
		header = whichkeys[key]
		sequence = fadict[header]
		outfile.write(">%s\n%s\n\n" %(header,sequence))
	outfile.close()



# prep system: create all folders
folders = [config.seqfolder,config.resfolder,config.imgfolder,config.dbfolder,config.orthfolder,config.compfolder,config.hmmfolder,config.seqpfolder,config.seqmfolder]
for f in folders:
	fn = config.mainfolder+'/'+f
	if os.path.exists(fn):
		continue
	else:
		os.system("mkdir %s" %fn)


# for all files in folder with ensembl downloads:
ensembldata="%s/%s" %(config.mainfolder,config.ensemblsource)
for enf in os.listdir(ensembldata):
	newname = "%s-%s_seq.fasta" %(config.idr,enf[:4])
	out.write("%s," %enf[:4])
	# copy the file and rename to [idr]-[abbr]_seq.fasta.gz
	os.system("cp %s/%s %s/%s.gz" %(ensembldata,enf,ensembldata,newname))
	# unzip the file and move it to sequences
	os.system("gunzip %s/%s.gz" %(ensembldata,newname))
	os.system("mv %s/%s %s/%s" %(ensembldata,newname,config.mainfolder,config.seqpfolder))
	# only save single isoform per gene
	replace_db("%s/%s/%s" %(config.mainfolder,config.seqpfolder,newname))
	# run HMMer
	command = "%s/%s/hmmsearch -o %s/%s/%s-%s_hmmsearch.txt --incT 3.0 %s/%s %s/%s/%s" %(config.mainfolder,config.hmmerbin,config.mainfolder,config.hmmfolder,config.idr,enf[:4],config.mainfolder,config.pfamc2h2,config.mainfolder,config.seqpfolder,newname)
	os.system(command)

