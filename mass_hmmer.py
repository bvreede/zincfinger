#!/bin/python
import os

#individual groups
chor = ['drer','ggal','hsap','mmus','xtro'] #chor
arth = ['dmel','tcas','isca','dpul','smar','turt'] #arth
eani = ['nvec','mlei','tadh','sman','aque'] #eani
plan = ['atri','ppat','crei','atha','slyc','osat','smoe','cmer'] #plan
prot = ['pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam'] #prot
anim = ['hrob','spur','lgig','bmal','cele','drer','ggal','hsap','mmus','xtro','dmel','tcas','isca','dpul','smar','turt','nvec','mlei','tadh','sman','aque']

#all species
sppall = ['pfal','bnat','tthe','gthe','lmaj','ehux','pinf','glam','dmel','tcas','isca','dpul','smar','turt','atri','atha','crei','cmer','osat','ppat','smoe','slyc','aque','bmal','cele','hrob','lgig','mlei','nvec','sman','spur','tadh','drer','ggal','hsap','mmus','xtro']

#species in concatenated files
sppcombi = ['arth','chor','eani','anim','prot','plan']

'''
for sp in spp1:
	#command = "python findmotif2.py /home/barbara/Dropbox/shared_work/zinc_finger_data/newdata/sequences/150111-SM00355-%s_seq.fasta" %sp
	command = "/Users/BarbaraMaria/Downloads/hmmer-3.1b2-macosx-intel/binaries/hmmsearch -o ~/Dropbox/shared_work/zinc_finger_data/newdata/results/150111-SM00355-%s_hmmsearch.txt --incT 3.0 /Users/BarbaraMaria/Dropbox/shared_work/zinc_finger_data/newdata/sequences/zf-C2H2.hmm /Users/BarbaraMaria/Dropbox/shared_work/zinc_finger_data/newdata/sequences/150111-SM00355-%s_seq.fasta" %(sp,sp)
	os.system(command)
	#print "motifs identified for", sp
	print "hmm screen done for", sp
'''

count = 0
for sp in sppall:
	count += 1
	command = "python findmotif2.py /home/barbara/Dropbox/zinc_finger_data/newdata/sequences/150602-SM00355-%s_seq.fasta" %sp
	#command = "/Users/BarbaraMaria/Downloads/hmmer-3.1b2-macosx-intel/binaries/hmmsearch -o ~/Dropbox/zinc_finger_data/newdata/results/150602-SM00355-%s_hmmsearch.txt --incT 3.0 /Users/BarbaraMaria/Dropbox/zinc_finger_data/newdata/sequences/zf-C2H2.hmm /Users/BarbaraMaria/Dropbox/zinc_finger_data/newdata/sequences/150602-SM00355-%s_seq.fasta" %(sp,sp)
	os.system(command)
	print "%s: motifs identified for %s." %(count,sp)
	#print "hmm screen done for", sp
