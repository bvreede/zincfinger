#!/bin/python
import os

spp1 = ['dmel2','tcas2','dpul2','isca2','smar2','turt2','arth2']



spp2 = ['atri','atha','crei','cmer','osat','ppat','smoe','slyc','aque','bmal','cele','hrob','lgig','mlei','nvec','sman','spur','tadh','drer','ggal','hsap','mmus','xtro']

spp2 = ['bnat','ehux','glam','gthe','lmaj','pinf','pfal','tthe']


spp2 = ['dmel','tcas','isca','dpul','smar','turt']


sppall = ['dmel','tcas','isca','dpul','smar','turt','atri','atha','crei','cmer','osat','ppat','smoe','slyc','aque','bmal','cele','hrob','lgig','mlei','nvec','sman','spur','tadh','drer','ggal','hsap','mmus','xtro']


'''
for sp in spp1:
	#command = "python findmotif2.py /home/barbara/Dropbox/shared_work/zinc_finger_data/newdata/sequences/150111-SM00355-%s_seq.fasta" %sp
	command = "/Users/BarbaraMaria/Downloads/hmmer-3.1b2-macosx-intel/binaries/hmmsearch -o ~/Dropbox/shared_work/zinc_finger_data/newdata/results/150111-SM00355-%s_hmmsearch.txt --incT 3.0 /Users/BarbaraMaria/Dropbox/shared_work/zinc_finger_data/newdata/sequences/zf-C2H2.hmm /Users/BarbaraMaria/Dropbox/shared_work/zinc_finger_data/newdata/sequences/150111-SM00355-%s_seq.fasta" %(sp,sp)
	os.system(command)
	#print "motifs identified for", sp
	print "hmm screen done for", sp
'''


for sp in sppall:
	command = "python findmotif2.py /home/barbara/Dropbox/zinc_finger_data/newdata/sequences/150602-SM00355-%s_seq.fasta" %sp
	#command = "/Users/BarbaraMaria/Downloads/hmmer-3.1b2-macosx-intel/binaries/hmmsearch -o ~/Dropbox/zinc_finger_data/newdata/results/150602-SM00355-%s_hmmsearch.txt --incT 3.0 /Users/BarbaraMaria/Dropbox/zinc_finger_data/newdata/sequences/zf-C2H2.hmm /Users/BarbaraMaria/Dropbox/zinc_finger_data/newdata/sequences/150602-SM00355-%s_seq.fasta" %(sp,sp)
	os.system(command)
	print "motifs identified for", sp
	#print "hmm screen done for", sp
