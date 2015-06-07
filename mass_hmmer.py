#!/bin/python
import os

spp = ['dmel2','tcas2','dpul2','isca2','smar2','turt2','all2']

for sp in spp:
	command = "findmotif.py /home/barbara/Dropbox/shared_work/zinc_finger_data/newdata/sequences/150111-SM00355-%s_seq.fasta" %sp
	#command = "/Users/BarbaraMaria/Downloads/hmmer-3.1b2-macosx-intel/binaries/hmmsearch -o ~/Dropbox/shared_work/zinc_finger_data/newdata/results/150111-SM00355-%s_hmmsearch.txt --incT 3.0 /Users/BarbaraMaria/Dropbox/shared_work/zinc_finger_data/newdata/sequences/zf-C2H2.hmm /Users/BarbaraMaria/Dropbox/shared_work/zinc_finger_data/newdata/sequences/150111-SM00355-%s_seq.fasta" %(sp,sp)
	os.system(command)



