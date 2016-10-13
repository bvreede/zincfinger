#!/usr/bin/python

'''
This script runs all orthconservation scripts as a unit.

Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 10 August 2015
'''

import os,config


#first run orthconservation.py
for s in ["d700",'arth','sppall']:
	os.system("orthconservation.py %s" %s)

#then run orthconservation-part2.py
resfolder = "%s/%s" %(config.mainfolder,config.resfolder)
orthfilelist = [f for f in os.listdir(resfolder) if "-detail.csv" in f]
ranfilelist = [f for f in os.listdir(resfolder) if "-detail-random.csv" in f]
infilelist  = orthfilelist+ranfilelist
for f in infilelist:
	os.system("orthconservation-part2.py %s" %f)
