#!/bin/bash

gunzip $1
unz=${1%.gz}

python ~/github/zincfinger/tablenamechange.py $unz

gzip $unz			#rezips the previously unzipped file
unz2=${unz%.csv}"_hc.csv"	#name of file changed by python script
gzip $unz2			#zips new file
