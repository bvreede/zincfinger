#!/bin/bash

gunzip $1

python ~/github/zincfinger/tablenamechange.py ${1%.gz}

gzip ${1%.gz}


