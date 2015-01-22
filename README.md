zincfinger
==========

*1. Find C2H2 motifs in fasta files*

Use the script 'findmotif.py'. Customize the code to find the source files in your system. This script will function on a (series of) fasta file(s) with protein sequences, where the headers consist of '[geneID]|[genename]|[proteinID]'.

It needs three subfolders in the main folder: 

1. sequences; a folder for the input sequences, as well as the output sequences (it will create a fasta file with all hits per motif)
2. images; a folder where the heatmap (showing overlap in hits of different motifs) will be created
3. results; a folder for all other results: a csv database as well as strings describing the sequences of motifs found in each protein (in fasta format).
