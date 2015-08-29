# zincfinger

## 0. Before you start...
1. Set up your system. Create a main folder with the following subfolders: "sequences", "results", "images", "databases", "evolview", and "orthologs".
2. Download your data: protein sequences from Ensembl Biomart, ensuring that the fasta headers consist of (in this order) >geneID|genename|proteinID. Name them with a four character species specification, followed by an underscore (and don't use further underscores in the filename!). _e.g. 150602-SM00355-xtro_seq.fasta_
3. Download HMMer (http://hmmer.janelia.org/) and the C2H2 Pfam model (http://pfam.xfam.org/family/PF00096/hmm via http://pfam.xfam.org/family/PF00096) and run HMMer on your downloaded Ensembl data. Save the result file in mainfolder/results.
4. Customize **config.py** to contain the path to your main folder, the hmm result file, and a list of the species you downloaded from ensembl.

## 1. Find C2H2 motifs in fasta files

Use the script **findmotif.py**. This script will function on a fasta file with protein sequences, where the headers consist of 'geneID|genename|proteinID'. It uses regular expressions to 

## 2. Identify orthologs between species

Use the script **compara2orths.py**.

## 3. Determine conservation between orthologs

Use the script **orthconservation.py**
Use the script **orthconservation-part2.py**
Use the script **orthconservation-part3.py**

## 4. Determine the order of motif types

Use the script **mlocation.py**

## Scripts to check:
evolview_motifcount
fmcompare
lireg

