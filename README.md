zincfinger
==========

*1. Find C2H2 motifs in fasta files*

Use the script 'findmotif.py'. Customize the code to find the source files in your system. This script will function on a (series of) fasta file(s) with protein sequences, where the headers consist of 'geneID|genename|proteinID'.

It needs three subfolders in the main folder: 

1. sequences; a folder for the input sequences, as well as the output sequences (it will create a fasta file with all hits per motif)
2. images; a folder where the heatmap (showing overlap in hits of different motifs) will be created
3. results; a folder for all other results: a csv database as well as strings describing the sequences of motifs found in each protein (in fasta format).

*2. Process C2H2 hits to check for GO enrichment with GOrilla*

!NB: only works for species that have proper GO annotation, and that can be checked with GOrilla (http://cbl-gorilla.cs.technion.ac.il/).

Use two scripts in succession: 'select4gorilla.py' and 'translate4gorilla.py'. The former selects for each motif those genes/proteins that have either 1. ONLY hits for this motif ('exclusive'), or 2. hits for this motif, irrespective of any other motifs in the same protein ('inclusive'). The script currently allows overlap (i.e. if motif A and motif B are both present in the same location, then the protein can still count as an exclusive for both motif A and B), but this can be removed in the script (indicated in the comments).

translate4gorilla.py is necessary as GOrilla does not accept all forms of gene/protein ID. This script requires a 'translation database', which contains the following columns:

Furthermore, translate4gorilla.py creates a text document, which can be used immediately as input for GOrilla.
