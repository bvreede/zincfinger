# zincfinger

## 0. Before you start...
1. Set up your system. Create a main folder with the following subfolders: "sequences", "results", "images", "databases", "evolview", and "orthologs".
2. Download your data: protein sequences from Ensembl Biomart, ensuring that the fasta headers consist of (in this order) >geneID|genename|proteinID.
3. Download HMMer (http://hmmer.janelia.org/) and run it on your downloaded Ensembl data. If you want to run it on multiple items, you can use **mass_hmmer.py** to speed up the process; make sure that the correct list of files/species names is in there!
4. Customize **config.py** to contain the path to your main folder, and a list of the species you downloaded from ensembl.

## 1. Find C2H2 motifs in fasta files

Use the script 'findmotif.py'. This script will function on a (series of) fasta file(s) with protein sequences, where the headers consist of 'geneID|genename|proteinID'.

It needs three subfolders in the main folder: 

1. sequences; a folder for the input sequences, as well as the output sequences (it will create a fasta file with all hits per motif)
2. images; a folder where the heatmap (showing overlap in hits of different motifs) will be created
3. results; a folder for all other results: a csv database as well as strings describing the sequences of motifs found in each protein (in fasta format).

*2. Process C2H2 hits to check for GO enrichment with GOrilla*

!NB: only works for species that have proper GO annotation, and that can be checked with GOrilla (http://cbl-gorilla.cs.technion.ac.il/).

Use two scripts in succession: 'select4gorilla.py' and 'translate4gorilla.py' (don't forget to customize paths). The former selects for each motif those genes/proteins that have either 1. ONLY hits for this motif ('exclusive'), or 2. hits for this motif, irrespective of any other motifs in the same protein ('inclusive'). The script currently allows overlap (i.e. if motif A and motif B are both present in the same location, then the protein can still count as an exclusive for both motif A and B), but this can be removed in the script (indicated in the comments).

translate4gorilla.py is necessary as GOrilla does not accept all forms of gene/protein ID. This script requires a 'translation database', which contains the following columns:
ENSEMBL Gene ID; Associated Gene Name; ENSEMBL Protein ID; RefSeq Protein ID; RefSeq mRNA; UniProt/SwissProt ID.
The database can easily be constructed using Ensembl BioMart.

Finally, translate4gorilla.py creates a text document, which can be used immediately as input for GOrilla.

In case all putative translations for the protein are empty, the gene symbol will be noted instead. However, all symbols will be noted only at the end of the document, so that they can be potentially left out of a GOrilla analysis. Symbols appear under the header 'SYMBOLS:' in the text document.
