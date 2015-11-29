# zincfinger

## 0. Before you start...
1. Set up your system. Create a main folder with the following subfolders: "sequences", "results", "images", "databases", "evolview", "compara", and "orthologs".
2. Download your data: protein sequences from Ensembl Biomart, ensuring that the fasta headers consist of (in this order) >geneID|genename|proteinID. Name them with a four character species specification, followed by an underscore (and don't use further underscores in the filename!). _e.g. 150602-SM00355-xtro_seq.fasta_
3. Download HMMer (http://hmmer.janelia.org/) and the C2H2 Pfam model (http://pfam.xfam.org/family/PF00096/hmm via http://pfam.xfam.org/family/PF00096) and run the binary hmmsearch on your downloaded Ensembl data. Save the result file in mainfolder/results, preferably [fileidentifier]_hmmsearch.txt. _e.g. 150602-SM00355-xtro_hmmsearch.fasta_
4. Customize **config.py** to contain the path to your main folder, the hmm result file, and a list of the species you downloaded from ensembl.

## 1. Find C2H2 motifs in fasta files

Use the script **findmotif.py**. This script will function on a fasta file with protein sequences, where the headers consist of 'geneID|genename|proteinID'. It uses regular expressions to identify possible distinct C2H2 variants in the protein sequence, and cross-matches them with the HMMer results to toss out false positives.

_Customize:_ 
- Add your HMMer results file to the script (called 'hmmfile'). You can opt to adjust this per species, or create a concatenated file with all HMMer results and call this from the script. The script is currently set to search for a file with the same identifier as the fasta file, and _hmmsearch.txt as suffix (see 0.3).
- Adjust the **config.py** file to determine the motifs the script will search for (called 'motiflist').

_Usage:_ findmotif.py path/to/fastafile

_Output:_
- In 'databases':
  - fasta files with the sequences of all motifs found ([fileidentifier]_hmmallmotifs.fa)
  - individual fasta files per motif
  - all motifs translated to their string notation ([fileidentifier]_hmmallmotifs.txt)
  - a csv database of each protein and locations of detected motifs ([fileidentifier]_hmmhitsdb.csv)
- In 'evolview':
  - a newick file with all motifs (to upload as a new "tree")
  - a heatmap file which shows the frequency of overlap per motif, can be used as heatmap annotation with the above newick file
- In 'images':
  - three bar graphs indicating how many motifs were found per motif type; 'stacked' splits up ambiguous and non-ambiguous; 'non-ambiguous_motifs' only shows non-ambiguous motifs.
  - two heatmaps showing motif overlap; 'singlenorm' only shows normalization over the x axis; 'doublenorm' shows normalization over both axes.
- In 'results':
  - a csv database showing motif overlaps ([fileidentifier]_hmmmotifstats.csv)
  - every time the script is run, two databases will be appended: hitcount_allspp.csv, showing the total numbers of hits per motif, with each file that was run as a new row, and hitcount_allspp-nonambg.csv, showing only non-ambiguous (i.e. non-overlapping) hits.

## 2. Identify orthologs between species

- Use the script **compara2orths.py**. This script mines ensembl compara to identify orthologs between the proteins in your database. It saves the identified orthologs to per-species databases in the 'orthologs' folder.
- Before you start, ensure that for each species, only one fasta file exists in the 'sequences' folder (species are identified by a four-character specification; see 0.2).

_Customize:_ 
- NB! When customizing species lists, make sure the list 'chor' in **config.py** contains all the vertebrates in your species list (i.e. those species that ensembl stores in their vertebrate database).
- The script allows download of the compara data to local files; this can be useful if you want to run the script multiple times (mining the database online can take time). In this case, under the heading 'options', set 'saving' to 1. Then, for further runs, set 'saving' and 'parseonline' to 0, and 'parselocal' to 1.
For normal runs, it is sufficient to have 'parseonline' to 1 and 'parselocal' and 'saving' to 0.

_Usage:_ findmotif.py path/to/fastafile

_Output:_
- In 'orthologs':
  - csv databases (one for each species) with all orthologs of proteins of this species found in the other species in use.
- In 'compara' (only if option 'saving' is run):
  - a local copy of the data in the ensemb compara api.

## 3. Determine conservation between orthologs

Use the script **orthconservation.py**
Use the script **orthconservation-part2.py**
Use the script **orthconservation-part3.py**

## 4. Determine the order of motif types

Use the script **mlocation.py**

## Scripts to check:
evolview_motifcount

