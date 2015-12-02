# zincfinger

## 0. Before you start...
1. Set up your system. Create a main folder with the following subfolders: "sequences", "results", "images", "databases", "evolview", "compara", "hmm", and "orthologs".
2. Download your data: protein sequences from Ensembl Biomart, ensuring that the fasta headers consist of (in this order) >geneID|genename|proteinID. Name them with a four character species specification, followed by an underscore (and don't use further underscores in the filename!). _e.g. 150602-SM00355-xtro_seq.fasta_. Ensure the first part of the identifier (in this case _150602-SM00355_) is the same with all databases you want to use. Save the databases in the 'sequences' folder.
3. Download HMMer (http://hmmer.janelia.org/) and the C2H2 Pfam model (http://pfam.xfam.org/family/PF00096/hmm via http://pfam.xfam.org/family/PF00096) and run the binary hmmsearch on your downloaded Ensembl data. Save the result file in mainfolder/hmm, preferably [fileidentifier]_hmmsearch.txt. (e.g. 150602-SM00355-xtro_hmmsearch.fasta). 
4. Customize **config.py** to contain the path to your main folder, the hmm result file, the identifier from step 2, and a list of the species you downloaded from ensembl in 'sppall'. If you downloaded vertebrates, make sure they are added to the 'chor' list.

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

(Only useful if in the previous step you have worked with multiple species...) Use the script **compara2orths.py**. This script mines ensembl compara (see http://ensemblgenomes.org/info/data/pan_compara) to identify orthologs between the proteins in your database. It saves the identified orthologs to per-species databases in the 'orthologs' folder.
- Before you start, ensure that for each species, only one fasta file exists in the 'sequences' folder (species are identified by a four-character specification; see 0.2).

_Customize:_ 
- NB! When customizing species lists, make sure the list 'chor' in **config.py** contains all the vertebrates in your species list (i.e. those species that ensembl stores in their vertebrate database).
- The script allows download of the compara data to local files; this can be useful if you want to run the script multiple times (mining the database online will take time). In this case, under the heading 'options', set 'saving' to 1 (NB, 'parselocal' should always be set to 0 if you haven't run the 'saving' option yet). Then, for further runs, set 'saving' and 'parseonline' to 0, and 'parselocal' to 1.
For normal runs, it is sufficient to have 'parseonline' to 1 and 'parselocal' and 'saving' to 0.

_Usage:_ compara2orths.py

_Output:_
- In 'orthologs':
  - csv databases (one for each species) with all orthologs of proteins of this species found in the other species in use.
- In 'compara' (only if option 'saving' is run):
  - a local copy of the data in the ensemb compara api.

## 3. Determine conservation between orthologs
### 3.1 Separate ortholog pairs into categories.
To do this, use the script **orthconservation.py**. This script requires the python package 'jellyfish' to be installed.
- The script defines five categories of ortholog pairs:
  - those with identical motif structures
  - those with substitutions only
  - those with structural differences only
  - those with additions/deletions between the pair members
  - other differences (i.e. combinations of the above).
- In addition to separating ortholog pairs into these categories, the script also generates a random motif sequence with each pair (the random sequence is based on one of the members of the pair, but with randomized motifs — for details, see the accompanying paper) and runs the comparison again, now separating this random pair into one of the above categories.

_Customize:_
- Specify the **identifier** ('idr') prior to running the script! This is the text preceding your species identification in the original database (see 0.2). For example, with input file _150602-SM00355-xtro_seq.fasta_ the 'idr' text should be _150602-SM00355_. If you don't do this, the script won't be able to locate your files.
- You can customize the species selection in this script, which will cause the script to only use ortholog pairs from the species you selected. The best way to do this is to specify your selection in **config.py** and to refer the variable 'spp' to that selection.
- If you want to compare the results for various species combinations, the script needs to be run separately for every combination.
- In the "options" section, you can adjust the output name as well. Don't forget to do this with each different species combination!

_Usage:_ orthconservation.py

_Output:_
- The results will be printed in the terminal, as well as saved in an output file.
- Some pairs as identified by compara (see 2.) won't exist in your database. These pairs can't be categorized, and are not part of the analysis. These are printed on the terminal as 'not counted'.
- In the folder 'results':
  - A .csv file with the results as printed in the terminal (named _[identifier]-[speciescombination]_orthcomp.csv_)
  - A .csv file with names and motif sequences of the 'identical' and 'substitution' pairs (named _[identifier]-[speciescombination]_orthcomp-detail.csv_)
  - A .csv file with names and motif sequences of the 'identical' and 'substitution' pairs of the random model (named _[identifier]-[speciescombination]_orthcomp-detail-random.csv_)

### 3.2 Compare individual motifs
To do this, use the script **orthconservation-part2.py**.
- This script uses the 'detailed' (both 'detail.csv' and 'detail-random.csv') output files from the previous step to identify motifs that can be directly compared. This happens in two ways:
  - Motifs that remain unsubstituted between orthologs are compared amino acid by amino acid.
  - Motifs that are substituted between orthologs are categorized to determine which motifs get substituted by which.

_Customize:_ 
- In the amino acid comparison of unsubstituted motifs, it is possible to only regard ortholog pairs spanning a certain evolutionary distance. This can be specified in the 'list of species for comparison', by defining two groups of species that will only be compared between (i.e. members of group1 with members of group2) but not among (i.e. members of group1 with other members of group1) the groups.

_Usage:_ orthconservation-part2.py path/to/inputfile (run separately for each file generated in the previous step!)

_Output:_
- In 'evolview':
  - a heatmap source file that can be used as heatmap annotation with the motif "tree" (generated during step 1), showing which motif types are evolutionarily substituted by which.
  


### 3.3 Do something else
- Use the script **orthconservation-part3.py**

## 4. Determine the order of motif types
- Use the script **mlocation.py**

## Scripts to check:
evolview_motifcount

