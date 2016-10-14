# zincfinger


## 0. Before you start...
1. Download your data: protein sequences from Ensembl Biomart ([Ensembl](http://www.ensembl.org/biomart/martview/0fea912a097e984e080910d0e481bc04), [Metazoa](http://metazoa.ensembl.org/biomart/martview/41b7027ad573051d5ae6042363b08980), [Plants](http://plants.ensembl.org/biomart/martview/8e64fd02b779de836a8ad9e145548997), [Protists](http://protists.ensembl.org/biomart/martview/85097007194679c2b70acd1ceb4b668b)), ensuring that the fasta headers consist of (in this order) >geneID|genename|proteinID. Download as a zipped file, and make sure that the filename starts with a four character species specification. Optional: use the filter for SMART-domain 'SM00355' on the data prior to download.
2. Save all data in a single folder inside a main folder.
3. Download [HMMer](http://hmmer.org/) and the [C2H2 Pfam model](http://pfam.xfam.org/family/PF00096/hmm) (via [Pfam]( http://pfam.xfam.org/family/PF00096)). Save the binary and model in a second folder inside the main folder.
4. Customize **config.py** to contain the path to your main folder, the folder with downloaded ensembl data, the path to the hmm binaries and pfam model, and a list of the species you downloaded from ensembl (in 'sppall'). If you downloaded data from the main Ensembl (i.e. not Metazoa, Plants, or Protists), make sure they are added to the 'chor' list.
5. Run the script `prepdata.py` to prepare the protein sequences and runn HMMer. This script will clean up the data to contain only a single isoform per gene (the longest).

## 1. Find C2H2 motifs in fasta files
Use the script **findmotif-collected.py**. This script will use another script (**findmotif.py**) and functions on fasta files with protein sequences, where the headers consist of 'geneID|genename|proteinID'. It uses regular expressions to identify possible distinct C2H2 variants in the protein sequence, and cross-matches them with the HMMer results to toss out false positives. It will run the actual motif search twice; first to score the frequency of all motifs, and subsequently to identify a selection of motifs to use in further analysis.

_Customize:_ 
- The criteria for motif selection can be customized in the script.

_Usage:_ `findmotif-collected.py`

_Output:_
- In 'sequences/motifs':
  - fasta file of protein sequences translated to motif sequences ([fileidentifier]_hmmprotstring.fa).
- In 'data':
  - fasta file of selected proteins ([fileidentifier]_seq.fa).
  - fasta files with the sequences of all motifs found ([fileidentifier]_hmmallmotifs.fa)
  - individual fasta files per motif ([fileidentifier]-motseq-[motif].fa)
  - all motifs found translated to their string notation ([fileidentifier]_hmmallmotifs.txt), as a total collection of what was found, including ambiguous motifs.
  - a csv database of each protein and locations of detected motifs ([fileidentifier]_hmmhitsdb.csv)
- In 'images':
  - three bar graphs indicating how many motifs were found per motif type; 'stacked' splits up ambiguous and non-ambiguous; 'non-ambiguous_motifs' only shows non-ambiguous motifs.
  - a heatmap file which shows the frequency of overlap per motif.
- In 'results':
  - a csv database showing motif overlaps ([fileidentifier]_hmmmotifstats.csv)
  - a txt file that counts raw occurrences and non-ambiguous (i.e. non-overlapping, abbreviated as 'na') occurrences per motif, per original input fasta file ([fileidentifier]_motifcount.txt)

## 2. Identify orthologs between species
- _NB: this is only useful if in the previous step you have worked with multiple species..._
- Use the script **compara2orths.py**. This script mines ensembl compara (see http://ensemblgenomes.org/info/data/pan_compara) to identify orthologs between the proteins in your database. It saves the identified orthologs to per-species databases in the 'orthologs' folder.
- Before you start, ensure that for each species, only one fasta file exists in the 'sequences' folder (species are identified by a four-character specification; see 0.2).

_Customize:_ 
- NB! When customizing species lists, make sure the list 'chor' in **config.py** contains all the vertebrates in your species list (i.e. those species that ensembl stores in their vertebrate database).
- The script allows download of the compara data to local files; this can be useful if you want to run the script multiple times (mining the database online will take time). In this case, under the heading 'options', set 'saving' to 1 (NB, 'parselocal' should always be set to 0 if you haven't run the 'saving' option yet). Then, for further runs, set 'saving' and 'parseonline' to 0, and 'parselocal' to 1.
For normal runs, it is sufficient to have 'parseonline' to 1 and 'parselocal' and 'saving' to 0.

_Usage:_ `compara2orths.py`

_Output:_
- In 'orthologs':
  - csv databases (one for each species) with all orthologs of proteins of this species found in the other species in use.
- In 'compara' (only if option 'saving' is run):
  - a local copy of the data in the ensemb compara api.

## 3. Determine conservation between orthologs
NB: all the following scripts can be run as a unit with `orthconservation-all.py`.

### 3.1 Separate ortholog pairs into categories.
- To do this, use the script **orthconservation.py**. This script requires the python package 'jellyfish' to be installed.
- The script defines five categories of ortholog pairs:
  - those with identical motif structures
  - those with substitutions only
  - those with structural differences only
  - those with additions/deletions between the pair members
  - other differences (i.e. combinations of the above).
- In addition to separating ortholog pairs into these categories, the script also generates a random motif sequence with each pair (the random sequence is based on one of the members of the pair, but with randomized motifs — for details, see the accompanying paper) and runs the comparison again, now separating this random pair into one of the above categories.

_Customize:_
- You can customize the species selection in this script, which will cause the script to only use ortholog pairs from the species you selected. Currently, three selections are specified: 'arth', for all arthropod species; 'd700', for species with a minimum evolutionary distance of 700my, and 'sppall', for all species available. When you run the script, you will need to specify one of these groups. The default is 'sppall'.
- If you want to compare the results for various species combinations, the script needs to be run separately for every combination.

_Usage:_ `orthconservation.py [species selection]`

_Output:_
- The results will be printed in the terminal, as well as saved in an output file.
- Some pairs as identified by compara (see 2.) won't exist in your database. These pairs can't be categorized, and are not part of the analysis. These are printed on the terminal as 'not counted'.
- In the folder 'results':
  - A .csv file with the results as printed in the terminal (named _[identifier]-[speciescombination]_orthcomp.csv_)
  - A .csv file with names and motif sequences of the 'identical' and 'substitution' pairs (named _[identifier]-[speciescombination]_orthcomp-detail.csv_)
  - A .csv file with names and motif sequences of the 'identical' and 'substitution' pairs of the random model (named _[identifier]-[speciescombination]_orthcomp-detail-random.csv_)
- In the folder 'images':
  - Two .png files of piecharts, one for the ortholog data and one for the random model. 

### 3.2 Compare individual motifs
- To do this, use the script **orthconservation-part2.py**.
- This script uses the 'detailed' ('...-detail.csv') output files from the previous step to identify motifs that can be directly compared. This happens in two ways:
  - Motifs that remain unsubstituted between orthologs are compared amino acid by amino acid.
  - Motifs that are substituted between orthologs are categorized to determine which motifs get substituted.

_Customize:_ 
- In the amino acid comparison of unsubstituted motifs, it is possible to only regard ortholog pairs spanning a certain evolutionary distance. This can be specified in the 'list of species for comparison', by defining two groups of species that will only be compared between (i.e. members of group1 with members of group2) but not among (i.e. members of group1 with other members of group1) the groups.

_Usage:_ `orthconservation-part2.py path/to/inputfile`

_Output:_
- In 'results':
  - a bar graph source file that can be used to generate a bar graph (recognizable by 'bardata' in the name), assessing per motif how often it is found conserved (i.e. in direct comparisons the motif is unchanged), or non-conserved (i.e. different motifs are found in the same location between orthologs).
- In 'images':
  - heatmaps per motif where amino acid substitution is scored. Top row is substitution between amino acids of orthologous motifs; bottom row is comparisons between randomly paired motifs of the same kind. _NB: if this does not appear, check that the species you used are separated into group1 and group2; these comparisons are only made with species combinations where one is part of group1 and the other of group2._


## 4. Determine the order of motif types
- Use the script **mlocation.py** to determine whether certain motifs are found in the beginning, middle or end of a connected series.

_Usage:_ `mlocation.py species-set` (the species set is e.g. arth, for arthropodes, or sppall, for all species. Check the script for more options!)

_Output:_
- In 'images':
  - Four .svg images with bar charts plotting the spacing between key amino acids (total, C-C, C-H, or H-H) to the location of motifs with these properties: either in the beginning or at the end of connected series.
  - Three corresponding .csv files with the absolute counts of motifs that were used as the base of this chart (the chart shows percentages only).
 
- Printed on the terminal:
  - results from a chi-squared goodness-of-fit test for the locations of motifs: significant results indicate that the selection of motifs (e.g. 5 residues between CC) deviate significantly from the average distribution.
