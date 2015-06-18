# open all ortholog files and import combinations as frozensets
# add option to use only limited species here
# make it into a set, so that duplicates are removed

# generate 10 categories:
	# Ortholog-identical
	# Ortholog-substitution
	# Ortholog-structure
	# Ortholog-addition/subtraction
	# Ortholog-other
	# Random-identical
	# Random-substitution
	# Random-structure
	# Random-addition/subtraction
	# Random-other

# for each ortholog combination:
	# retrieve motif sequence
	# generate random sequence based on randomly picked motif sequence + species-specific motif frequencies
	# for both sequence combinations:
		# calculate lowest levenshtein distance
		# if distance is 0, the motif sequences are identical
			# add to 'identical', and continue
		# if lengths of sequences and number of Z (and their indices) are the same, substitution explains the difference
			# add to 'substitution'
			# if this was ortholog combo, save their names for further processing.
			# continue
		# remove Z and calculate levenshtein distance again
		# if distance is 0, structure explains the difference
			# add to 'structure', and continue
		# if distance is equal to difference in length, addition of motifs explains the difference
			# add to add/subtr, and continue
		# else:
			# add to other/combination

# For ortholog and for random: make a pie chart with the results (or just output them as numbers and manually make a pie chart, whatever)
