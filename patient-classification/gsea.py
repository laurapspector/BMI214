'''
October 30, 2018

This script takes as input an expression file for multiple patients, a sample file indicating
their true group assignments, and a file of KEGG gene sets. It performs gene set enrichment
analysis (GSEA), writes the true enrichment scores for each KEGG gene set to a file, and prints
the number of significantly enriched pathways. It uses 100 permutations to generate the 
background distribution.

Usage: python3 gsea.py expfile sampfile keggfile

'''
# Import modules
import sys
import pandas
import numpy as np
import math
from operator import itemgetter

#### ------- CLASSES ------- ####
class ExpressionData(object):
	'''
	Reads in and stores the expression data file as a pandas DataFrame object.
	'''
	def __init__(self):
		# Initialize the dataframe variable and a list for holding the ordered
		# list of patients from the input file.
		self.df = None
		self.gene_list = []

	def read_data(self, file_name):
		'''
		Reads in an expression file with rows=patient IDs and columns=genes
		and creates a pandas DataFrame. Saves patient IDs in order to a list.
		Inputs:
			file_name = expression file name (a command line argument)
		'''
		f = open(file_name, 'r')
		file_dataframe = pandas.read_csv(f, sep='\t', index_col=0)
		f.close()

		# Update the class attributes
		self.df = file_dataframe
		self.gene_list = list(file_dataframe.index)

class EnrichmentScore(object):
	'''
	Performs operations needed to compute enrichment score for each gene set.
	Also writes enrichment scores to an output file.
	'''
	def __init__(self):
		# Initialize objects to hold numpy arrays for vectorized operations
		# and a dictionary for mapping genes to rank.
		self.DE_array = None
		self.ranked_gene_dict = {}
		self.ES_array = None

	def rank_genes(self, df, samples_dict):
		'''
		Computes fold change in expression between patients and samples
		and ranks genes by their fold change.
		Inputs:
			df = pandas DataFrame containing gene expression data
			samples_dict = dictionary mapping patient ID to true label
		'''
		# Initialize empty pandas Series for storing sum of expression values
		# and computing differential expression
		control_array = pandas.Series(0., index=df.index)
		patient_array = pandas.Series(0., index=df.index)

		# Get N for each population
		n_controls = len(samples_dict) - sum(samples_dict.values())
		n_patients = sum(samples_dict.values())
		
		# Iterate over DataFrame, sum a series with the initialized Series that
		# corresponds to its group label (0 or 1)
		for key, value in df.iteritems():
			if samples_dict[key] == 0: # control
				control_array += value
			else: # patient
				patient_array += value

		# Get the mean for each Series
		control_array /= n_controls
		patient_array /= n_patients

		# Compute fold change and rank genes
		# Note that expression values are already log2 normalized
		self.DE_array = patient_array - control_array
		self.DE_array = self.DE_array.sort_values(ascending=False)

		# Fill dictionary that maps a gene name to its rank
		self.ranked_gene_dict = {i[1]:(i[0]) for i in list(enumerate(self.DE_array.index))}

	def calculate_ES(self, kegg_dict, kegg_list):
		'''
		Performs a Brownian bridge and computes the enrichment score
		as the positive supremum of the Brownian bridge.
		Inputs:
			kegg_dict = dictionary mapping KEGG gene set name to genes contained within
			kegg_list = ordered list of KEGG gene set names from keggfile
		'''
		# Initialize 1D array to store erichment scores for each gene set
		self.ES_array = np.zeros(len(kegg_dict))

		# Loop through KEGG gene sets, compute the up_step and down_step penalities
		# and compute a running sum for a Brownian bridge until the supremum is reached.
		N_T = len(self.ranked_gene_dict)
		for key, value in kegg_dict.items():
			N_G = len(value)
			supremum = 0
			brownian_bridge = 0
			up_step = math.sqrt( (N_T-N_G)/float(N_G) )
			down_step = -1 * math.sqrt( float(N_G)/(N_T-N_G) )
			prev_rank = 0

			# Sort genes by rank; for each ranked gene, add to running sum 
			# down_step * number of genes between this and previous 
			# ranked gene (rank-wise) plus 1 * up_step
			for gene_rank in sorted( [self.ranked_gene_dict[i] for i in value] ):
				brownian_bridge += ((gene_rank-prev_rank)*down_step + up_step)
				prev_rank = gene_rank + 1
				if brownian_bridge > supremum:
					supremum = brownian_bridge

			# Store the supremum in ES_array
			self.set_ES(supremum, kegg_list.index(key))

	def set_ES(self, supremum, row):
		'''
		Stores the supremum for a given KEGG gene set in an array.
		Inputs:
			supremum = positive supremum of the Brownian bridge
			row = row of array to set as supremum, rows are ordered according to input file
		'''
		self.ES_array[row] = supremum

	def write_ES(self, kegg_list):
		'''
		Writes list of enrichment scores to an output file.
		Inputs:
			kegg_list = ordered list of KEGG gene set names from keggfile
		Outputs:
			kegg_enrichment_scores.txt = tab delimited file with ranked gene sets and their
				enrichment scores
		'''
		output_file=open('kegg_enrichment_scores.txt', 'w')
		for key, value in sorted( zip(kegg_list,self.ES_array), key=itemgetter(1), reverse=True):
			output_file.write('\t'.join([key, '{:.2f}'.format(value) +'\n'])) # check that this works
		output_file.close()

class EnrichmentScoreArray(object):
	'''
	Stores a 2D array used to compute the nominal p-value after generating the 
	background distribution.
	'''
	def __init__(self):
		# Initialize an array object for storing enrichment scores, dictionary for mapping 
		# KEGG gene sets to their p-values, and a running sum of significant gene sets.
		self.array = None
		self.p_values = {}
		self.sig_gene_sets = 0

	def initialize_array(self, array):
		'''
		Initializes the 2D array with a 1D array consisting of the KEGG gene sets'
		true enrichment scores.
		Inputs:
			array = array of true enrichment scores (ES_array from EnrichmentScore object)
		'''
		self.array = array

	def save_array(self, new_array):
		'''
		Appends a new array onto the existing array. Each column contains ordered
		enrichment scores (according to input file, not rank) for a different gene
		set with the first column being the true enrichment scores and subsequent 
		columns corresponding to enrichment scores after sample label permutation.
		Inputs:
			new_array = an ES_array from the EnrichmentScore object
		'''
		self.array = np.vstack( (self.array, new_array) )

	def transpose_array(self):
		'''
		Transposes the array with all enrichment scores for ease of computing p-values.
		Simply returns array formed with numpy vstack function to rows=gene sets and
		columns=permutations.
		'''
		self.array = np.transpose(self.array)

	def calculate_p_value(self, kegg_list, n_permutations):
		'''
		Calculates nominal p-value across each gene set. Performs Bonferroni correction
		and computes the number of significant pathways.
		Inputs:
			kegg_list = ordered list of KEGG gene set names from keggfile
		Outputs:
			Prints to the console the number of significant pathways.
		'''
		fuzzy_equals_diff = 10**(-6)
		for row in self.array:
			num_greater = len(np.where(row[1:] > (row[0] - fuzzy_equals_diff))[0])
			nominal_p_value = num_greater/float(n_permutations)

			# Perform Bonferroni correction
			corrected_p_value = nominal_p_value * len(self.array)
			if corrected_p_value < 0.05:
				self.sig_gene_sets += 1

		print("Number significant pathways = " + str(self.sig_gene_sets))


class PermutationTest(object):
	'''
	Permutes sample labels and stores permuted sample labels for iterating.
	Achieved by permuting the sample IDs and assigning them to the non-permuted
	sample labels, then storing the new association in a dictionary.
	'''
	def __init__(self, expr_df, samples_dict, n_permutations):
		'''
		Generates list of dictionaries to mirror format of samples_dict
		where each patient ID maps to a new, random label.
		Inputs:
			expr_df = the expression file DataFrame
			samples_dict = dictionary mapping patient IDs to their labels
			n_permutations = number of permutations
		'''
		self.all_permutations = []
		ID_list = list(expr_df) # ordered list of sample IDs
		label_list = [samples_dict[i] for i in ID_list] # ordered list of sample labels
		self.permutations = []
		for i in range(n_permutations):
			np.random.shuffle(ID_list) # shuffle the sample IDs
			# assign randomly ordered sample ID a label by zipping together list of
			# randomly ordered sample IDs and ordered sample labels
			permuted_labels_dict = dict( zip(ID_list, label_list) )
			self.all_permutations.append(permuted_labels_dict)

#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.

    Inputs:
        a = a number
        b = another number
    Returns:
        True or False
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)

#### ------- MAIN METHODS ------- ####
def main():
	'''
	Runs all the functions!
	'''
	expression_file = sys.argv[1]
	samples_file = sys.argv[2]
	kegg_file = sys.argv[3]
	n_permutations = 100
	
	# Save patient IDs and true labels to a dictionary
	samples_dict = {}
	g = open(samples_file, 'r')
	input_lines = [line.strip() for line in g if line.strip()]
	for line in input_lines:
		line = line.split('\t')
		sample_name = line[0]
		sample_label = int(line[1]) # '0'=control or '1'=patient
		samples_dict[sample_name] = sample_label
	g.close()

	# Read in expression file
	exprfile = ExpressionData()
	exprfile.read_data(expression_file)

	# Generate list of all genes to exclude genes from KEGG gene sets that
	# do not also exist in the list of genes with expression data
	all_genes = exprfile.gene_list

	# Generate dictionary to map KEGG gene sets to the genes contained within
	# and a list of ordered KEGG gene set names
	kegg_dict = {}
	kegg_list = []
	h = open(kegg_file, 'r')
	for line in h:
		if line.strip():
			line = line.strip().split('\t')
			kegg_list.append(line[0])
			kegg_dict[line[0]] = [i for i in line[2:] if i in all_genes]
	h.close()

	# Compute true enrichment scores
	ES_matrix = EnrichmentScoreArray()
	enrichment_score = EnrichmentScore()
	enrichment_score.rank_genes(exprfile.df, samples_dict)
	enrichment_score.calculate_ES(kegg_dict, kegg_list)
	ES_matrix.initialize_array(enrichment_score.ES_array)

	# Write ranked gene sets and enrichment scores to file
	enrichment_score.write_ES(kegg_list)

	# Generate background distribution
	permuted_labels = PermutationTest(exprfile.df, samples_dict, n_permutations)
	for i in permuted_labels.all_permutations:
	 	permuted_enrichment_score = EnrichmentScore()
	 	permuted_enrichment_score.rank_genes(exprfile.df, i)
	 	permuted_enrichment_score.calculate_ES(kegg_dict, kegg_list)
	 	ES_matrix.save_array(permuted_enrichment_score.ES_array)
	ES_matrix.transpose_array()

	# Compute corrected p-values and print number of significant pathways
	ES_matrix.calculate_p_value(kegg_list, n_permutations)

if __name__=="__main__":
    main()
