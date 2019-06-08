'''
October 30, 2018

This script performs K-Nearest Neighbors classification using leave-one-out cross-validation.
It prints the sensitivity and specificity values to the console and writes the predicted
sample labels to an output file.

Usage: python3 knn.py expfile sampfile k p

'''
# Import modules
import sys
import pandas
from operator import itemgetter
import numpy as np

#### ------- CLASSES ------- ####

class DistanceMatrix(object):
	'''
	Stores a 2D array of size number samples x number samples with each element being
	the euclidean distance between the row and column sample.
	'''

	def __init__(self, numpy_array):
		'''
		Initialize objects for storing distance values
		Inputs:
			numpy_array = a 2D numpy array containing expression data, rows=genes and columns=samples
		'''
		self.ncol = np.size(numpy_array, 1)
		self.matrix = np.zeros( (self.ncol, self.ncol) )
		self.data_array = numpy_array

	def initialize_matrix(self):
		'''
		Pre-compute the euclidean distance from each sample to every other sample
		and store in the 2D array.
		Inputs:
			None
		'''
		for i in range(self.ncol):

			# Fill the diagonal with arbitrary large values to avoid having a sample be a closest neighbor with itself
			self.matrix[i][i] = 10**9

			# compute distances across upper half of the diagonal
			for j in range(i+1, self.ncol):
				first_array = self.data_array[:,i]
				second_array = self.data_array[:,j]

				# compute pair-wise distance between columns of the input array
				distance_array = first_array - second_array
				distance = np.sqrt(np.sum(np.square(distance_array)))

				# set value=distance at P_1 x P_2 and P_2 x P_1 in the new array
				self.matrix[i][j] = distance
				self.matrix[j][i] = distance


class KNNClassifier(object):
	'''
	Uses leave-one-out cross-validation and K-Nearest Neighbors classification to predict sample labels.
	'''
	def __init__(self, k, p):
		'''
		Initialize attributes for storing command line inputs and lists for storing true
		and predicted class labels.
		Inputs:
			k = number nearest neighbors
			p = fraction of nearest neighbors needed for positive classification
		'''
		self.k = k
		self.p = p
		self.y_true = []
		self.y_pred = []
		self.output_list = []

	def predict(self, distance_matrix, samples_dict, samples_list):
		'''
		Predicts test sample label using K-Nearest Neighbors algorithm.
		Inputs:
			distance_matrix = 2D array storing distances between all combinations of samples
			samples_dict = dictionary mapping sample ID to label (whether true or permuted)
			samples_list = ordered list of sample IDs from expression file
		'''
		# True labels, ordered as read in from the expression file
		self.y_true = [samples_dict[i] for i in samples_list]

		# Make sure the value of k from the command line is valid for leave-one-out cross-validation
		if self.k >= len(self.y_true): 
			self.k = len(self.y_true) - 1

		# Find k samples with smallest distance to the test sample
		# If greater than p of the k samples belong to the patient group,
		# Assign test sample to the patient group; otherwise to the healthy group
		for i in range(len(distance_matrix.matrix)):
			k_indices = np.argpartition(distance_matrix.matrix[i], self.k)[:self.k]
			k_labels = [samples_dict[samples_list[g]] for g in k_indices]
			num_patients = sum( k_labels )
			if (self.k * self.p) < num_patients:
				assignment = 1
			else:
				assignment = 0

			# Add predicted label to a new list to be used for generating confusion matrix
			self.y_pred.append(assignment)
			# Store the sample ID and its predicted label
			self.output_list.append( (samples_list[i], assignment))

	def write_assignments(self):
		''''
		Writes sample IDs and their predicted labels to an output file.
		Inputs:
			None
		Outputs:
			sample_assignments.txt = tab delimited file with sample ID and predicted label
				sorted by sample ID
		'''
		h=open("sample_assignments.txt", 'w')
		for i in sorted(self.output_list, key=itemgetter(0)):
			h.write(i[0] + '\t' + str(i[1]) + '\n')
		h.close()
		

class ConfusionMatrix(object):
	'''
	Stores a confusion matrix object for calculating sensitivity and specificity.
	'''
	def __init__(self, y_true, y_pred):
		'''
		Initializes confusion matrix, stores ordered lists of true and predicted labels
		Inputs:
			y_true = ordered list of true sample labels
			y_pred = ordered list of predicted sample labels from K Nearest Neighbors classifier
		'''
		self.cm = [[0,0], [0,0]]
		self.y_true = y_true
		self.y_pred = y_pred

		# Generate the confusion matrix, a 2x2 array
		# Format is [[TN, FP], [FN, TP]]
		for (true, pred) in zip(y_true, y_pred):
			self.update_cm(true, pred)

	def update_cm(self, true, pred):
		'''
		Updates a specific entry in the confusion matrix
		Inputs:
			true = true label for given sample
			pred = predicted label for given sample
		'''
		self.cm[true][pred] += 1

	def print_stats(self):
		'''
		Prints sensitivity and specificity values to the console
		Inputs:
			None
		Outputs:
			Prints sensitivity and specificity to the console
		'''
		sensitivity = self.cm[1][1] / float(self.cm[1][1] + self.cm[1][0])
		specificity = self.cm[0][0] / float(self.cm[0][0] + self.cm[0][1])
		print("Sensitivity: " + '{:.2f}'.format(sensitivity))
		print("Specificity: " + '{:.2f}'.format(specificity))

#### ------ USEFUL FUNCTIONS ------- ####

def read_expression_files(file_name, index_col):
	'''
	Reads in an expression file and saves it as a pandas DataFrame object
	Inputs:
		file_name = file name taken as a command line argument
		index_col = the column to assign as the row labels
	Outputs:
		Returns a pandas DataFrame object with the expression data, rows=genes and columns=samples
	'''
	f = open(file_name, 'r')
	file_dataframe = pandas.read_csv(f, sep='\t', index_col=index_col)
	f.close()
	return file_dataframe

#### ------- MAIN METHODS ------- ####

def main():
	'''
	Runs all the functions!
	'''
	expression_file = sys.argv[1]
	samples_file = sys.argv[2]
	k = int(sys.argv[3])
	p = float(sys.argv[4])

	# Read in sample info file and make dictionary mapping sample name to true label
	samples_dict = {}
	g = open(samples_file, 'r')
	input_lines = [line.strip() for line in g if line.strip()]
	for line in input_lines:
		line = line.split('\t')
		sample_name = line[0]
		sample_label = int(line[1]) # '0'=control or '1'=patient
		samples_dict[sample_name] = sample_label
	g.close()

	# Read in expression file and create a dataframe
	expression_df = read_expression_files(expression_file, 0)
	# Generate ordered list of samples from expression file
	samples_list = [i for i in expression_df.columns.values]
	# Convert dataframe to numpy 2D array
	expression_array = expression_df.values

	# Pre-compute the distance matrix
	distance_matrix = DistanceMatrix(expression_array)
	distance_matrix.initialize_matrix()

	# Generate a K-Nearest Neighbors classifier object
	knn = KNNClassifier(k, p)
	# Predict the labels for all test samples using leave-one-out cross-validation
	# and the K-Nearest Neighbors algorithm
	knn.predict(distance_matrix, samples_dict, samples_list)
	# Write predicted labels to an output file
	knn.write_assignments()

	# Generate a confusion matrix and print the sensitivity and specificity of our model
	confusion_matrix = ConfusionMatrix(knn.y_true, knn.y_pred)
	confusion_matrix.print_stats()


if __name__=="__main__":
    main()
