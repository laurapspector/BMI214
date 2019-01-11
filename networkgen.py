'''
December 4, 2018

This script takes in a CSV file of drug fingerprints based on SMILES codes, a CSV file of those drugs with the 
proteins to which they bind, and a CSV file with a selection of proteins and the indications for which those
proteins are targeted by drugs. It generates a network edgelist for protein pairs with significant Tanimoto
summary values that can be read by the networkx package to draw a network.

Usage: python3 networkgen.py <drugs.csv> <targets.csv> <protein_nodes.csv>
'''

import sys
from chemoUtils import Tanimoto # import Tanimoto class
import itertools
from numpy import random

#### ------- CLASSES ------- ####
class NetworkGen(object):
	'''
	Object to compute and store an edgelist for a network.
	'''
	def __init__(self, node_file):
		'''
		Constructor
		Inputs:
			node_file = name of file with proteins to be included in network
		'''
		self.edgelist = []

		f = open(node_file, 'r')
		f.readline()
		nodes = [line.strip().split(',')[0] for line in f]
		f.close()
		
		# Generator with pairwise nodes (no self-pairs or reciprocal pairs)
		self.pairwise_nodes = itertools.combinations(nodes, 2)

	def compute_edges(self, my_tanimoto):
		'''
		Computes bootstrap p-value for protein pairs to determine whether they should have an edge
		between them. Uses default 500 samples.
		Inputs:
			my_tanimoto = instance of Tanimoto class
		Returns:
			None
		'''
		fuzzy_equals_diff = 10**(-16)

		for proteinA, proteinB in self.pairwise_nodes:
			b_pvalue = my_tanimoto.compute_bootstrap_pvalue(500, proteinA, proteinB)

			# if bootstrap p-value is less than or equal to 0.05, add protein pair to edgelist
			if b_pvalue <= (0.05 + fuzzy_equals_diff):
				self.edgelist.append([proteinA, proteinB])

	def write_to_file(self):
		'''
		Writes edgelist to file
		Inputs:
			None
		Returns:
			None
		'''
		output_file = open('network_edgelist.txt', 'w')
		for pair in self.edgelist:
			output_file.write(' '.join(pair))
			output_file.write('\n')
		output_file.close()

#### ------- MAIN METHODS ------- ####
def main():
	# Runs the script!
	my_tanimoto = Tanimoto()
	my_tanimoto.read_files(sys.argv[1], sys.argv[2])
	network = NetworkGen(sys.argv[3])
	# Set seed for pseudorandom number generator with default value
	random.seed(214)
	network.compute_edges(my_tanimoto)
	network.write_to_file()

if __name__=="__main__":

	main()

