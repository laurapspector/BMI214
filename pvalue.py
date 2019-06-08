'''
December 4, 2018

This script generates a bootstrap p-value for the Tanimoto summary between two proteins. It takes as input
a CSV file of drug fingerprints based on SMILES codes, a CSV file of those drugs with their protein targets,
and two proteins (uniprot_accession ID). 

Optional arguments: '-n' [number of samples for bootstrap] '-r' [seed for pseudorandom number generator]

Usage: python3 pvalue.py -n <INT> -r 214 <drugs.csv> <targets.csv> <proteinA> <proteinB>
'''
import sys
import argparse
from chemoUtils import Tanimoto # import Tanimoto class
from numpy import random

#### ------- USEFUL FUNCTIONS ------- ####
def fill_parser(parser, inputs):
	'''
	Parses command line arguments, positional and optional
	Inputs:
		parser = an argparse ArugmentParser object
		inputs = a list of command line arguments
	Returns:
		parser.parse_args() = ArgumentParser object with arguments stored in Namespace object
	'''
	parameter_list = ['-n', '-r', 'drug_file', 'target_file', 'proteinA', 'proteinB']
	parameter_default = [500, 214, None, None, None, None]
	parameter_type = [int, int, str, str, str, str]

	for parameter, p_default, p_type in zip(parameter_list, parameter_default, parameter_type):
		parser.add_argument(parameter, default=p_default, type=p_type)

	return parser.parse_args()

#### ------- MAIN METHODS ------- ####
def main():
	# Runs the script!
	parser = argparse.ArgumentParser()
	args = fill_parser(parser, sys.argv)
	my_tanimoto = Tanimoto()
	my_tanimoto.read_files(args.drug_file, args.target_file)
	# Set seed for pseudorandom number generator
	random.seed(args.r)
	# Print to the console the bootstrap pvalue for given protein pair
	print(my_tanimoto.compute_bootstrap_pvalue(args.n, args.proteinA, args.proteinB))

if __name__=="__main__":

	main()