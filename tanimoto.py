'''
December 4, 2018

This script takes in a CSV file of drug fingerprints based on SMILES codes and CSV file of those drugs with the 
proteins to which they bind and generates an output file of Tanimoto coefficients between all pairs of drugs as 
well as an indication whether each pair has a shared protein target (output file name is a required argument).

Usage: python3 tanimoto.py <drugs.csv> <targets.csv> <outputfile.csv>
'''
import sys
from chemoUtils import Tanimoto # import Tanimoto class
import itertools

#### ------- USEFUL FUNCTIONS ------- ####
def all_Tc(tanimoto):
	'''
	Generates list of all pairwise Tanimoto coefficients and indication whether drug pairs share a protein target.
	Inputs:
		tanimoto = instance of Tanimoto class from chemoUtils module
	Returns:
		output_list = list for writing to file with format 'drug1, drug2, Tanimoto coefficient, shared/not shared'
	'''
	# Generator of all pairwise drugs (no self-pairs or reciprocal pairs)
	all_pairwise_drugs = itertools.combinations(tanimoto.enumerate_drugs.keys(), 2)	

	output_list = []
	for d1, d2 in all_pairwise_drugs:
		# Convert db_id to int index for computing T_coefficient
		d1_int = tanimoto.enumerate_drugs[d1]
		d2_int = tanimoto.enumerate_drugs[d2]
		T_coefficient = tanimoto.tanimoto_coefficient(d1_int, d2_int)
		shared = '0'
		if len( set(tanimoto.drug_to_proteins[d1]) & set(tanimoto.drug_to_proteins[d2]) ) > 0:
			shared = '1'
		output_list.append([d1, d2, '{:.6f}'.format(T_coefficient), shared])

	return output_list

def write_scores(output_file, output_list):
	'''
	Writes returned value of all_Tc to a file.
	Inputs:
		output_file = output file name
		output_list = list with drugs pairs, Tanimoto coefficients and shared/not shared status
	Returns:
		None
	'''
	out = open(output_file, 'w')

	for i in output_list:
		out.write(','.join(i))
		out.write('\n')

	out.close()

#### ------- MAIN METHODS ------- ####
def main():
	# Runs the script!
	my_tanimoto = Tanimoto()
	my_tanimoto.read_files(sys.argv[1], sys.argv[2])
	output_list = all_Tc(my_tanimoto)
	write_scores(sys.argv[3], output_list)

if __name__=="__main__":

	main()
	