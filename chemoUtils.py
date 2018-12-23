'''
December 4, 2018

This is a module to accompany Project 4. It performs widely used functions including reading files,
computing a Tanimoto coefficient, computing a Tanimoto summary, and computing a bootstrap p-value. 
All functions are contained within the Tanimoto class.

Common usage (in main script): from chemoUtils import Tanimoto
'''

import sys
from collections import defaultdict
import itertools
import numpy

class Tanimoto(object):
	'''
	Class to store data from input files and perform computations relating to Tanimoto
	coefficients, summaries, and bootstrap p-values.
	'''
	def __init__(self):
		# Initialize storage structures
		self.fingerprints = {}
		self.enumerate_drugs = {}
		self.protein_to_drugs = defaultdict(list)
		self.drug_to_proteins = defaultdict(list)
		self.Tc_array = None
		self.pairwise_drugs = None
		self.all_ligands = None
		
	def read_files(self, drug_file, target_file):
		'''
		Converts input CSV files to dictionaries.
		Inputs:
			drug_file = CSV with columns db_id, generic_name, maccs
			target_file = CSV with columns db_id, uniprot_accession, uniprot_ID
		Returns:
			None
		'''
		count = 0
		drugs = open(drug_file, 'r')
		drugs.readline()
		for line in drugs:
			line = line.strip().split(',')
			self.fingerprints[count] = set(line[2].split()) # db_id int index: set, fingerprints
			self.enumerate_drugs[line[0]] = count # db_id: db_id int index
			count += 1
		drugs.close()

		drug_targets = open(target_file, 'r')
		drug_targets.readline()
		for line in drug_targets:
			line = line.strip().split(',')
			self.protein_to_drugs[line[1]].append(self.enumerate_drugs[line[0]]) # uniprot_accession: list, db_id int index
			self.drug_to_proteins[line[0]].append(line[1]) # db_id: list, uniprot_accession		
		drug_targets.close()
		
		# Fill an array for hashing Tanimoto coefficients (Tc) with value -1.0
		# This will be filled as Tc are computed to avoid repeating computations
		self.Tc_array = numpy.full( (len(self.fingerprints), len(self.fingerprints)), -1.0 )
		self.all_ligands = numpy.arange(len(self.fingerprints)) # Make sure ligand list is ordered for reproducibility

	def tanimoto_coefficient(self, d1, d2):
		'''
		Computes Tanimoto coefficient (Tc) between two drugs.
		Inputs:
			d1 = int, drug
			d2 = int, paired drug
		Returns:
			T_coefficient = float, Tc
		'''
		T_coefficient = 0.0
		
		# Sort drug indexes so only top half of Tc_array need be filled
		(smallest, largest) = sorted((d1, d2))
		
		# Check if this Tc has not been computed yet (== -1.0)
		# If not, compute Tc and store value in Tc_array
		if self.Tc_array[smallest][largest] == -1.0:
			d1_fingerprint = self.fingerprints[d1]
			d2_fingerprint = self.fingerprints[d2]
			T_coefficient = len(d1_fingerprint & d2_fingerprint) / float(len(d1_fingerprint | d2_fingerprint))
			self.Tc_array[smallest][largest] = T_coefficient

		# If already computed, take Tc directly from Tc_array			
		else:
			T_coefficient = self.Tc_array[smallest][largest]

		return T_coefficient

	def tanimoto_summary(self):
		'''
		Computes a Tanimoto summary between two proteins
		Inputs:
			None
		Returns:
			T_summary = float, Tanimoto summary
		'''
		T_summary = 0.0

		for d1, d2 in self.pairwise_drugs:
			T_coefficient = self.tanimoto_coefficient(d1, d2)
			# Add Tc to T_summary if Tc is above cutoff of 0.5
			if T_coefficient > 0.5:
				T_summary += T_coefficient

		return T_summary

	def compute_bootstrap_pvalue(self, n, proteinA, proteinB):
		'''
		Computes a bootstrap p-value for the Tanimoto summary of two proteins
		by random sampling with replacement for <n> samples.
		Inputs:
			n = number of iterations for sampling
			proteinA = a protein (uniprot_accession)
			proteinB = paired protein (uniprot_accession)
		Returns:
			b_pvalue/float(n) = float, bootstrap p-value
		'''
		b_pvalue = 0
		fuzzy_equals_diff = 10**(-16)

		# Generate pairwise combinations of drugs from the set for each protein
		self.pairwise_drugs = itertools.product(self.protein_to_drugs[proteinA], self.protein_to_drugs[proteinB])		
		# Compute T_summary for the proteins
		# allow for some precision errors for random samples with fuzzy_equals_diff
		real_T_summary = self.tanimoto_summary() - fuzzy_equals_diff

		# Will be sampling lists of the same size as each drug set from the set of all drugs
		n_A = len(self.protein_to_drugs[proteinA])
		n_B = len(self.protein_to_drugs[proteinB])

		# Sample (with replacement) and tally sampled T_summary that are greater than
		# or equal to the real T_summary. Then divide by number of sampling iterations.
		for i in range(n):
			random_A = numpy.random.choice(self.all_ligands, size=n_A)
			random_B = numpy.random.choice(self.all_ligands, size=n_B)
			self.pairwise_drugs = itertools.product(random_A, random_B)
			sample = self.tanimoto_summary()
			if sample >= real_T_summary:
				b_pvalue += 1

		return b_pvalue/float(n)
