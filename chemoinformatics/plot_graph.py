'''
December 4, 2018

This script draws a networkx network using a provided edgelist file. It takes as arguments a network edgelist, 
a CSV file with a selection of proteins that will be nodes in the network, and a file path to a network figure.

Usage: python3 plot_graph.py <network_edgelist.txt> <protein_nodes.csv> <file_path_to_output_figure.png>
'''
import sys, pandas, networkx, matplotlib

#### ------- CLASSES------- ####
class Graph(object):
	'''
	Object to store, create, and modify the network graph.
	'''
	def __init__(self, args):
		'''
		Constructor. Reads in a edgelist to a networkx graph object and a CSV file
		with columns uniprot_accession, uniprot_id, indications into a pandas dataframe.
		Inputs:
			args = list of command line arguments
		Returns:
			None
		'''
		with open(args[1], 'rb') as edgelist:
			self.graph = networkx.read_edgelist(edgelist)

		with open(args[2], 'r') as proteins:
			self.protein_df = pandas.read_csv(proteins, sep=',', index_col=0)

		self.output_file = args[3]			

	def make_graph(self):
		'''
		Relabels graph nodes by uniprot_ID, re-colors nodes, draws the graph, and saves it to file.
		Inputs:
			None
		Returns:
			None
		'''
		# Relabel nodes from uniprot_accession (edgelist value) to uniprot_id
		df_dict = self.protein_df.to_dict()['uniprot_id']
		mapping = {key:value for key,value in df_dict.items() if key in self.graph.nodes}
		graph_relabeled = networkx.relabel_nodes(self.graph, mapping)

		# Draw graph while re-coloring nodes
		networkx.draw_networkx(graph_relabeled, node_color=self.color_nodes(), font_size=10)

		# Set size of figure and save to file
		matplotlib.pyplot.gcf().set_size_inches(8,8)
		matplotlib.pyplot.savefig(self.output_file, dpi=150)

	def color_nodes(self):
		'''
		Generates a list of colors for re-coloring nodes by 'indication'
		Inputs:
			None
		Returns:
			color_list = list, ordered colors for re-coloring nodes of graph
		'''
		color_dict = {'bp': 'red', 'bp;cholesterol': 'green', 'bp;cholesterol;diabetes': 'blue', 'bp;diabetes': 'purple'}
		color_list = [color_dict[self.protein_df.loc[node, self.protein_df.columns[1]]] for node in self.graph.nodes]
		return color_list

#### ------- MAIN METHODS ------- ####
def main():
	# Runs the script!
	graph = Graph(sys.argv)
	graph.make_graph()

if __name__=="__main__":

	main()
