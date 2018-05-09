import k_clique_networkX
import mcl
import networkx as nx
import time as tt
import re
import math

en_clust_scores = []
en_clust_from = []

def have_that_edge(G, n1, n2):
	if G.has_edge(n1,n2):
		return True
	return False

def graph_updater(G, c):
	for node1 in c:
		for node2 in c:
			if have_that_edge(G, node1, node2):
				G.remove_edge(node1, node2)	
	return G


def starter(filename):
	f = open(filename, 'r')

	data = f.read()
	f.close()
	proteins = re.findall(r"[.\w]+", data)
	print(len(proteins))

	indexed_proteins = {}

	index = 0

	for pr in proteins:
	    if indexed_proteins.get(pr) == None:
	        indexed_proteins[pr] = index
	    
	        index = index + 1

	print("Total Nodes: ", str(len(indexed_proteins))) 

	data = data.split('\n')
	data = data[:len(data)-1]
	new_data = list()
	for element in data:
		add = element.split('\t')
		add[0] = indexed_proteins[add[0]]
		add[1] = indexed_proteins[add[1]]
		new_data.append(add)


	return new_data

def ensembler(filename):

	print("Starting Ensemble")
	en_clust = []

	indexed_proteins = starter(filename)

	print(len(indexed_proteins))

	G = nx.Graph()


	for con in indexed_proteins:
		proteinA = con[0]
		proteinB = con[1]
		G.add_node(proteinA)
		G.add_node(proteinB)
		G.add_edge(proteinA, proteinB)
		G.add_edge(proteinA, proteinA)
		G.add_edge(proteinB, proteinB)

	
	while G.number_of_edges()>3000:
		
		print('number of nodes: ' + str(G.number_of_nodes()) )
		print("number of edges: ", str(G.number_of_edges()))
		
		comm = []
		k = 3

		
		while True:
			commk = k_clique_networkX.k_clique_community_finder(k, G)
			comm = comm + commk
			k += 1
			if not commk:
				break
			
		size_k_clique = len(comm)

		#print(nx.adjacency_matrix(G_mcl))



		#B = nx.adjacency_matrix(G)
		#print(B.todense())
		GtoMCL = nx.to_numpy_matrix(G)
		#print("Matrix row len: ", str(len(GtoMCL.getrow(0).todense().tolist()[0])))
		#print("Matrix col len: ", str(len(GtoMCL.getcol(0).todense().tolist())))
		print("MCL")

		comm = comm + mcl.mcl_algorithm(GtoMCL, G.number_of_nodes())
		
		if not comm:
			break

		#print('comm: ')
		#print(comm)
		#print('here4')
		#choosing the best one
		
		max_clust_val = -1
		index = 0
		for c in comm:
			if index < size_k_clique+1:
				sc = metric1(c, G, 'k_clique')
			else: 
				sc = metric1(c, G, 'mcl')
			if sc > max_clust_val:
				max_clust_val = sc
				best_clust = c
				best_clust_index = index
			index += 1


			
		


		#print('elapsed: ** \n')
		#print(end-start)

		en_clust.append(best_clust)
		en_clust_scores.append(max_clust_val)

		if best_clust_index < size_k_clique + 1:
			en_clust_from.append('k_clique')
		else:
			en_clust_from.append('mcl')


		G = graph_updater(G, best_clust)


		filename = 'w.txt'
		nx.write_edgelist(G, filename)


	return(en_clust)


def metric1(clust, G, fr):
	#print(clust)
	counter = 0
	for i in clust:
		for j in clust:
			if have_that_edge(G, i, j) or have_that_edge(G, j, i):
				if i != j:
					counter = counter + 1
				else:
					counter = counter + 0.3


	maxx = len(clust)*(len(clust)-1) + 0.7*len(clust)
	if fr == 'mcl':
		maxx = maxx*100

	return (len(clust)*counter)/maxx



res = ensembler('wi2007.txt')


for i in range(len(res)):
	print(str(i) + 'th CLUSTER: ' + str(res[i]) + 'FROM: ' + str(en_clust_from[i]) +  '\n SCORE: ' + str(en_clust_scores[i]) + ' \n' )
