import numpy
from scipy.sparse.csgraph import connected_components
import sys
import networkx as nx
import collections
from networkx.algorithms.approximation.clique import max_clique


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Mert Inan
# Date: 04/05/2018 -
# Description:  This code runs the k-clique algorithm on the protein-protein
#               interaction dataset using the following papers algorithm:
#               "Uncovering the overlapping community structure of complex
#               networks in nature and society"
#               Input:
#               -k: size of clique
#               -filename: name of the tab-separated file
#               Output:
#               -communities: List of lists
# Version: 1.0
# Changes:
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def k_clique_community_finder(k,filename):

    # Initialize a graph structure
    G = nx.Graph()

    ## 1. Read TSV
    f = open(filename, 'r')
    for line in f:
        proteinA = line.split()[0]
        proteinB = line.split()[1]
        G.add_node(proteinA)
        G.add_node(proteinB)
        G.add_edge(proteinA,proteinB)
    f.close()

    #print(G.number_of_nodes())
    #print(G.number_of_edges())

    ## 4. Find cliques
    cliques = nx.find_cliques(G)
    cliques = [frozenset(c) for c in cliques if len(c) >= k]

    # Keep a dictionary of proteins and their cliques
    protein_dict = collections.defaultdict(list)
    for clique in cliques:
        for node in clique:
            protein_dict[node].append(clique)

    ## 5. Create a Clique-Clique Adjacency Matrix
    clique_clique = numpy.zeros((len(cliques), len(cliques)))
    # Fill the matrix with the number of common nodes between cliques
    for i in xrange(len(cliques)):
        clique1 = list(cliques)[i]
        #Just look at the adjacent cliques
        adjacent_cliques = set()
        for n in clique1:
            for adj_clique in protein_dict[n]:
                if clique1 != adj_clique:
                    adjacent_cliques.add(adj_clique)
        for clique2 in adjacent_cliques:
            j = list(cliques).index(clique2)
            if i == j:
                if len(set(clique1).intersection(set(clique2))) >= k:
                    clique_clique[i][j] = 1
            else:
                if len(set(clique1).intersection(set(clique2))) >= (k-1):
                    clique_clique[i][j] = 1

    ## 7. Print the communities
    n_comm, comm = connected_components(clique_clique)
    #print comm
    communities = []
    for comm_idx in xrange(n_comm-1):
        #print str(comm_idx)
        comm_set = set()
        for i in xrange(len(cliques)):
            if comm[i] == comm_idx:
                frozen_set =  list(cliques)[i]
                for elm in frozen_set:
                    comm_set.add(elm)
        #print comm_set
        communities.append(list(comm_set))

    #print(communities)
    return communities

if __name__ == "__main__":
    k = 3
    while True:
        comm = k_clique_community_finder(k,sys.argv[1])
        if not comm:
            break
        print("Running k-clique for k = " + str(k))
        print(comm)
        print("\n")
        k = k + 1
