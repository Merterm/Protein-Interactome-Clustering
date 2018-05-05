import numpy
from scipy.sparse import csr_matrix
import csv

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Mert Inan
# Date: 04/05/2018 -
# Description:  This code runs the k-clique algorithm on the protein-protein
#               interaction dataset using the following papers algorithm:
#               "Uncovering the overlapping community structure of complex
#               networks in nature and society"
#               Input:
#               -k
#               -name of the tab-separated file
#               Output:
#               -communities: List of lists
# Version: 1.0
# Changes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# find() equivalent of MATLAB
def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]

# Checks whether a is a member of b
def ismember(a, b):
    bind = {}
    print a, b
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, None) for itm in a]

# Recursive transferring nodes function
def transfer_nodes(A, B, size, all_nodes, cliques):
    #print cliques
    # If A U B or B is already a part of a clique
    found_AB = False
    found_A = False
    #print len(cliques)
    for c in xrange(len(cliques)):
        #print "inside"
        for cc in xrange(len(list(cliques)[c])):
            #print cliques[c]
            if all(ismember(A, list(cliques)[c])):
                found_A = True
            if all(ismember(set.union(A,B), list(cliques)[c])):
                found_AB = True
                break
    #print(len(A))
    #print found_A, found_AB
    if found_AB or ((len(A) != size) and (len(B) != 0)):
        #print("inside found_AB")
        return {}
    elif len(A) == size:
        #print(len(A))
        if found_A:
            print("inside found_A if")
            return {}
        else:
            print("inside found_A else")
            return A
    else:
        #print("inside main else")
        if indices(list(B), lambda x: x >= max(A)) is not []:
            return {}
        C = {}
        max_idx = indices(list(B), lambda x: x >= max(A))
        for w in xrange(max_idx[0],len(B)):
            temp_B = B
            temp_A = A
            temp_A = temp_A.add(list(B)[w])
            con_of_B ={}
            for idx in len(B):
                if all_nodes[list(B)[w]][idx] == 1:
                    con_of_B.add(idx)
            set_Bw = {}
            set_Bw.add(list(B)[w])
            temp_B = set.difference(con_of_B,set_Bw)
            return C.add(transfer_nodes(temp_A,temp_B,size,all_nodes,cliques))


# Initialize a set to store all the unique nodes
all_nodes = set()

## 1. Read TSV
f = open('wi.txt', 'r')
for line in f:
    proteinA = line.split()[0]
    proteinB = line.split()[1]
    all_nodes.add(proteinA)
    all_nodes.add(proteinB)
all_nodes = list(all_nodes)
#print all_nodes
f.close()

## 2. Create an Adjacency Matrix
# Initialize the adjacency matrix
adj = numpy.zeros((len(all_nodes),len(all_nodes)))
#print adj

# Add the values to the adjacency matrix
f = open('wi.txt', 'r')
for line in f:
    proteinA = line.split()[0]
    proteinB = line.split()[1]
    #print proteinA, proteinB
    adj[all_nodes.index(proteinA)][all_nodes.index(proteinB)] = 1
    adj[all_nodes.index(proteinB)][all_nodes.index(proteinA)] = 1
f.close()

# Add 1s to the diagonal of adjacency
for i in xrange(len(all_nodes)):
    for j in xrange(len(all_nodes)):
        if i is j:
            adj[i][j] = 1

#csvfile = open('adj.csv', 'w')
#adjwriter = csv.writer(csvfile)
#for i in xrange(len(all_nodes)):
#    adjwriter.writerow(adj[i])

zerocnt = 0
for i in xrange(len(all_nodes)):
    for j in xrange(len(all_nodes)):
        if adj[i][j] == 0:
            zerocnt = zerocnt + 1

zero_ratio = float(zerocnt)/float((len(all_nodes)*len(all_nodes)))
print zero_ratio*100, "% sparse"

## 3. Find the largest clique size
# Find the degree sequence of adjacency
all_degrees = (numpy.sum(adj,1))-1
degree_sequence = numpy.sort(all_degrees)[::-1]
max_clique_size = 0
for i in xrange(len(degree_sequence)):
    if degree_sequence[i] >= i:
        max_clique_size = i + 1
    else:
        break
print "Max Clique Size: ", max_clique_size

## 4. Find cliques
cliques = set()
for size in xrange(max_clique_size, 3, -1):
    # Modified Bron-Kerbosch Algorithm to find Cliques
    tempAdj = adj
    for n in xrange(len(all_nodes)):
        A = n
        found = numpy.where(tempAdj[n] == 1)[0]
        #print found
        set_A = set()
        set_A.add(A)
        B = set.difference(set(numpy.unique(found)),set_A)
        #print B
        C = transfer_nodes(set_A, B, size, tempAdj, cliques)
        #print C
        if len(C) != 0:
            print C
            cliques.add(list(C)[len(C)-1])
        # Remove processed nodes
        for j in xrange(len(all_nodes)):
            tempAdj[n][j] = 0
            tempAdj[j][n] = 0
    print cliques



## 5. Create a Clique-Clique Adjacency Matrix
## 6. Remove elements less than k from Clique-Clique Matrix
## 7. Print the edges in the communities
