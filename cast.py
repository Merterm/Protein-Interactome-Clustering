import networkx as nx

def read_graph(filename):
	GG = nx.Graph()
	f = open('dataset.txt', 'r')
	for line in f:
		proteinA = line.split()[0]
		proteinB = line.split()[1]
		GG.add_node(proteinA)
		GG.add_node(proteinB)
		GG.add_edge(proteinA,proteinB)
	f.close()

	G = dict()
	for edge in GG.edges():
		x, y = (edge[0]), (edge[1])
		if x not in G: G[x] = set()
		if y not in G: G[y] = set()
		G[x].add(y)
		G[y].add(x)
	return G


def biggest_degree(G):
	max_deg_node = ''
	max_deg = 0
	for i in G:
		if len(G[i]) > max_deg:
			max_deg_node = i
			max_deg = len(G[i])

	return max_deg_node 


def distancer(clust, node):
	return 1.4

def how_many_edge_it_brings(clust, node):
	#print('hwnsdhwndhwndwh')
	#print(clust)
	counter = 0
	for i in clust:
		for j in G[i]:
	#		print(str(i) + 'hehe boi '+str(j))
			if j == node:
				counter += 1
	return counter

def update_graph(G, clust):
	print(clust)
	for i in clust:
		for j in clust:
			for k in G[i]:
				if k == j:
					G[i][j] = ''

	for i in clust:
		flag = True
		for j in G[i]:
			if j != '':
				flag = False
		if flag == True:
			print('za')
			del G[i]
	return G

def cast(G, n, e):
	all_clust = []
	while(n > 0):
		print('zazaza')
		cur_clust = []
		cur_clust_n = 0 #current cluster's number of nodes
		cur_clust_e = 0 #current cluster's number of edges
		v = biggest_degree(G)
		cur_clust.append(v)
		cur_clust_n += 1

		print(G['Veli'])
		cur_clust1 = cur_clust
		for i in cur_clust1:
			print(i)
			print(1)
			print(len(cur_clust))
			for j in list(G[i]):

				if distancer(i, j) <= 1.5:
					cur_clust.append(j)
					cur_clust_n += 1
					#print(cur_clust_n)
					cur_clust_e += how_many_edge_it_brings(cur_clust, j)

		
		print('qqqqq')
		all_clust.append(cur_clust)
		G = update_graph(G, cur_clust)
		n = len(G)
		print('-------------------------')

	return all_clust

G = read_graph('wi.txt')
print(cast(G, 8, 8))

