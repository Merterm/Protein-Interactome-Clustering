class Louvain:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges
        self.communities = [n for n in nodes]
        self.actual_partition = []
        self.w = [0 for n in nodes]
        self.edges_of_node = {}
        
        # from the paper, m is the total link weight in the network
        # k_i is the total link weight attached to node i
        self.m = 0
        self.k_values = [0 for n in nodes]
        
        for edge in edges: 
            self.m += edge[1] # add weight of the edge (1) to m
            self.k_values[edge[0][0]] += edge[1] # first node
            self.k_values[edge[0][1]] += edge[1] # second node
            
            # set edges_of_node for first node
            if edge[0][0] not in self.edges_of_node:
                self.edges_of_node[edge[0][0]] = [edge]
            else:
                self.edges_of_node[edge[0][0]].append(edge)
                
            # set edges_of_node for second node
            if edge[0][1] not in self.edges_of_node:
                self.edges_of_node[edge[0][1]] = [edge]
            elif edge[0][0] != edge[0][1]:
                self.edges_of_node[edge[0][1]].append(edge)
                
    def initial_partition(self, network): # network = (nodes, edges)
        partition = [[node] for node in network[0]]
        self.s_in = [0 for node in network[0]]
        self.s_tot = [self.k_values[node] for node in network[0]]
        
        # apply for self-loops
        for edge in network[1]:
            if edge[0][0] == edge[0][1]:
                self.s_in[edge[0][0]] += edge[1]
                self.s_in[edge[0][1]] += edge[1]
                
        return partition
        
    def get_neighbors_of_node(self, node):
        neighbors = []
        
        for edge in self.edges_of_node[node]:
            if edge[0][0] == edge[0][1]:
                continue
            if edge[0][0] == node:
                neighbors.append(edge[0][1])
            if edge[0][1] == node:
                neighbors.append(edge[0][0])
        
        return neighbors
    
    def modularity(self, partition):
        res = 0
        for i in range(len(partition)):
            res += self.s_in[i] / (self.m * 2) - (self.s_tot[i] / (self.m * 2)) ** 2
        return res
        
    def modularity_gain(self, node, com, k_i_in):
        # gain of putting node in com, k_i_in represents the total link weights from node to nodes in comm
        return 2 * k_i_in - self.s_tot[com] * self.k_values[node] / self.m

    def perform_first_phase(self, network): # network = (nodes, edges)
        best_partition = self.initial_partition(network)
        
        while 1:
            improvement = 0
            
            for node in network[0]:
                best_community = self.communities[node]
                best_gain = 0
                
                # remove node
                best_partition[self.communities[node]].remove(node)
                best_shared_links = 0
                
                for edge in self.edges_of_node[node]:
                    if edge[0][0] == edge[0][1]:
                        continue
                    if self.communities[edge[0][0]] == self.communities[node] and edge[0][1] == node or self.communities[edge[0][1]] == self.communities[node] and edge[0][0] == node:
                        best_shared_links += edge[1]
                
                self.s_in[self.communities[node]] -= 2 * (best_shared_links + self.w[node])
                self.s_tot[self.communities[node]] -= self.k_values[node]
                self.communities[node] = -1
                
                communities = {}
                
                for neighbor in self.get_neighbors_of_node(node):
                    community = self.communities[neighbor]
                    if community in communities:
                        continue
                        
                    communities[community] = 1
                    shared_links = 0
                    
                    for edge in self.edges_of_node[node]:
                        if edge[0][0] == edge[0][1]:
                            continue
                        if self.communities[edge[0][0]] == community and edge[0][1] == node or self.communities[edge[0][1]] == community and edge[0][0] == node:
                            shared_links += edge[1]
                            
                    # compute modularity gain by moving node to neighbor
                    gain = self.modularity_gain(node, community, shared_links)
                    if gain > best_gain:
                        best_community = community
                        best_gain = gain
                        best_shared_links = shared_links
                    
                # move node to the community that maximizes modularity gain
                best_partition[best_community].append(node)
                self.communities[node] = best_community
                self.s_in[best_community] += 2 * (best_shared_links + self.w[node])
                self.s_tot[best_community] += self.k_values[node]  
                
                if self.communities[node] != best_community:
                    improvement = 1
            
            if not improvement:
                break
                
        return best_partition

    def perform_second_phase(self, network, partition): # network = (nodes, edges)
        new_nodes = [n for n in range(len(partition))]
        
        # re-index communities and edges
        new_communities = []
        
        a = {}
        i = 0
        for community in self.communities:
            if community in a:
                new_communities.append(a[community])
            else:
                a[community] = i
                new_communities.append(i)
                i += 1
                
        self.communities = new_communities
        
        new_edges = {}
        
        for edge in network[1]:
            ci = self.communities[edge[0][0]]
            cj = self.communities[edge[0][1]]
            try:
                new_edges[(ci, cj)] += edge[1]
            except KeyError:
                new_edges[(ci, cj)] = edge[1]
                
        new_edges = [(k, v) for k, v in new_edges.items()]
        
        # update k_values
        self.edges_of_node = {}
        self.w = [0 for n in new_nodes]
        self.k_values = [0 for n in new_nodes]
        
        for edge in new_edges:
            self.k_values[edge[0][0]] += edge[1]
            self.k_values[edge[0][1]] += edge[1]
            if edge[0][0] == edge[0][1]:
                self.w[edge[0][0]] += edge[1]
            
            if edge[0][0] not in self.edges_of_node:
                self.edges_of_node[edge[0][0]] = [edge]
            else:
                self.edges_of_node[edge[0][0]].append(edge)
            
            if edge[0][1] not in self.edges_of_node:
                self.edges_of_node[edge[0][1]] = [edge]
            elif edge[0][0] != edge[0][1]:
                self.edges_of_node[edge[0][1]].append(edge)
                
        self.communities = [n for n in new_nodes]
        return (new_nodes, new_edges)

    def apply_louvain(self):
        network = (self.nodes, self.edges)
        best_partition = [[node] for node in network[0]]
        best_q = -1
        i = 1
        
        while 1:
            print("==================================== Pass %d ====================================" % i)
            i += 1
            
            partition = self.perform_first_phase(network)
            q = self.modularity(partition)
            partition = [c for c in partition if c]
            print("%s (modularity = %.8f)" % (partition, q))
            
            if self.actual_partition:
                actual = []
                for p in partition:
                    part = []
                    for n in p:
                        part.extend(self.actual_partition[n])
                    actual.append(part)
                self.actual_partition = actual
            else:
                self.actual_partition = partition
            if q == best_q:
                break
            network = self.perform_second_phase(network, partition)
            best_partition = partition
            best_q = q
        
        return (self.actual_partition, best_q)
