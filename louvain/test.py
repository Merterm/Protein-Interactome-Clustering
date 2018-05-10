from louvain import Louvain

def main():
    nodes = {} # dictionary of integers
    edges = [] # list of ((node, node), weight)
    
    def sort_graph(nodes, edges):
        # sort based on node ids
        nodes = list(nodes.keys())
        nodes.sort()
    
        sorted_nodes = []
        sorted_edges = []
    
        i = 0
        a = {}
    
        for node in nodes:
            sorted_nodes.append(i)
            a[node] = i
            i += 1
    
        for edge in edges:
            sorted_edges.append(((a[edge[0][0]], a[edge[0][1]]), edge[1]))
        
        return (sorted_nodes, sorted_edges)
    
    f = open("data/wi2007.txt", 'r')
    lines = list(f)
    
    for line in lines:
        # each line in our dataset is in the format below
        # first_node second_node
        node_list = line.split()
        nodes[node_list[0]] = 1
        nodes[node_list[1]] = 1
        weight = 1
        
        edges.append(((node_list[0], node_list[1]), weight))
        
    sorted_nodes, sorted_edges = sort_graph(nodes, edges)
    print("Nodes: %d, Edges: %d" % (len(sorted_nodes), len(sorted_edges)))
    
    partition, modularity = Louvain(sorted_nodes, sorted_edges).apply_louvain()

if __name__ == "__main__":
    main()
    
