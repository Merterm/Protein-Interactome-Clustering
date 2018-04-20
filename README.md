# Protein Interactome Clustering 
Bilkent CS425: Algorithms for Web-Scale Data Project

Community detection algorithms such as Infomap, Louvain modularity and clique-based methods can be used to detect protein groupings and subunits. We propose to implement those algorithms on the C. elegans worm protein interactome data, or viral protein and host protein interactome data. Viral interactome can provide insights about the attachability of the viral proteins to the host proteins, hence, detecting communities of proteins may help tackling the pathway of viral attacks. On the other hand, C. elegans is the worm that has all of its neural networks mapped, hence community detection on its proteins may provide insights into neural functionality.

In this project we are trying to implement an ensemble model combining InfoMap, Louvain Modularity and k-clique methods to cluster overlapping protein interactome data from http://interactome.dfci.harvard.edu/
