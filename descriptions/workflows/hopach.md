Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH)
===============================================================

The HOPACH clustering algorithm builds a hierarchical tree of clusters by recursively partitioning a data set, while
ordering and possibly collapsing clusters at each level. The algorithm uses the Mean/Median Split Silhouette (MSS) criteria
to identify the level of the tree with maximally homogeneous clusters. It also runs the tree down to produce a final
ordered list of the elements. The non-parametric bootstrap allows one to estimate the probability that each element
belongs to each cluster (fuzzy clustering).