GAT: Genomic Association Tester
==============================================

A common question in genomic analysis is whether two sets of genomic intervals overlap significantly.
This question arises, for example, in the interpretation of ChIP-Seq or RNA-Seq data. The Genomic
Association Tester (GAT) is a tool for computing the significance of overlap between multiple sets of
genomic intervals. GAT estimates significance based on simulation.

Gat implemements a sampling algorithm. Given a chromosome (workspace) and segments of interest, for
example from a ChIP-Seq experiment, gat creates randomized version of the segments of interest falling
into the workspace. These sampled segments are then compared to existing genomic annotations.

The sampling method is conceptually simple. Randomized samples of the segments of interest are created
in a two-step procedure. Firstly, a segment size is selected from to same size distribution as the
original segments of interest. Secondly, a random position is assigned to the segment. The sampling stops
when exactly the same number of nucleotides have been sampled. To improve the speed of sampling, segment
overlap is not resolved until the very end of the sampling procedure. Conflicts are then resolved by
randomly removing and re-sampling segments until a covering set has been achieved. Because the size of
randomized segments is derived from the observed segment size distribution of the segments of interest,
the actual segment sizes in the sampled segments are usually not exactly identical to the ones in the
segments of interest. This is in contrast to a sampling method that permutes segment positions within
the workspace.