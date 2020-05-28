Motif Finding with HOMER with random background regions
---------------------------------------------------

HOMER contains a novel motif discovery algorithm that was designed for regulatory element analysis
in genomics applications (DNA only, no protein). It is a differential motif discovery algorithm,
which means that it takes two sets of sequences and tries to identify the regulatory elements that
are specifically enriched in on set relative to the other. It uses ZOOPS scoring (zero or one
occurrence per sequence) coupled with the hypergeometric enrichment calculations (or binomial) to
determine motif enrichment. HOMER also tries its best to account for sequenced bias in the dataset.
It was designed with ChIP-Seq and promoter analysis in mind, but can be applied to pretty much any
nucleic acids motif finding problem.

Here is how we generate background for Motifs Analysis
-------------------------------------
1. Take input file with regions in a form of “chr" “start" “end"
2. Sort and remove duplicates from this regions file
3. Extend each region in 20Kb into both directions
4. Merge all overlapped extended regions
5. Subtract not extended regions from the extended ones
6. Randomly distribute not extended regions within the regions
    that we got as a result of the previous step
7. Get fasta file from these randomly distributed regions (from the previous step). Use it as background
   
For more information please refer to:
-------------------------------------
[Official documentation](http://homer.ucsd.edu/homer/motif/)