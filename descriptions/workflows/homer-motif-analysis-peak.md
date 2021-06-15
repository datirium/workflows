Motif Finding with HOMER with target and background regions from peaks
---------------------------------------------------

HOMER contains a novel motif discovery algorithm that was designed for regulatory element analysis
in genomics applications (DNA only, no protein). It is a differential motif discovery algorithm,
which means that it takes two sets of sequences and tries to identify the regulatory elements that
are specifically enriched in on set relative to the other. It uses ZOOPS scoring (zero or one
occurrence per sequence) coupled with the hypergeometric enrichment calculations (or binomial) to
determine motif enrichment. HOMER also tries its best to account for sequenced bias in the dataset.
It was designed with ChIP-Seq and promoter analysis in mind, but can be applied to pretty much any
nucleic acids motif finding problem.

For more information please refer to:
-------------------------------------
[Official documentation](http://homer.ucsd.edu/homer/motif/)