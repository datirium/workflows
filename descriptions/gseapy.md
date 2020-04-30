GSEAPY: Gene Set Enrichment Analysis in Python
==============================================

Gene Set Enrichment Analysis is a computational method that determines whether an a priori
defined set of genes shows statistically significant, concordant differences between two
biological states (e.g. phenotypes).

GSEA requires as input an expression dataset, which contains expression profiles for multiple samples.
While the software supports multiple input file formats for these datasets, the tab-delimited GCT format
is the most common. The first column of the GCT file contains feature identifiers (gene ids or symbols in
the case of data derived from RNA-Seq experiments). The second column contains a description of the feature;
this column is ignored by GSEA and may be filled with “NA”s. Subsequent columns contain the expression
values for each feature, with one sample's expression value per column. It is important to note that there
are no hard and fast rules regarding how a GCT file's expression values are derived. The important point is
that they are comparable to one another across features within a sample and comparable to one another
across samples. Tools such as DESeq2 can be made to produce properly normalized data (normalized counts)
which are compatible with GSEA.