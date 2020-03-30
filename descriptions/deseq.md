Differential gene expression analysis
=====================================

Differential gene expression analysis based on the negative binomial distribution

Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.

DESeq1
------

High-throughput sequencing assays such as RNA-Seq, ChIP-Seq or barcode counting provide quantitative readouts
in the form of count data. To infer differential signal in such data correctly and with good statistical power,
estimation of data variability throughout the dynamic range and a suitable error model are required.
Simon Anders and Wolfgang Huber propose a method based on the negative binomial distribution, with variance and mean
linked by local regression and present an implementation, [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html),
as an R/Bioconductor package

DESeq2
------

In comparative high-throughput sequencing assays, a fundamental task is the analysis of count data,
such as read counts per gene in RNA-seq, for evidence of systematic changes across experimental conditions.
Small replicate numbers, discreteness, large dynamic range and the presence of outliers require a
suitable statistical approach. [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html),
a method for differential analysis of count data,
using shrinkage estimation for dispersions and fold changes to improve stability and interpretability of estimates.
This enables a more quantitative analysis focused on the strength rather than the mere presence of differential expression.