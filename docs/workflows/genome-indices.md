Creates indices for:
* [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
* [bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016)

It performs the following steps:

1. `STAR --runMode genomeGenerate` to generate indices, based on [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) and [GTF](http://mblab.wustl.edu/GTF2.html) input files, returns results as an array of files
2. Outputs indices as [Direcotry](http://www.commonwl.org/v1.0/CommandLineTool.html#Directory) data type
3. Separates *chrNameLength.txt* file from Directory output
4. `bowtie-build` to generate indices requires genome [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) file as input, returns results as a group of main and secondary files