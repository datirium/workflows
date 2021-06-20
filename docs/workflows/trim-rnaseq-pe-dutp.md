Modified original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
**RNA-Seq** basic analysis for a **pair-end** experiment.
A corresponded input [FASTQ](http://maq.sourceforge.net/fastq.shtml) file has to be provided.

Current workflow should be used only with the single-end RNA-Seq data. It performs the following steps:
1. Trim adapters from input FASTQ files
2. Use STAR to align reads from input FASTQ files according to the predefined reference indices; generate unsorted BAM file and alignment statistics file
3. Use fastx_quality_stats to analyze input FASTQ files and generate quality statistics files
4. Use samtools sort to generate coordinate sorted BAM(+BAI) file pair from the unsorted BAM file obtained on the step 1 (after running STAR)
5. Generate BigWig file on the base of sorted BAM file
6. Map input FASTQ files to predefined rRNA reference indices using Bowtie to define the level of rRNA contamination; export resulted statistics to file
7. Calculate isoform expression level for the sorted BAM file and GTF/TAB annotation file using GEEP reads-counting utility; export results to file