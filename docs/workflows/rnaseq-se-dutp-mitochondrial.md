Slightly changed original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
**RNA-Seq** basic analysis for **strand specific single-read** experiment.
An additional steps were added to map data to mitochondrial chromosome only and then merge the output.

Experiment files in [FASTQ](http://maq.sourceforge.net/fastq.shtml) format either compressed or not can be used.

Current workflow should be used only with single-read strand specific RNA-Seq data. It performs the following steps:
1. `STAR` to align reads from input FASTQ file according to the predefined reference indices; generate unsorted BAM file and alignment statistics file
2. `fastx_quality_stats` to analyze input FASTQ file and generate quality statistics file
3. `samtools sort` to generate coordinate sorted BAM(+BAI) file pair from the unsorted BAM file obtained on the step 1 (after running STAR)
5. Generate BigWig file on the base of sorted BAM file
6. Map input FASTQ file to predefined rRNA reference indices using Bowtie to define the level of rRNA contamination; export resulted statistics to file
7. Calculate isoform expression level for the sorted BAM file and GTF/TAB annotation file using `GEEP` reads-counting utility; export results to file