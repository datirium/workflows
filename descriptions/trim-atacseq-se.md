This ATAC pipeline is based on original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
**ChIP-Seq** basic analysis workflow for a **single-read** experiment with Trim Galore. The pipeline was adapted for
ATAC-Seq single-read data analysis by updating genome coverage step.
### Data Analysis Steps
SciDAP starts from the .fastq files which most DNA cores and commercial NGS companies return. Starting from raw data allows us to ensure that all experiments have been processed in the same way and simplifies the deposition of data to GEO upon publication. The data can be uploaded from users computer, downloaded directly from an ftp server of the core facility by providing a URL or from GEO by providing SRA accession number.
Our current pipelines include the following steps:
1. Trimming the adapters with TrimGalore. This step is particularly important when the reads are long and the fragments are short as in ATAC -resulting in sequencing adapters at the end of read. If adapter is not removed the read will not map. TrimGalore can recognize standard adapters, such as Nexterra/Tn5 adapters.
2. QC
3. (Optional) trimming adapters on 5' or 3' end by the specified number of bases.
4. Mapping reads with BowTie. Only uniquely mapped reads with less than 3 mismatches are used in the downstream analysis. Results are saved as a .bam file.
5. Reads mapping to chromosome M are removed. Since there are many copies of chromosome M in the cell and it is not protected by histones, some ATAC libraries have up to 50% of reads mapping to chrM. We recommend using OMNI-ATAC protocol that reduces chrM reads and provides better specificity. 
6.  (Optional) Removal of duplicates (reads/pairs of reads mapping to exactly same location). This step is used to remove reads overamplified in PCR. Unfortunately, it may also remove "good" reads. We usually do not remove duplicates unless the library is heavily duplicated. Please note that MACS2 will remove 'excessive' duplicates during peak calling ina smart way (those not supported by other nearby reads).
7.  Peakcalling by MACS2. (Optionally), it is possible to specify read extension length for MACS2 to use if the length determined automatically is wrong. 
8.  Generation of BigWig coverage files for display on the browser. Since the cuts by the Tn5 transposome are 9bp apart, we show coverage by 9bp reads rather than fragments as in ChIP-Seq. The coverage shows the number of fragments at each base in the genome normalized to the number of millions of mapped reads. This way the peak of coverage will be located at the most accessible site. 

### Details
_Trim Galore_ is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt)
and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to consistently
apply adapter and quality trimming to FastQ files, with extra functionality for RRBS data.

In outputs it returns coordinate sorted BAM file alongside with index BAI file, quality
statistics of the input FASTQ file, reads coverage in a form of BigWig file, peaks calling
data in a form of narrowPeak or broadPeak files, islands with the assigned nearest genes and
region type, data for average tag density plot (on the base of BAM file).

Workflow starts with step *fastx\_quality\_stats* from FASTX-Toolkit
to calculate quality statistics for input FASTQ file.

At the same time `bowtie` is used to align
reads from input FASTQ file to reference genome *bowtie\_aligner*. The output of this step
is unsorted SAM file which is being sorted and indexed by `samtools sort` and `samtools index`
*samtools\_sort\_index*.

Based on workflow’s input parameters indexed and sorted BAM file
can be processed by `samtools rmdup` *samtools\_rmdup* to get rid of duplicated reads.
If removing duplicates is not required the original input BAM and BAI
files return. Otherwise step *samtools\_sort\_index\_after\_rmdup* repeat `samtools sort` and `samtools index` with BAM and BAI files.

Right after that `macs2 callpeak` performs peak calling *macs2\_callpeak*. On the base of returned outputs the next step
*macs2\_island\_count* calculates the number of islands and estimated fragment size. If the last
one is less that 80bp (hardcoded in the workflow) `macs2 callpeak` is rerun again with forced fixed
fragment size value (*macs2\_callpeak\_forced*). If at the very beginning it was set in workflow
input parameters to force run peak calling with fixed fragment size, this step is skipped and the
original peak calling results are saved.

In the next step workflow again calculates the number of islands and estimates fragment size (*macs2\_island\_count\_forced*)
for the data obtained from *macs2\_callpeak\_forced* step. If the last one was skipped the results from *macs2\_island\_count\_forced* step
are equal to the ones obtained from *macs2\_island\_count* step.

Next step (*macs2\_stat*) is used to define which of the islands and estimated fragment size should be used
in workflow output: either from *macs2\_island\_count* step or from *macs2\_island\_count\_forced* step. If input
trigger of this step is set to True it means that *macs2\_callpeak\_forced* step was run and it returned different
from *macs2\_callpeak* step results, so *macs2\_stat* step should return [fragments\_new, fragments\_old, islands\_new],
if trigger is False the step returns [fragments\_old, fragments\_old, islands\_old], where sufix "old" defines
results obtained from *macs2\_island\_count* step and sufix "new" - from *macs2\_island\_count\_forced* step.

The following two steps (*bamtools\_stats* and *bam\_to\_bigwig*) are used to calculate coverage on the base
of input BAM file and save it in BigWig format. For that purpose bamtools stats returns the number of
mapped reads number which is then used as scaling factor by bedtools genomecov when it performs coverage
calculation and saves it in BED format. The last one is then being sorted and converted to BigWig format by
bedGraphToBigWig tool from UCSC utilities. To adapt the pipeline for ATAC-Seq data analysis we calculate genome
coverage using only the first 9 bp from every read.

Step *get\_stat* is used to return a text file with statistics
in a form of [TOTAL, ALIGNED, SUPRESSED, USED] reads count.

Step *island\_intersect* assigns genes and regions to the islands obtained from *macs2\_callpeak\_forced*.
Step *average\_tag\_density* is used to calculate data for average tag density plot on the base of BAM file.