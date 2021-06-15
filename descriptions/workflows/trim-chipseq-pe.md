This ChIP-Seq pipeline is based on the  original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
**ChIP-Seq** basic analysis workflow for a **paired-end** experiment with Trim Galore.
### Data Analysis
SciDAP starts from the .fastq files which most DNA cores and commercial NGS companies return. Starting from raw data allows us to ensure that all experiments have been processed in the same way and simplifies the deposition of data to GEO upon publication. The data can be uploaded from users computer, downloaded directly from an ftp server of the core facility by providing a URL or from GEO by providing SRA accession number.
Our current pipelines include the following steps:
1. Trimming the adapters with TrimGalore. This step is particularly important when the reads are long and the fragments are short-resulting in sequencing adapters at the end of read. If adapter is not removed the read will not map. TrimGalore can recognize standard adapters, such as Illumina or Nexterra/Tn5 adapters.
2. QC
3. (Optional) trimming adapters on 5' or 3' end by the specified number of bases.
4. Mapping reads with BowTie. Only uniquely mapped reads with less than 3 mismatches are used in the downstream analysis. Results are saved as a .bam file.
5.  (Optional) Removal of duplicates (reads/pairs of reads mapping to exactly same location). This step is used to remove reads overamplified in PCR. Unfortunately, it may also remove "good" reads. We usually do not remove duplicates unless the library is heavily duplicated. Please note that MACS2 will remove 'excessive' duplicates during peak calling ina smart way (those not supported by other nearby reads).
6.  Peakcalling by MACS2. (Optionally), it is possible to specify read extension length for MACS2 to use if the length determined automatically is wrong. 
7.  Generation of BigWig coverage files for display on the browser. The coverage shows the number of fragments at each base in the genome normalized to the number of millions of mapped reads. In the case of PE sequencing the fragments are real, but in the case of single reads the fragments are estimated by extending reads to the average fragment length found by MACS2 or specified by the user in 6.

### Details
_Trim Galore_ is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt)
and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to consistently
apply adapter and quality trimming to FastQ files, with extra functionality for RRBS data.

A [FASTQ](http://maq.sourceforge.net/fastq.shtml) input file has to be provided.


In outputs it returns coordinate sorted BAM file alongside with index BAI file,
quality statistics for both the input FASTQ files, reads coverage in a form of BigWig file,
peaks calling data in a form of narrowPeak or broadPeak files, islands with the assigned nearest
genes and region type, data for average tag density plot (on the base of BAM file).

Workflow starts with running fastx_quality_stats (steps fastx_quality_stats_upstream and
fastx_quality_stats_downstream) from FASTX-Toolkit to calculate quality statistics for both upstream
and downstream input FASTQ files. At the same time Bowtie is used to align reads from input FASTQ
files to reference genome (Step bowtie_aligner). The output of this step is unsorted SAM file which
is being sorted and indexed by samtools sort and samtools index (Step samtools_sort_index).
Depending on workflow’s input parameters indexed and sorted BAM file
could be processed by samtools rmdup (Step samtools_rmdup) to remove all possible read duplicates.
In a case when removing duplicates is not necessary the step returns original input BAM and BAI
files without any processing. If the duplicates were removed the following step
(Step samtools_sort_index_after_rmdup) reruns samtools sort and samtools index with BAM and BAI files,
if not - the step returns original unchanged input files. Right after that macs2 callpeak performs
peak calling (Step macs2_callpeak). On the base of returned outputs the next step
(Step macs2_island_count) calculates the number of islands and estimated fragment size. If the last
one is less that 80 (hardcoded in a workflow) macs2 callpeak is rerun again with forced fixed
fragment size value (Step macs2_callpeak_forced). If at the very beginning it was set in workflow
input parameters to force run peak calling with fixed fragment size, this step is skipped and the
original peak calling results are saved. In the next step workflow again calculates the number
of islands and estimated fragment size (Step macs2_island_count_forced) for the data obtained from
macs2_callpeak_forced step. If the last one was skipped the results from macs2_island_count_forced step
are equal to the ones obtained from macs2_island_count step.
Next step (Step macs2_stat) is used to define which of the islands and estimated fragment size should be used
in workflow output: either from macs2_island_count step or from macs2_island_count_forced step. If input
trigger of this step is set to True it means that macs2_callpeak_forced step was run and it returned different
from macs2_callpeak step results, so macs2_stat step should return [fragments_new, fragments_old, islands_new],
if trigger is False the step returns [fragments_old, fragments_old, islands_old], where sufix "old" defines
results obtained from macs2_island_count step and sufix "new" - from macs2_island_count_forced step.
The following two steps (Step bamtools_stats and bam_to_bigwig) are used to calculate coverage on the base
of input BAM file and save it in BigWig format. For that purpose bamtools stats returns the number of
mapped reads number which is then used as scaling factor by bedtools genomecov when it performs coverage
calculation and saves it in BED format. The last one is then being sorted and converted to BigWig format by
bedGraphToBigWig tool from UCSC utilities. Step get_stat is used to return a text file with statistics
in a form of [TOTAL, ALIGNED, SUPRESSED, USED] reads count. Step island_intersect assigns genes and regions
to the islands obtained from macs2_callpeak_forced. Step average_tag_density is used to calculate data for
average tag density plot on the base of BAM file.
