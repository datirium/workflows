The original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
**ChIP-Seq** basic analysis workflow for a **paired-end** experiment.
A [FASTQ](http://maq.sourceforge.net/fastq.shtml) input file has to be provided.

The pipeline produces a sorted BAM file alongside with index BAI file, quality
statistics of the input FASTQ file, coverage by estimated fragments as a BigWig file, peaks calling
data in a form of narrowPeak or broadPeak files, islands with the assigned nearest genes and
region type, data for average tag density plot.

Workflow starts with step *fastx\_quality\_stats* from FASTX-Toolkit
to calculate quality statistics for input FASTQ file.

At the same time `bowtie` is used to align
reads from input FASTQ file to reference genome *bowtie\_aligner*. The output of this step
is an unsorted SAM file which is being sorted and indexed by `samtools sort` and `samtools index`
*samtools\_sort\_index*.

Depending on workflow’s input parameters indexed and sorted BAM file
can be processed by `samtools rmdup` *samtools\_rmdup* to get rid of duplicated reads.
If removing duplicates is not required the original BAM and BAI
files are returned. Otherwise step *samtools\_sort\_index\_after\_rmdup* repeat `samtools sort` and `samtools index` with BAM and BAI files without duplicates.

Next `macs2 callpeak` performs peak calling *macs2\_callpeak* and the next step
reports *macs2\_island\_count*  the number of islands and estimated fragment size. If the latter
is less that 80bp (hardcoded in the workflow) `macs2 callpeak` is rerun again with forced fixed
fragment size value (*macs2\_callpeak\_forced*). It is also possible to force MACS2 to use pre set  fragment size in the first place.

Next step (*macs2\_stat*) is used to define which of the islands and estimated fragment size should be used
in workflow output: either from *macs2\_island\_count* step or from *macs2\_island\_count\_forced* step. If input
trigger of this step is set to True it means that *macs2\_callpeak\_forced* step was run and it returned different
from *macs2\_callpeak* step results, so *macs2\_stat* step should return [fragments\_new, fragments\_old, islands\_new],
if trigger is False the step returns [fragments\_old, fragments\_old, islands\_old], where sufix "old" defines
results obtained from *macs2\_island\_count* step and sufix "new" - from *macs2\_island\_count\_forced* step.

The following two steps (*bamtools\_stats* and *bam\_to\_bigwig*) are used to calculate coverage from BAM file and save it in BigWig format. For that purpose bamtools stats returns the number of
mapped reads  which is then used as scaling factor by bedtools genomecov when it performs coverage
calculation and saves it as a BEDgraph file whichis then  sorted and converted to BigWig format by
bedGraphToBigWig tool from UCSC utilities. Step *get\_stat* is used to return a text file with statistics
in a form of [TOTAL, ALIGNED, SUPRESSED, USED] reads count.

Step *island\_intersect* assigns nearest genes and regions to the islands obtained from *macs2\_callpeak\_forced*.
Step *average\_tag\_density* is used to calculate data for average tag density plot from the BAM file.