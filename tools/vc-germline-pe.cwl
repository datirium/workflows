cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-gatk4:v1.0.0


inputs:

  threads:
    type: int
    inputBinding:
      prefix: "-t"
    doc: |
      Number of threads for parallel processing.

  genome_folder:
    type: Directory
    inputBinding:
      prefix: "-r"
    doc: |
      BWA index directory containing reference genome FASTA with associated indices.

  ploidy:
    type: int
    inputBinding:
      prefix: "-p"
    doc: |
      Ploidy of reference organism (Default 2).

  snpeffdb:
    type: string
    inputBinding:
      prefix: "-s"
    doc: |
      Name of SNPEFF database to use for SNP effect annotation.

  read1file:
    type: File
    inputBinding:
      prefix: "-a"
    doc: |
      Read1 data in FASTQ format, received after paired end sequencing.

  read2file:
    type: File
    inputBinding:
      prefix: "-b"
    doc: |
      Read2 data in FASTQ format, received after paired end sequencing.

  snp_QD:
    type: float
    inputBinding:
      prefix: "-c"
    doc: |
      This is the variant confidence (from the QUAL field) divided by the unfiltered
      depth of non-hom-ref samples. This annotation is intended to normalize the
      variant quality in order to avoid inflation caused when there is deep coverage.
      For filtering purposes it is better to use QD than either QUAL or DP directly.
      Default < 2.0

  snp_FS:
    type: float
    inputBinding:
      prefix: "-d"
    doc: |
      This is the Phred-scaled probability that there is strand bias at the site.
      Strand Bias tells us whether the alternate allele was seen more or less often
      on the forward or reverse strand than the reference allele. When there is
      little to no strand bias at the site, the FS value will be close to 0.
      Default > 60.0

  snp_MQ:
    type: float
    inputBinding:
      prefix: "-e"
    doc: |
      This is the root mean square mapping quality over all the reads at the site.
      Instead of the average mapping quality of the site, this annotation gives the
      square root of the average of the squares of the mapping qualities at the site.
      It is meant to include the standard deviation of the mapping qualities.
      Including the standard deviation allows us to include the variation in the
      dataset. A low standard deviation means the values are all close to the mean,
      whereas a high standard deviation means the values are all far from the mean.
      When the mapping qualities are good at a site, the MQ will be around 60.
      Default < 40.0

  snp_SOR:
    type: float
    inputBinding:
      prefix: "-f"
    doc: |
      This is another way to estimate strand bias using a test similar to the
      symmetric odds ratio test. SOR was created because FS tends to penalize
      variants that occur at the ends of exons. Reads at the ends of exons tend to
      only be covered by reads in one direction and FS gives those variants a bad
      score. SOR will take into account the ratios of reads that cover both alleles.
      Default > 4.0

  snp_MQRankSum:
    type: float
    inputBinding:
      prefix: "-g"
    doc: |
      Compares the mapping qualities of the reads supporting the reference allele
      and the alternate allele.
      Default < -12.5

  snp_ReadPosRankSum:
    type: float
    inputBinding:
      prefix: "-i"
    doc: |
      Compares whether the positions of the reference and alternate alleles are
      different within the reads. Seeing an allele only near the ends of reads is
      indicative of error, because that is where sequencers tend to make the most
      errors.
      Default < -8.0

  indel_QD:
    type: float
    inputBinding:
      prefix: "-j"
    doc: |
      This is the variant confidence (from the QUAL field) divided by the unfiltered
      depth of non-hom-ref samples. This annotation is intended to normalize the
      variant quality in order to avoid inflation caused when there is deep coverage.
      For filtering purposes it is better to use QD than either QUAL or DP directly.
      Default < 2.0

  indel_FS:
    type: float
    inputBinding:
      prefix: "-k"
    doc: |
      This is the Phred-scaled probability that there is strand bias at the site.
      Strand Bias tells us whether the alternate allele was seen more or less often
      on the forward or reverse strand than the reference allele. When there is
      little to no strand bias at the site, the FS value will be close to 0.
      Default > 200.0

  indel_SOR:
    type: float
    inputBinding:
      prefix: "-l"
    doc: |
      This is another way to estimate strand bias using a test similar to the
      symmetric odds ratio test. SOR was created because FS tends to penalize
      variants that occur at the ends of exons. Reads at the ends of exons tend to
      only be covered by reads in one direction and FS gives those variants a bad
      score. SOR will take into account the ratios of reads that cover both alleles.
      Default > 10.0


outputs:

  sorted_dedup_bam:
    type: File
    outputBinding:
      glob: sorted_dedup_reads.bam
    secondaryFiles: sorted_dedup_reads.bam.bai
    doc: |
      sorted deduplicated bam

  chrom_length_tsv:
    type: File
    outputBinding:
      glob: chrom_length.tsv
    doc: |
      Tab delimited chromosome length file: <chromName><TAB><chromSize>

  bqsr2_indels_vcf:
    type: File
    outputBinding:
      glob: bqsr2_indels.vcf.gz
    secondaryFiles: bqsr2_indels.vcf.idx
    doc: |
      indels called after filtering and recalibration

  bqsr2_snps_vcf:
    type: File
    outputBinding:
      glob: bqsr2_snps.vcf.gz
    secondaryFiles: bqsr2_snps.vcf.idx
    doc: |
      snps called after filtering and recalibration

  bqsr2_all_ann_vcf:
    type: File?
    outputBinding:
      glob: bqsr_all.ann.vcf.gz
    secondaryFiles: bqsr_all.ann.vcf.idx
    doc: |
      snps called after filtering and recalibration with effect annotations

  raw_indels_vcf:
    type: File
    outputBinding:
      glob: raw_indels.vcf.gz
    secondaryFiles: raw_indels.vcf.idx
    doc: |
      indels called from gatk HaplotypeCaller using sorted_dedup_reads.bam

  raw_snps_vcf:
    type: File
    outputBinding:
      glob: raw_snps.vcf.gz
    secondaryFiles: raw_snps.vcf.idx
    doc: |
      snps called from gatk HaplotypeCaller using sorted_dedup_reads.bam

  overview:
    type: File
    outputBinding:
      glob: overview.md
    doc: |
      overview of inputs and read, alignment, and variant metrics

  insert_size_histogram:
    type: File
    outputBinding:
      glob: insert_size_histogram.pdf
    doc: |
      histrogram of read insert sizes

  recalibration_plots:
    type: File
    outputBinding:
      glob: recalibration_plots.pdf
    doc: |
      detailed information about recalibration data

  snpEff_summary:
    type: File?
    outputBinding:
      glob: snpEff_summary.html
    doc: |
      summary of SNPEFF annotations

  snpEff_genes:
    type: File?
    outputBinding:
      glob: snpEff_genes.txt
    doc: |
      text file containing gene details, required for link in html summary

  log_file_stdout:
    type: stdout

  log_file_stderr:
    type: stderr


baseCommand: ["run_vc_germlinepe.sh"]
stdout: vc-germline_stdout.log
stderr: vc-germline_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "vc-germline"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/vc-germline.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  Shell wrapper for the Broad Institute's best practices gatk4 germline variant calling pipeline.

  Primary Output files:
  - bqsr2_indels.vcf, filtered and recalibrated indels
  - bqsr2_snps.ann.vcf, filtered and recalibrated snps with effect annotations
  Secondary Output files:
  - raw_indels.vcf, first pass indel calls
  - raw_snps.vcf, first pass snp calls
  Reports:
  - overview.md (input list, alignment metrics, variant counts)
  - insert_size_histogram.pdf
  - recalibration_plots.pdf
  - snpEff_summary.html

  PARAMS:
    SECTION 1: general
    -h  help	show this message
    -t  INT	number of threads
    -r  DIR   path to bwa indices folder of genome reference
    -p  INT    ploidy, number of copies per chromosome (default should be 2)
    -s  DIR    path to SNPEFF database (see SNPEFFDB section below)
    -a  FILE   path to read1 fastq file
    -b  FILE   path to read2 fastq file
    SECTION 2: SNP filters (see Step 6 Notes: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)
    -c  FLOAT  QD
    -d  FLOAT  FS
    -e  FLOAT  MQ
    -f  FLOAT  SOR
    -g  FLOAT  MQRankSum
    -i  FLOAT  ReadPosRankSum
    SECTION 3: Indel filters (see Step 7 Notes: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)
    -j  FLOAT  QD
    -k  FLOAT  FS
    -l  FLOAT  SOR

  SNPEFFDB:
    # snpeff databases
    docker run --rm -ti gatk4-dev /bin/bash
    java -jar $SNPEFF_JAR databases
    # use first column as SNPEFF -v input (e.g. "hg38")
    hg38, Homo_sapiens (USCS), http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg38.zip
    mm10, Mus_musculus, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_mm10.zip
    dm6.03, Drosophila_melanogaster, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_dm6.03.zip
    Rnor_6.0.86, Rattus_norvegicus, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_Rnor_6.0.86.zip
    R64-1-1.86, Saccharomyces_cerevisiae, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_R64-1-1.86.zip


  ____________________________________________________________________________________________________
  References:
    https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
    https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
