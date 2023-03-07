cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement


'sd:upstream':
  alignment_index: "bwa-index.cwl"


inputs:

  alias:
    type: string
    label: "Sample short name/Alias:"
    sd:preview:
      position: 1

  condition:
    type: string?
    label: "Experimental condition:"
    sd:preview:
      position: 2

  cells:
    type: string
    label: "Cells:"
    sd:preview:
      position: 3

  catalog:
    type: string?
    label: "Catalog No.:"
    sd:preview:
      position: 4

  index_directory:
    type: Directory
    'sd:upstreamSource': "alignment_index/index_directory"
    label: "BWA index:"
    'sd:localLabel': true
    doc: |
      BWA index directory of the reference genome, to be used for alignment and variant calling.
      An 'index' sample needs to be added to your project that will populate this dropdown.
    sd:preview:
      position: 5

  snpeffdb:
    type:
    - "null"
    - type: enum
      name: "SNPEFF database"
      symbols:
      - hg38
      - mm10
      - Drosophila_melanogaster
      - Rnor_6.0.86
      - R64-1-1.86
    label: "SNPEFF database:"
    'sd:localLabel': true
    doc: |
      SNPEFF database to use for predicting effect of detected SNPs (hg38 (human), mm10 (mouse), Drosophila_melanogaster (fly),
      Rnor_6.0.86 (rat), R64-1-1.86 (yeast), please contact support if you do not see your organism represented in the dropdown).
    sd:preview:
      position: 6

  fastq_file_R1:
    type:
      - File
      - type: array
        items: File
    label: "Read 1 file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: |
      Read 1 FASTQ file from a paired-end sequencing run.
    sd:preview:
      position: 7

  fastq_file_R2:
    type:
      - File
      - type: array
        items: File
    label: "Read 2 file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: |
      Read 2 FASTQ file that pairs with the input R1 file.
    sd:preview:
      position: 8

  ploidy:
    type: int?
    default: 2
    label: "Ploidy:"
    'sd:localLabel': true
    doc: |
      Number of sets of chromosomes in a cell of the reference organism. (Default: 2)
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Threads:"
    'sd:localLabel': true
    doc: |
      Number of threads to use for steps that support multithreading.
    'sd:layout':
      advanced: true

  snp_QD:
    type: float?
    default: 2.01
    label: "SNP - QD <:"
    'sd:localLabel': true
    doc: |
      This is the variant confidence (from the QUAL field) divided by the unfiltered
      depth of non-hom-ref samples. This annotation is intended to normalize the
      variant quality in order to avoid inflation caused when there is deep coverage.
      For filtering purposes it is better to use QD than either QUAL or DP directly.
      Default < 2.01
    'sd:layout':
      advanced: true

  snp_FS:
    type: float?
    default: 59.99
    label: "SNP - FS >:"
    'sd:localLabel': true
    doc: |
      This is the Phred-scaled probability that there is strand bias at the site.
      Strand Bias tells us whether the alternate allele was seen more or less often
      on the forward or reverse strand than the reference allele. When there is
      little to no strand bias at the site, the FS value will be close to 0.
      Default > 59.99
    'sd:layout':
      advanced: true

  snp_MQ:
    type: float?
    default: 40.01
    label: "SNP - MQ <:"
    'sd:localLabel': true
    doc: |
      This is the root mean square mapping quality over all the reads at the site.
      Instead of the average mapping quality of the site, this annotation gives the
      square root of the average of the squares of the mapping qualities at the site.
      It is meant to include the standard deviation of the mapping qualities.
      Including the standard deviation allows us to include the variation in the
      dataset. A low standard deviation means the values are all close to the mean,
      whereas a high standard deviation means the values are all far from the mean.
      When the mapping qualities are good at a site, the MQ will be around 60.
      Default < 40.01
    'sd:layout':
      advanced: true

  snp_SOR:
    type: float?
    default: 3.99
    label: "SNP - SOR >:"
    'sd:localLabel': true
    doc: |
      This is another way to estimate strand bias using a test similar to the
      symmetric odds ratio test. SOR was created because FS tends to penalize
      variants that occur at the ends of exons. Reads at the ends of exons tend to
      only be covered by reads in one direction and FS gives those variants a bad
      score. SOR will take into account the ratios of reads that cover both alleles.
      Default > 3.99
    'sd:layout':
      advanced: true

  snp_MQRankSum:
    type: float?
    default: -12.50
    label: "SNP - MQRankSum <:"
    'sd:localLabel': true
    doc: |
      Compares the mapping qualities of the reads supporting the reference allele
      and the alternate allele.
      Default < -12.50
    'sd:layout':
      advanced: true

  snp_ReadPosRankSum:
    type: float?
    default: -8.01
    label: "SNP - ReadPosRankSum <:"
    'sd:localLabel': true
    doc: |
      Compares whether the positions of the reference and alternate alleles are
      different within the reads. Seeing an allele only near the ends of reads is
      indicative of error, because that is where sequencers tend to make the most
      errors.
      Default < -8.01
    'sd:layout':
      advanced: true

  indel_QD:
    type: float?
    default: 2.01
    label: "Indel - QD <:"
    'sd:localLabel': true
    doc: |
      This is the variant confidence (from the QUAL field) divided by the unfiltered
      depth of non-hom-ref samples. This annotation is intended to normalize the
      variant quality in order to avoid inflation caused when there is deep coverage.
      For filtering purposes it is better to use QD than either QUAL or DP directly.
      Default < 2.01
    'sd:layout':
      advanced: true

  indel_FS:
    type: float?
    default: 199.99
    label: "Indel - FS >:"
    'sd:localLabel': true
    doc: |
      This is the Phred-scaled probability that there is strand bias at the site.
      Strand Bias tells us whether the alternate allele was seen more or less often
      on the forward or reverse strand than the reference allele. When there is
      little to no strand bias at the site, the FS value will be close to 0.
      Default > 199.99
    'sd:layout':
      advanced: true

  indel_SOR:
    type: float?
    default: 9.99
    label: "Indel - SOR >:"
    'sd:localLabel': true
    doc: |
      This is another way to estimate strand bias using a test similar to the
      symmetric odds ratio test. SOR was created because FS tends to penalize
      variants that occur at the ends of exons. Reads at the ends of exons tend to
      only be covered by reads in one direction and FS gives those variants a bad
      score. SOR will take into account the ratios of reads that cover both alleles.
      Default > 9.99
    'sd:layout':
      advanced: true


outputs:

  fastx_statistics_R1:
    type: File
    label: "FASTQ R1 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated quality statistics file for R1"
    outputSource: fastx_quality_stats_R1/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ R1 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ R1 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  fastx_statistics_R2:
    type: File
    label: "FASTQ R2 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated quality statistics file for R2"
    outputSource: fastx_quality_stats_R2/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ R2 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ R2 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  trim_report_R1:
    type: File
    label: "TrimGalore report for R1"
    doc: "TrimGalore generated log for FASTQ R1"
    outputSource: trim_fastq/report_file

  trim_report_R2:
    type: File
    label: "TrimGalore report for R2"
    doc: "TrimGalore generated log for FASTQ R2"
    outputSource: trim_fastq/report_file_pair

  snpEff_summary:
    type: File?
    label: "snpEff summary PDF"
    doc: "snpEff summary PDF"
    outputSource: call_germline_variants/snpEff_summary
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  insert_size_histogram:
    type: File
    label: "insert size histogram PDF"
    doc: "insert size histogram PDF"
    outputSource: call_germline_variants/insert_size_histogram
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  recalibration_plots:
    type: File
    label: "recalibration plots PDF"
    doc: "recalibration plots PDF"
    outputSource: call_germline_variants/recalibration_plots
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  sorted_dedup_bam:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "depulicated sorted alignments bam file"
    doc: "depulicated sorted alignments bam file"
    outputSource: call_germline_variants/sorted_dedup_bam

  chrom_length_tsv:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "tsv containing chromosome headers and length"
    doc: "Tab delimited chromosome length file: <chromName><TAB><chromSize>"
    outputSource: call_germline_variants/chrom_length_tsv

  bqsr2_indels:
    type: File
    format: "http://edamontology.org/format_3016"
    label: "indels called after filtering and recalibration"
    doc: "indels called after filtering and recalibration"
    outputSource: call_germline_variants/bqsr2_indels
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'variant'
        name: "Final Indel calls"
        height: 40
        displayMode: 'COLLAPSED'

  bqsr2_snps:
    type: File
    format: "http://edamontology.org/format_3016"
    label: "snps called after filtering and recalibration"
    doc: "snps called after filtering and recalibration"
    outputSource: call_germline_variants/bqsr2_snps
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'variant'
        name: "Final SNP calls"
        height: 40
        displayMode: 'COLLAPSED'

  bqsr2_snps_ann:
    type: File?
    format: "http://edamontology.org/format_3016"
    label: "snps called after filtering and recalibration with effect annotations"
    doc: "snps called after filtering and recalibration with effect annotations"
    outputSource: call_germline_variants/bqsr2_snps_ann

  raw_indels:
    type: File
    format: "http://edamontology.org/format_3016"
    label: "indels called from gatk HaplotypeCaller using sorted_dedup_reads.bam"
    doc: "indels called from gatk HaplotypeCaller using sorted_dedup_reads.bam"
    outputSource: call_germline_variants/raw_indels

  raw_snps:
    type: File
    format: "http://edamontology.org/format_3016"
    label: "snps called from gatk HaplotypeCaller using sorted_dedup_reads.bam"
    doc: "snps called from gatk HaplotypeCaller using sorted_dedup_reads.bam"
    outputSource: call_germline_variants/raw_snps

  overview:
    type: File
    format: "http://edamontology.org/format_3835"
    label: "variant calling metrics"
    doc: "markdown parsed metrics from run_vc_germlinepe.sh script run by vc-germline-pe.cwl"
    outputSource: call_germline_variants/overview
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  call_germline_variants_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    doc: "captures standard output from vc-germline-pe.cwl"
    outputSource: call_germline_variants/log_file_stdout

  call_germline_variants_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    doc: "captures standard error from vc-germline-pe.cwl"
    outputSource: call_germline_variants/log_file_stderr

  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "BigWig Track"
        height: 120


steps:

  extract_fastq_R1:
    label: "Loading unmapped sequence data for read 1"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_R1
      output_prefix:
        default: "merged_R1"
    out: [fastq_file]

  extract_fastq_R2:
    label: "Loading unmapped sequence data for read 2"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_R2
      output_prefix:
        default: "merged_R2"
    out: [fastq_file]

  trim_fastq:
    label: "Adapter trimming"
    doc: |
      For libraries sequenced on the Illumina platform it’s recommended to remove adapter sequences
      from the reads. If adapters are not trimmed there is a high risk of reads being unmapped to a
      reference genome. This becomes particularly important when the reads are long and the fragments
      are short - resulting in sequencing adapters at the end of read. If adapter trimming will cause
      all the reads become too short (<30bp), this step will be skipped.
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq_R1/fastq_file
      input_file_pair: extract_fastq_R2/fastq_file
      dont_gzip:
        default: true
      length:
        default: 30
      trim1:
        default: false
      paired:
        default: true
    out:
      - trimmed_file
      - trimmed_file_pair
      - report_file
      - report_file_pair

  bypass_trim:
    run: ../tools/bypass-trimgalore-pe.cwl
    in:
      original_fastq_file_1: extract_fastq_R1/fastq_file
      trimmed_fastq_file_1: trim_fastq/trimmed_file
      trimming_report_file_1: trim_fastq/report_file
      original_fastq_file_2: extract_fastq_R2/fastq_file
      trimmed_fastq_file_2: trim_fastq/trimmed_file_pair
      trimming_report_file_2: trim_fastq/report_file_pair
      min_reads_count:
        default: 100  # any small number should be good, as we are catching the case when trimgalore discarded all reads
    out:
      - selected_fastq_file_1
      - selected_report_file_1
      - selected_fastq_file_2
      - selected_report_file_2

  rename_R1:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_1
      target_filename:
        source: extract_fastq_R1/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  rename_R2:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_2
      target_filename:
        source: extract_fastq_R2/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  fastx_quality_stats_R1:
    label: "Quality control of unmapped sequence data for read 1"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_R1/target_file
    out: [statistics_file]

  fastx_quality_stats_R2:
    label: "Quality control of unmapped sequence data for read 2"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_R2/target_file
    out: [statistics_file]

  call_germline_variants:
    label: "Run pipeline for calling germline variants, and read, alignment, and variant stat generation"
    doc: |
      Calls shell wrapper for the Broad Institute's best practices gatk4 germline variant calling pipeline.
    run: ../tools/vc-germline-pe.cwl
    in:
      genome_folder: index_directory
      snpeffdb:
        source: snpeffdb
        valueFrom: $(self)
      ploidy: ploidy
      threads: threads
      read1file: rename_R1/target_file
      read2file: rename_R2/target_file
      snp_QD: snp_QD
      snp_FS: snp_FS
      snp_MQ: snp_MQ
      snp_SOR: snp_SOR
      snp_MQRankSum: snp_MQRankSum
      snp_ReadPosRankSum: snp_ReadPosRankSum
      indel_QD: indel_QD
      indel_FS: indel_FS
      indel_SOR: indel_SOR
    out: [sorted_dedup_bam, chrom_length_tsv, bqsr2_indels, bqsr2_snps, bqsr2_snps_ann, raw_indels, raw_snps, overview, log_file_stdout, log_file_stderr, insert_size_histogram, recalibration_plots, snpEff_summary]

  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: call_germline_variants/sorted_dedup_bam
      chrom_length_file: call_germline_variants/chrom_length_tsv
      pairchip:
        default: true
    out: [bigwig_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Variant calling germline paired-end"
label: "Variant calling germline paired-end"
s:alternateName: "Variant calling pipeline for germline paired-end samples"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/kraken2-classify-pe.cwl
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
  A workflow for the Broad Institute's best practices gatk4 germline variant calling pipeline.

  ## __Outputs__
  #### Primary Output files:
    - bqsr2_indels.vcf, filtered and recalibrated indels (IGV browser)
    - bqsr2_snps.vcf, filtered and recalibrated snps (IGV browser)
    - bqsr2_snps.ann.vcf, filtered and recalibrated snps with effect annotations
  #### Secondary Output files:
    - sorted_dedup_reads.bam, sorted deduplicated alignments (IGV browser)
    - raw_indels.vcf, first pass indel calls
    - raw_snps.vcf, first pass snp calls
  #### Reports:
    - overview.md (input list, alignment metrics, variant counts)
    - insert_size_histogram.pdf
    - recalibration_plots.pdf
    - snpEff_summary.html

  ## __Inputs__
  #### General Info
   - Sample short name/Alias: unique name for sample
   - Experimental condition: condition, variable, etc name (e.g. "control" or "20C 60min")
   - Cells: name of cells used for the sample
   - Catalog No.: vender catalog number if available
   - BWA index: BWA index sample that contains reference genome FASTA with associated indices.
   - SNPEFF database: Name of SNPEFF database to use for SNP effect annotation.
   - Read 1 file: First FASTQ file (generally contains "R1" in the filename)
   - Read 2 file: Paired FASTQ file (generally contains "R2" in the filename)
  #### Advanced
   - Ploidy: number of copies per chromosome (default should be 2)
   - SNP filters: see Step 6 Notes: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
   - Indel filters: see Step 7 Notes: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/

  #### SNPEFF notes:
    Get snpeff databases using `docker run --rm -ti gatk4-dev /bin/bash` then running `java -jar $SNPEFF_JAR databases`.
    Then, use the first column as SNPEFF input (e.g. "hg38").
     - hg38, Homo_sapiens (USCS), http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg38.zip
     - mm10, Mus_musculus, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_mm10.zip
     - dm6.03, Drosophila_melanogaster, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_dm6.03.zip
     - Rnor_6.0.86, Rattus_norvegicus, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_Rnor_6.0.86.zip
     - R64-1-1.86, Saccharomyces_cerevisiae, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_R64-1-1.86.zip

  ### __Data Analysis Steps__
  1. Trimming the adapters with TrimGalore.
     - This step is particularly important when the reads are long and the fragments are short - resulting in sequencing adapters at the ends of reads. If adapter is not removed the read will not map. TrimGalore can recognize standard adapters, such as Illumina or Nextera/Tn5 adapters.
  2. Generate quality control statistics of trimmed, unmapped sequence data
  3. Run germline variant calling pipeline, custom wrapper script implementing Steps 1 - 17 of the Broad Institute's best practices gatk4 germline variant calling pipeline (https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)

  ### __References__
  1. https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
  2. https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
  3. https://software.broadinstitute.org/software/igv/VCF
