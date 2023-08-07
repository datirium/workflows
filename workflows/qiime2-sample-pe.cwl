cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement


inputs:

  alias:
    type: string
    label: "Sample short name/Alias:"
    doc: |
      Used for samplename in downstream analyses. Ensure this is the same name used in the metadata samplesheet.
    sd:preview:
      position: 1

  environment:
    type: string?
    label: "Environment:"
    doc: |
      Optional input.
    sd:preview:
      position: 2

  catalog:
    type: string?
    label: "Catalog No.:"
    doc: |
      Optional input.
    sd:preview:
      position: 3

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
      position: 11

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
      position: 12

  trimLeftF:
    type: float?
    default: 0
    label: "Trim 5' of R1:"
    'sd:localLabel': true
    doc: |
      Recommended if adapters are still on the input sequences. Trims the first J bases from the 5' end of each forward read.

  trimLeftR:
    type: float?
    default: 0
    label: "Trim 5' of R2:"
    'sd:localLabel': true
    doc: |
      Recommended if adapters are still on the input sequences. Trims the first K bases from the 5' end of each reverse read.

  truncLenF:
    type: float?
    default: 1000
    label: "Truncate 3' of R1:"
    'sd:localLabel': true
    doc: |
      Recommended if quality drops off along the length of the read. Clips the remaining bases starting a INT from the 5' end of the forward read.

  truncLenR:
    type: float?
    default: 1000
    label: "Truncate 3' of R2:"
    'sd:localLabel': true
    doc: |
      Recommended if quality drops off along the length of the read. Clips the remaining bases starting a INT from the 5' end of the reverse read.

  threads:
    type: int?
    default: 4
    label: "Threads:"
    'sd:localLabel': true
    doc: |
      Number of threads to use for steps that support multithreading.


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

  overview:
    type: File
    format: "http://edamontology.org/format_3835"
    label: "summary of inputs"
    doc: "summary of inputs"
    outputSource: qiime_pipeline/fastq_summary
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  fastq_summary:
    type: File?
    label: "summary of input read data"
    doc: "summary of input read data"
    outputSource: qiime_pipeline/alpha_rarefaction
    'sd:visualPlugins':
    - qiime2:
        tab: 'Overview'
        target: "_blank"

  alpha_rarefaction:
    type: File?
    doc: "plot of OTU rarefaction"
    outputSource: qiime_pipeline/alpha_rarefaction
    'sd:visualPlugins':
    - qiime2:
        tab: 'Overview'
        target: "_blank"

  taxa_bar_plots:
    type: File?
    doc: "bar plot for exploring the taxonomic composition of the sample"
    outputSource: qiime_pipeline/taxa_bar_plots
    'sd:visualPlugins':
    - qiime2:
        tab: 'Overview'
        target: "_blank"


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
      are short - resulting in sequencing adapters at the end of a read. If adapter trimming will cause
      all the reads to become too short (<30bp), this step will be skipped.
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
    out: [target_file]

  rename_R2:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_2
      target_filename:
        source: extract_fastq_R2/fastq_file
        valueFrom: $(self.basename)
    out: [target_file]

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

  qiime_pipeline:
    label: "Run pipeline for processing a single 16S metagenomic sample using qiime2"
    doc: |
      Calls shell wrapper for QIIME2's 16S metagenomic processing pipeline.
    run: ../tools/qiime2-sample-pe.cwl
    in:
      samplename: alias
      read1file: rename_R1/target_file
      read2file: rename_R2/target_file
      trimLeftF: trimLeftF
      trimLeftR: trimLeftR
      truncLenF: truncLenF
      trimLeftR: trimLeftR
      threads: threads
    out:
      - overview
      - fastq_summary
      - alpha_rarefaction
      - taxa_bar_plots
      - log_file_stdout
      - log_file_stderr


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "16S metagenomic paired-end with QIIME2"
label: "16S metagenomic paired-end with QIIME2"
s:alternateName: "16S metagenomic paired-end pipeline using QIIME2 for single sample analysis"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/qiime2-sample-pe.cwl
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
