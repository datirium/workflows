cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement

'sd:metadata':
  - "../metadata/rnaseq-header.cwl"

'sd:upstream':
  kallisto_index: "kallisto-index.cwl"


inputs: 

  kallisto_index:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Kallisto index sample to use for pseudo-alignment:"
    'sd:upstreamSource': "kallisto_index/index_file"
    'sd:localLabel': true
    doc: |
      Kallisto index sample to use for pseudo-alignment, generated from the 'Kallisto index pipeline'.

  input_annotation_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Annotation file (tsv):"
    'sd:upstreamSource': "kallisto_index/input_annotation_file"
    doc: "TSV file containing gene annotations for the reference genome (from kallisto index upstream).\n\n
      Required columns (include headers as row 1 of TSV):\n
      \t1. RefseqId
      \t2. GeneId
      \t3. Chrom (gene/transcript id/name)
      \t4. TxStart
      \t5. TxEnd
      \t6. Strand\n\n
      NOTE: Sequence names (string after the '>') in the transcriptome FASTA must match column 3 (Chrom) of the annotation TSV."

  fastq_file_R1:
    type:
      - File
      - type: array
        items: File
    label: "Read 1 FASTQ file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: |
      Read 1 FASTQ file from a paired-end sequencing run.
    sd:preview:
      position: 5

  fastq_file_R2:
    type:
      - File
      - type: array
        items: File
    label: "Read 2 FASTQ file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: |
      Read 2 FASTQ file that pairs with the input R1 file.
    sd:preview:
      position: 6

  threads:
    type: int?
    default: 10
    label: "Threads:"
    'sd:localLabel': true
    doc: |
      Number of threads to use for steps that support multithreading.
    sd:preview:
      position: 20


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

  rpkm_isoforms:
    type: File
    format: "http://edamontology.org/format_3752"
    label: "kallisto estimated counts per transcript"
    doc: "Quantitation by isoform. NOT ACTUALLY RPKM, output name required for DESeq compatibility, these are kallisto esimate counts per transcript. The na values for unannotated genes have been removed."
    outputSource: kallisto_quant/transcript_counts

  rpkm_genes:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Quantitation by gene name"
    doc: "Quantitation by gene name. NOT ACTUALLY RPKM, output name required for DESeq compatibility, these are derived from kallisto esimate counts per transcript. The na values for unannotated genes have been removed."
    outputSource: group_isoforms/genes_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene Expression'
        Title: 'RPKM, grouped by gene name'

  rpkm_common_tss:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Quantitation by common TSS"
    doc: "Quantitation by common TSS. NOT ACTUALLY RPKM, output name required for DESeq compatibility, these are derived from kallisto esimate counts per transcript. The na values for unannotated genes have been removed."
    outputSource: group_isoforms/common_tss_file

  transcript_counts_all:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "kallisto estimated counts per transcript, with na for unannotated genes"
    doc: "These are kallisto esimate counts per transcript. This file contains na where annotations not available."
    outputSource: kallisto_quant/transcript_counts_standard
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Transcript Counts'

  overview_file:
    type: File
    format: "http://edamontology.org/format_3835"
    label: "summary of inputs"
    doc: "summary of inputs"
    outputSource: kallisto_quant/overview
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  pie_chart_data:
    type: File?
    label: "aligned transcript metrics"
    format: "http://edamontology.org/format_2330"
    doc: "tsv file containing transcript read statistics for sample level pie chart alignment summary"
    outputSource: kallisto_quant/pie_stats
    'sd:preview':
      'sd:visualPlugins':
      - pie:
          colors: ['#b3de69', '#99c0db', '#fdc381', '#fb8072']
          data: [$2, $3, $4, $5]

  kallisto_abundance_tsv:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "raw kallisto count estimates file"
    doc: "raw kallisto count estimates file"
    outputSource: kallisto_quant/kallisto_abundance_file


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
        default: "extracted_R1"
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
        default: "extracted_R2"
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

  kallisto_quant:
    run: ../tools/kallisto-quant-pe.cwl
    in:
      kallisto_index: kallisto_index
      annotation_tsv: input_annotation_file
      fastq_R1: rename_R1/target_file
      fastq_R2: rename_R2/target_file
      threads: threads
    out: [overview, pie_stats, kallisto_abundance_file, kallisto_runinfo_file, transcript_counts, transcript_counts_standard, log_file_stdout, log_file_stderr]

  group_isoforms:  
    run: ../tools/group-isoforms.cwl
    in:
      isoforms_file: kallisto_quant/transcript_counts
    out:
      - genes_file
      - common_tss_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Kallisto transcript quant pipeline paired end"
label: "Kallisto transcript quant pipeline paired end"
s:alternateName: "Kallisto transcript quant pipeline paired end"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/kallisto-quant.cwl
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
  This workflow runs paired end RNA-Seq reads using the kallisto quant tool against a kallisto index reference genome (see "Kallisto index pipeline").
  The kallisto transcript-level quantified samples are then compatible with the DESeq and GSEA downstream workflows.

  ### __Inputs__
   - Kallisto index sample (of experimental organism)
   - R1/R2 FASTQ files of RNA-Seq read data
   - number of threads to use for multithreading processes
  
  ### __Outputs__
   - kallisto quant file (transcript estimate tsv)

  ### __Data Analysis Steps__
  1. cwl calls dockercontainer robertplayer/scidap-kallisto to pseudo align reads using `kallisto quant`.
  2. abundance tsv is formatted, and additional files are produced for gene and common TSS counts for use in differential expression analysis
  3. read and alignment metrics are calculated for the sample piechart, and output to the overview.md file

  ### __References__
    -   Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525-527(2016), doi:10.1038/nbt.3519
