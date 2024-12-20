cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


'sd:upstream':
  database: "kraken2-databases.cwl"


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

  k2db:
    type: Directory
    'sd:upstreamSource': "database/k2db"
    label: "Reference genome database for metagenomic sequence classification:"
    'sd:localLabel': true
    doc: "Pre-built kraken2 reference genome database for taxonomic classification of sequencing reads. A 'database' sample needs to be added to your project that will populate this dropdown."
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
    doc: "Read 1 FASTQ file from a paired-end sequencing run"
    sd:preview:
      position: 4

  fastq_file_R2:
    type:
      - File
      - type: array
        items: File
    label: "Read 2 file:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: "Read 2 FASTQ file that pairs with the input R1 file"
    sd:preview:
      position: 5


outputs:

  fastx_statistics_r1:
    type: File
    label: "FASTQ 1 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ 1 quality statistics file"
    outputSource: extract_and_trim/statistics_file_r1
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ 1 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 1 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  fastx_statistics_r2:
    type: File
    label: "FASTQ 2 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ 2 quality statistics file"
    outputSource: extract_and_trim/statistics_file_r2
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ 2 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 2 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  trim_report_r1:
    type: File
    label: "TrimGalore report for read 1"
    doc: "TrimGalore generated log for FASTQ 1"
    outputSource: extract_and_trim/trimmed_fastq_r1_report

  trim_report_r2:
    type: File
    label: "TrimGalore report for read 2"
    doc: "TrimGalore generated log for FASTQ 2"
    outputSource: extract_and_trim/trimmed_fastq_r2_report

  k2_classified_reads_R1:
    type:
      - "null"
      - File
    format: "http://edamontology.org/format_1930"
    label: "classified r1 FASTQ file"
    doc: "classified r1 FASTQ file"
    outputSource: kraken2_classify/k2_classified_R1

  k2_classified_reads_R2:
    type:
      - "null"
      - File
    format: "http://edamontology.org/format_1930"
    label: "classified r2 FASTQ file"
    doc: "classified r2 FASTQ file"
    outputSource: kraken2_classify/k2_classified_R2

  k2_unclassified_reads_R1:
    type:
      - "null"
      - File
    format: "http://edamontology.org/format_1930"
    label: "unclassified r1 FASTQ file"
    doc: "unclassified r1 FASTQ file"
    outputSource: kraken2_classify/k2_unclassified_R1

  k2_unclassified_reads_R2:
    type:
      - "null"
      - File
    format: "http://edamontology.org/format_1930"
    label: "unclassified r2 FASTQ file"
    doc: "unclassified r2 FASTQ file"
    outputSource: kraken2_classify/k2_unclassified_R2

  k2_output:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "kraken2 raw output file"
    doc: "raw per read taxonomic classifications from kraken2"
    outputSource: kraken2_classify/k2_output

  k2_report:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "kraken2 report file"
    doc: "summary of all read taxonomic classifications from kraken2"
    outputSource: kraken2_classify/k2_report

  k2_report_tsv:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "kraken2 report file tsv"
    doc: "summary of all read taxonomic classifications from kraken2 formatted as a tsv"
    outputSource: kraken2_classify/k2_report_tsv
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Kraken Report'
        Title: 'Summary of Taxonomic Classifications'

  k2_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "parsed k2 stderr"
    doc: "markdown parsed standard error captured directly from kraken2 classify command in k2-classify-pe.cwl"
    outputSource: kraken2_classify/k2_stderr
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  krona_plot_link:
    type: File
    format: "http://edamontology.org/format_2331"
    label: "Krona plot - hierarchical visualization of taxonomic classifications"
    doc: "hierarchical visualization of taxonomic classifications"
    outputSource: kraken2_classify/krona_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"


steps:

  extract_and_trim:
    label: "Extract compressed reads if needed, concatentate multi-files per pair, trim each paired read file, then collect qc statistics"
    doc: |
      Extract compressed reads if needed, concatentate multi-files per pair, trim each paired read file, then collect qc statistics
    run: ../tools/extractandtrim-pe.cwl
    in:
      read1file: fastq_file_R1
      read2file: fastq_file_R2
    out:
      - trimmed_fastq_r1
      - trimmed_fastq_r2
      - trimmed_fastq_r1_report
      - trimmed_fastq_r2_report
      - statistics_file_r1
      - statistics_file_r2

  kraken2_classify:
    label: "Kraken2 taxonomic classification of sequence reads"
    doc: |
      Assigns taxonomy to each sequence in the input paired end read files, and reports raw
      classificaiton as well as a summary report.
    run: ../tools/k2-classify-pe.cwl
    in:
      k2db: k2db
      read1file: extract_and_trim/trimmed_fastq_r1
      read2file: extract_and_trim/trimmed_fastq_r2
    out:
      - k2_classified_R1
      - k2_classified_R2
      - k2_unclassified_R1
      - k2_unclassified_R2
      - k2_output
      - k2_report
      - k2_report_tsv
      - k2_stderr
      - krona_html

label: "Kraken2 Metagenomic pipeline paired-end"
doc: |
  This workflow taxonomically classifies paired-end sequencing reads in FASTQ format, that have been 
  adapter trimmed with trimgalore, using Kraken2 with a user-selected pre-built database from a list of
  [genomic index files](https://benlangmead.github.io/aws-indexes/k2).

  ### __Inputs__
  Kraken2 database for taxonomic classification:
    - [Viral (0.5 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz), all refseq viral genomes
    - [MinusB (8.7 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20221209.tar.gz), standard minus bacteria (archaea, viral, plasmid, human1, UniVec_Core)
    - [PlusPFP-16 (15.0 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz), standard (archaea, bacteria, viral, plasmid, human1, UniVec_Core) + (protozoa, fungi & plant) capped at 16 GB (shrunk via random kmer downselect)
    - [EuPathDB46 (34.1 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20201113.tar.gz), eukaryotic pathogen genomes with contaminants removed (https://veupathdb.org/veupathdb/app)
    - [16S_gg_13_5 (73 MB)](https://genome-idx.s3.amazonaws.com/kraken/16S_Greengenes13.5_20200326.tgz), Greengenes 16S rRNA database ([release 13.5](https://greengenes.secondgenome.com/?prefix=downloads/greengenes_database/gg_13_5/), 20200326)\n
    - [16S_silva_138 (112 MB)](https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz), SILVA 16S rRNA database ([release 138.1](https://www.arb-silva.de/documentation/release-1381/), 20200827)
  Read 1 file:
    - FASTA/Q input R1 from a paired end library
  Read 2 file:
    - FASTA/Q input R2 from a paired end library

  ### __Data Analysis Steps__
  1. Trimming the adapters with TrimGalore.
      - This step is particularly important when the reads are long and the fragments are short - resulting in sequencing adapters at the ends of reads. If adapter is not removed the read will not map. TrimGalore can recognize standard adapters, such as Illumina or Nextera/Tn5 adapters.
  2. Generate quality control statistics of trimmed, unmapped sequence data
  3. Classify trimmed reads with kraken2 and selected kraken database

  ### __References__
    - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0