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

  fastx_statistics_upstream:
    type: File
    label: "FASTQ 1 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ 1 quality statistics file"
    outputSource: fastx_quality_stats_upstream/statistics_file
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

  fastx_statistics_downstream:
    type: File
    label: "FASTQ 2 statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ 2 quality statistics file"
    outputSource: fastx_quality_stats_downstream/statistics_file
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

  trim_report_upstream:
    type: File
    label: "TrimGalore report Upstream"
    doc: "TrimGalore generated log for FASTQ 1"
    outputSource: trim_fastq/report_file

  trim_report_downstream:
    type: File
    label: "TrimGalore report Downstream"
    doc: "TrimGalore generated log for FASTQ 2"
    outputSource: trim_fastq/report_file_pair

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

  kraken2_classify_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    doc: "captures standard output from k2-classify-pe.cwl"
    outputSource: kraken2_classify/stdout_log

  kraken2_classify_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    doc: "captures standard error from k2-classify-pe.cwl"
    outputSource: kraken2_classify/stderr_log

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

  rename_upstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_1
      target_filename:
        source: extract_fastq_R1/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  rename_downstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_2
      target_filename:
        source: extract_fastq_R2/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  fastx_quality_stats_upstream:
    label: "Quality control of unmapped sequence data for read 1"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_upstream/target_file
    out: [statistics_file]

  fastx_quality_stats_downstream:
    label: "Quality control of unmapped sequence data for read 2"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_downstream/target_file
    out: [statistics_file]

  kraken2_classify:
    label: "Kraken2 taxonomic classification of sequence reads"
    doc: |
      Assigns taxonomy to each sequence in the input paired end read files, and reports raw
      classificaiton as well as a summary report.
    run: ../tools/k2-classify-pe.cwl
    in:
      k2db: k2db
      read1file: rename_upstream/target_file
      read2file: rename_downstream/target_file
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
      - stdout_log
      - stderr_log

$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Kraken2 Metagenomic pipeline paired-end"
label: "Kraken2 Metagenomic pipeline paired-end"
s:alternateName: "Taxonomic Read Classification Workflow with Kraken2 for a paired-end experiment with Trim Galore"

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
  This workflow taxonomically classifies paired-end sequencing reads in FASTQ format, that have been optionally
  adapter trimmed with trimgalore, using Kraken2 and a user-selected pre-built database from a list of
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
  Advanced Inputs Tab (Optional):
     - Number of bases to clip from the 3p end
     - Number of bases to clip from the 5p end

  ### __Outputs__
   - k2db, an upstream database used by kraken2 classifier

  ### __Data Analysis Steps__
  1. Trimming the adapters with TrimGalore.
      - This step is particularly important when the reads are long and the fragments are short - resulting in sequencing adapters at the ends of reads. If adapter is not removed the read will not map. TrimGalore can recognize standard adapters, such as Illumina or Nextera/Tn5 adapters.
  2. Generate quality control statistics of trimmed, unmapped sequence data
  3. (Optional) Clipping of 5' and/or 3' end by the specified number of bases.
  4. Mapping reads to primary genome index with Bowtie.

  ### __References__
    - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0