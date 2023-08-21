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
      Where the sample was collected. Optional input.
    sd:preview:
      position: 2

  catalog:
    type: string?
    label: "Catalog No.:"
    doc: |
      If available. Optional input.
    sd:preview:
      position: 3

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
      position: 11

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
      position: 12

  trimLeftF:
    type: int?
    default: 0
    label: "Trim 5' of R1:"
    'sd:localLabel': true
    doc: |
      Recommended if adapters are still on the input sequences. Trims the first J bases from the 5' end of each forward read.

  trimLeftR:
    type: int?
    default: 0
    label: "Trim 5' of R2:"
    'sd:localLabel': true
    doc: |
      Recommended if adapters are still on the input sequences. Trims the first K bases from the 5' end of each reverse read.

  truncLenF:
    type: int
    default: 0
    label: "Truncate 3' of R1:"
    'sd:localLabel': true
    doc: |
      Clips the forward read starting M bases from the 5' end (before trimming). If base quality is OK for entire read, value should be set to the expected number of Illumina cycles for R1.

  truncLenR:
    type: int
    default: 0
    label: "Truncate 3' of R2:"
    'sd:localLabel': true
    doc: |
      Clips the reverse read starting N bases from the 5' end (before trimming).  If base quality is OK for entire read, value should be set to the expected number of Illumina cycles for R2.

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
    outputSource: qiime_pipeline/overview
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  fastq_summary:
    type: File?
    label: "Summary of input FASTQ reads"
    doc: "summary of input read data"
    outputSource: qiime_pipeline/fastq_summary
    'sd:visualPlugins':
    - qiime2:
        tab: 'Overview'
        target: "_blank"

  alpha_rarefaction:
    type: File?
    label: "Alpha rarefaction curve"
    doc: "plot of OTU rarefaction"
    outputSource: qiime_pipeline/alpha_rarefaction
    'sd:visualPlugins':
    - qiime2:
        tab: 'Overview'
        target: "_blank"

  taxa_bar_plots:
    type: File?
    label: "Taxonomic classifications bar plot"
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
      truncLenR: truncLenR
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

s:name: "16S metagenomic paired-end QIIME2 Sample"
label: "16S metagenomic paired-end QIIME2 Sample"
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
  A workflow for processing a single 16S sample via a QIIME2 pipeline.

  ## __Outputs__
  #### Output files:
    - overview.md, list of inputs
    - demux.qzv, summary visualizations of imported data
    - alpha-rarefaction.qzv, plot of OTU rarefaction
    - taxa-bar-plots.qzv, relative frequency of taxomonies barplot

  ## __Inputs__
  #### General Info
   - Sample short name/Alias: Used for samplename in downstream analyses. Ensure this is the same name used in the metadata samplesheet.
   - Environment: where the sample was collected
   - Catalog No.: catalog number if available (optional)
   - Read 1 FASTQ file: Read 1 FASTQ file from a paired-end sequencing run.
   - Read 2 FASTQ file: Read 2 FASTQ file that pairs with the input R1 file.
   - Trim 5' of R1: Recommended if adapters are still on the input sequences. Trims the first J bases from the 5' end of each forward read.
   - Trim 5' of R2: Recommended if adapters are still on the input sequences. Trims the first K bases from the 5' end of each reverse read.
   - Truncate 3' of R1: Recommended if quality drops off along the length of the read. Clips the forward read starting M bases from the 5' end (before trimming).
   - Truncate 3' of R2: Recommended if quality drops off along the length of the read. Clips the reverse read starting N bases from the 5' end (before trimming).
   - Threads: Number of threads to use for steps that support multithreading.

  ### __Data Analysis Steps__
  1. Generate FASTX quality statistics for visualization of unmapped, raw FASTQ reads.
  2. Import the data, make a qiime artifact (demux.qza), and summary visualization process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.
  4. Generate a phylogenetic tree for diversity analyses and rarefaction processing and plotting.
  5. Taxonomy classification of amplicons. Performed using a Naive Bayes classifier trained on the Greengenes2 database "gg_2022_10_backbone_full_length.nb.qza".

  ### __References__
  1. Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9
    
