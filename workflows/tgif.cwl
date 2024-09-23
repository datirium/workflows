cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

'sd:upstream':
  genome_indices: "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  fastq_file:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ input file(s)"
    format: "http://edamontology.org/format_1930"
    doc: "WGS or NCATS enriched sequencing library reads data in a FASTQ format (ONT expected but PacBio/Illumina also accepted)"
    'sd:localLabel': true

  read_type:
    type:
    - "null"
    - type: enum
      name: "read type"
      symbols:
      - map-ont
      - map-pb
      - sr
    default: "map-ont"
    label: "Sequencing read type"
    doc: "Sequencing read type for minimap2: map-ont (Oxford Nanopore), map-pb (PacBio), sr (short read/Illumina)"

  reference_fasta:
    type: File?
    format: "http://edamontology.org/format_1929"
    'sd:upstreamSource': "genome_indices/fasta_output"
    label: "Reference Genome FASTA"
    doc: "Reads will be aligned against this reference. Select from the dropdown if not uploading your own custom reference genome FASTA."

  plasmid_fasta:
    type: File?
    format: "http://edamontology.org/format_1929"
    default: null
    label: "Vector/Plasmid FASTA"
    doc: "Reads will be aligned against this sequence along with the reference genome. File can only contain 2 lines maximum, a standard FASTA header for line 1, and a single unbroken DNA sequence string for line 2."

  yn_igv_output:
    type:
    - "null"
    - type: enum
      name: "igvoutput"
      symbols:
      - y
      - n
    default: "y"
    label: "Include alignment output for IGV?"
    doc: "Output sorted bam of downselected reads for IGV? Default: y"
    'sd:layout':
      advanced: true

  yn_plot_output:
    type:
    - "null"
    - type: enum
      name: "plot output"
      symbols:
      - y
      - n
    default: "y"
    label: "Include pileup plots output?"
    doc: "Output read pileup plot per probable insertion site? Default: y"
    'sd:layout':
      advanced: true

  threads_count:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  tgif_insertions_all:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Insertion Sites (all)"
    doc: "All found insertion sites converted to tsv format for scidap table"
    outputSource: run_tgif_ncats/tgif_insertions_all

  tgif_insertions_filtered:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Insertion Sites (all)"
    doc: "All found insertion sites converted to tsv format for scidap table"
    outputSource: run_tgif_ncats/tgif_insertions_filtered
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Insertion Sites (filtered)"
        Title: "Insertion Sites (filtered)"

  insertion_site_plots:
    type: File
    label: "TAR with plots of probable insertion sites"
    doc: "TAR with plots of probable insertion sites"
    outputSource: run_tgif_ncats/insertion_site_plots

  alignment_files:
    type: File
    label: "Compressed TAR with sorted bam and bai files against reference and vector sequences"
    doc: "Compressed TAR with sorted bam and bai files against reference and vector sequences"
    outputSource: run_tgif_ncats/insertion_site_plots

  tgif_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout log for tgif_ncats tool"
    doc: "stdout log for tgif_ncats tool"
    outputSource: run_tgif_ncats/stdout_log

  tgif_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr log for tgif_ncats tool"
    doc: "stderr log for tgif_ncats tool"
    outputSource: run_tgif_ncats/stderr_log

  primer3_output:
    type: File
    label: "TAR containing primer3 output"
    doc: "TAR containing primer3 output"
    outputSource: run_tgif_primer3/primer3_output

  primer3_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout log for primer3 tool"
    doc: "stdout log for primer3 tool"
    outputSource: run_tgif_primer3/stdout_log

  primer3_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr log for primer3 tool"
    doc: "stderr log for primer3 tool"
    outputSource: run_tgif_primer3/stderr_log

  tgif_summary:
    type: File?
    label: "Markdown formatted table with summary stats"
    format: "http://edamontology.org/format_3835"
    doc: "Markdown formatted table with summary stats"
    outputSource: run_tgif_summary/summary_file
    'sd:visualPlugins':
    - markdownView:
        tab: "Overview"

  summary_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout log for summary tool"
    doc: "stdout log for summary tool"
    outputSource: run_tgif_summary/stdout_log

  summary_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr log for summary tool"
    doc: "stderr log for summary tool"
    outputSource: run_tgif_summary/stderr_log


steps:

  run_tgif_ncats:
    run: ../tools/tgif-ncats.cwl
    in:
      input_fastq: fastq_file
      read_type: read_type
      reference_fasta: reference_fasta
      plasmid_fasta: plasmid_fasta
      yn_igv_output: yn_igv_output
      yn_plot_output: yn_plot_output
      threads_count: threads_count
    out:
      - tgif_insertions_all
      - tgif_insertions_filtered
      - insertion_site_plots
      - alignment_files
      - script_log_file
      - stdout_log
      - stderr_log

  run_tgif_primer3:
    run: ../tools/tgif-primer3.cwl
    in:
      insertions_filtered: run_tgif_ncats/tgif_insertions_filtered
      reference_fasta: reference_fasta
      threads_count: threads_count
    out:
      - primer3_output
      - stdout_log
      - stderr_log

  run_tgif_summary:
    run: ../tools/tgif-summary.cwl
    in:
      tgif_ncats_log_file: run_tgif_ncats/script_log_file
      insertions_all: run_tgif_ncats/tgif_insertions_all
      insertions_filtered: run_tgif_ncats/tgif_insertions_filtered
    out:
      - summary_file
      - stdout_log
      - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "TgIF - Transgene Insertion Finder"
label: "TgIF - Transgene Insertion Finder"
s:alternateName: "TgIF - Transgene Insertion Finder"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/tgif.cwl
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
  TgIF (trans-gene insertion finder)
  ==============================================
  The TgIF algorithm returns a list of probable insertion sites in a target organism. It requires the user to provided a FASTQ file of ONT (Oxford Nanopore Technologies) reads (-f), the reference FASTA of the trans-gene (Tg) vector containing the insertion sequence (-i), and the reference FASTA of the target organism (-r). The algorithm is tailored for ONT reads from a modified nCATS[1] (nanopore Cas9-targeted sequencing) enriched library, however the algorithm will also produce informative results from a FASTQ derived from WGS (shotgun) sequencing libraries. The modified nCATS method is described here, and a brief overview can be found below.

  The basic workflow of TgIF is alignment (using minimap2[2]) of reads (-f) to a combined reference of the Tg vector (containing the desired insertion sequence) and target organism (ie. -i and -r are concatenated), and then searching for valleys (or gaps) in the resulting pileup of reads that map to both references at MAPQ>=30. A starting position (ps) of a valley is where the depth (d) at dp=0 and dp-1>0, an ending position (pe) of a valley is where the depth at dp=0 and dp+1>0, and a potential insertion scar is the gap between and including ps and pe.

  Primary Output files:
  - insertions_all.tsv, all probable insertion sites identified from the input fastq data
  - insertions_filtered.tgif, filtered sites that are most probable based on logic (4) above
  - reportsummary.md, summary of alignment metrics and insertion sites found

  Secondary Output files:
  - insertion_site_plots.tar, package of probable insertion site pileup plots
  - alignment_files.tar.gz, contains bam/bai for visualizing aligned reads to reference genome and vector sequence
  - primer3.tar, contains F/R primers for each filtered insertion site designed by primer3

  Documents
  ==============================================
  - github Page: https://github.com/jhuapl-bio/TgIF/tree/main

  References
  ==============================================
  - Gilpatrick, T. et al. Targeted nanopore sequencing with Cas9-guided adapter ligation. Nature Biotechnology 38, 433–438 (2020).
  - Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191
  - O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.
  - Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9 [19505943]
  - R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
  - H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
  - Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M and Rozen SG. Primer3--new capabilities and interfaces. Nucleic Acids Res. 2012 Aug 1;40(15):e115.