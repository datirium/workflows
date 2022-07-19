class: Workflow
cwlVersion: v1.1


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  species:
    type:
    - type: enum
      symbols:
      - "hsa"
      - "mmu"
      - "rat"
    label: "Species (organism), as specified in library file or taxon id"
    doc: |
      Species (organism), as specified in library file or taxon id

  starting_material:
    type:
    - type: enum
      symbols:
      - "rna"
      - "dna"
    label: "Type of starting material"
    doc: |
      Type of starting material. This affects which part of reference V
      gene sequences will be used for alignment, with or without intron

  receptor_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "xcr"
      - "tcr"
      - "bcr"
      - "tra"
      - "trb"
      - "trg"
      - "trd"
      - "igh"
      - "igk"
      - "igl"
    default: "xcr"
    label: "Dedicated receptor type for analysis. Use all chains by default"
    doc: |
      Dedicated receptor type for analysis. By default, all T- and B-cell
      receptor chains are analyzed. MiXCR has special aligner kAligner2,
      which is used when B-cell receptor type is selected. Possible values
      for --receptor-type are: xcr (all chains), tcr, bcr, tra, trb, trg,
      trd, igh, igk, igl.

  fastq_file:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ file(s) (optionally compressed)"
    doc: "FASTQ file(s) (optionally compressed)"

  min_trim_length:
    type: int?
    default: 30
    label: "Discard reads that became shorter than this value after adapter trimming"
    doc: |
      Discard reads that became shorter than length INT because of either
      quality or adapter trimming. A value of '0' effectively disables
      this behaviour.
    'sd:layout':
      advanced: true

  contig_assembly:
    type: boolean?
    default: false
    label: "Assemble full receptor sequences"
    doc: |
      Whether to assemble full receptor sequences (assembleContigs).
      This option may slow down the computation.
    'sd:layout':
      advanced: true

  impute_germline_on_export:
    type: boolean?
    default: false
    label: "Use germline segments for uncovered gene features"
    doc: |
      Use germline segments (printed with lowercase letters) for
      uncovered gene features.
    'sd:layout':
      advanced: true

  only_productive:
    type: boolean?
    default: false
    label: "Filter out-of-frame sequences and clonotypes with stop-codons"
    doc: |
      Filter out-of-frame sequences and clonotypes with stop-codons
      in clonal sequence export
    'sd:layout':
      advanced: true

  assemble_partial_rounds:
    type: int?
    default: 2
    label: "Number of times to assemble overlapping fragmented sequencing reads into long-enough CDR3 containing contigs"
    doc: |
      Number of consequent assemblePartial executions (assembles overlapping
      fragmented sequencing reads into long-enough CDR3 containing contigs)
    'sd:layout':
      advanced: true

  do_not_extend_alignments:
    type: boolean?
    default: false
    label: "Do not perform extension of good TCR alignments"
    doc: |
      Do not perform extension of good TCR alignments
    'sd:layout':
      advanced: true

  unweighted:
    type: boolean?
    default: false
    label: "Do not normalize all statistics by clonotype frequency"
    doc: |
      If not set, all statistics will be weighted by clonotype frequency
    'sd:layout':
      advanced: true

  amino_acid:
    type: boolean?
    default: false
    label: "Use CDR3 amino acid sequences for spectratype calculation instead of nucleotide ones"
    doc: |
      Use CDR3 amino acid sequences for spectratype calculation instead of nucleotide ones
    'sd:layout':
      advanced: true

  top:
    type: int?
    default: 5
    label: "Number of top clonotypes to visualize. Should not exceed 10"
    doc: "Number of top clonotypes to visualize. Should not exceed 10"
    'sd:layout':
      advanced: true

  intersect_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "strict"
      - "nt"
      - "ntV"
      - "ntVJ"
      - "aa"
      - "aaV"
      - "aaVJ"
      - "aa!nt"
    default: "strict"
    label: "Clonotype features to compare when checking if two clonotypes match"
    doc: |
      Specifies which clonotype features (CDR3 sequence, V/J segments, hypermutations)
      will be compared when checking if two clonotypes match.
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  fastqc_report:
    type: File
    outputSource: fastqc_fastq_file/html_file
    label: "FastqQC report for FASTQ file"
    doc: |
      FastqQC report for FASTQ file
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  trim_adapters_report:
    type: File
    label: "Adapter trimming report for FASTQ file"
    doc: "Adapter trimming report for FASTQ file"
    outputSource: trim_adapters/report_file

  mixcr_alignment_file:
    type: File
    label: "MiXCR alignment file"
    doc: "MiXCR alignment file"
    outputSource: mixcr_shotgun/alignment_file

  mixcr_assembled_clonotypes_file:
    type: File
    label: "MiXCR assembled clonotypes file"
    doc: "MiXCR assembled clonotypes file"
    outputSource: mixcr_shotgun/assembled_clonotypes_file

  mixcr_all_clonotypes_file:
    type: File
    label: "MiXCR all clonotypes file"
    doc: "MiXCR all clonotypes file"
    outputSource: mixcr_shotgun/all_clonotypes_file

  mixcr_report_file:
    type: File
    label: "MiXCR report file"
    doc: "MiXCR report file"
    outputSource: mixcr_shotgun/report_file

  vdj_file:
    type: File
    label: "All clonotypes file from MiXCR converted to VDJTools format"
    doc: "All clonotypes file from MiXCR converted to VDJTools format"
    outputSource: vdjtools_convert_mixcr/vdj_file

  vdj_basic_stats_file:
    type: File
    label: "A set of basic sample statistics, such as read counts, number of clonotypes, etc."
    doc: "A set of basic sample statistics, such as read counts, number of clonotypes, etc."
    outputSource: vdjtools_calc_basic_stats/basic_stats_file

  vdj_spectratype_insert_wt_file:
    type: File
    label: "VDJTools spectratype"
    doc: "VDJTools spectratype"
    outputSource: vdjtools_calc_spectratype/spectratype_insert_wt_file

  vdj_spectratype_ndn_wt_file:
    type: File
    label: "VDJTools spectratype"
    doc: "VDJTools spectratype"
    outputSource: vdjtools_calc_spectratype/spectratype_ndn_wt_file

  vdj_spectratype_aa_wt_file:
    type: File?
    label: "VDJTools spectratype"
    doc: "VDJTools spectratype"
    outputSource: vdjtools_calc_spectratype/spectratype_aa_wt_file

  vdj_spectratype_nt_wt_file:
    type: File?
    label: "VDJTools spectratype"
    doc: "VDJTools spectratype"
    outputSource: vdjtools_calc_spectratype/spectratype_nt_wt_file
  
  vdj_fancy_spectratype_file:
    type: File
    label: "VDJTools fancy spectratype"
    doc: "VDJTools fancy spectratype"
    outputSource: vdjtools_plot_fancy_spectratype/fancy_spectratype_file

  vdj_fancy_spectratype_plot:
    type: File
    label: "VDJTools fancy spectratype"
    doc: "VDJTools fancy spectratype"
    outputSource: vdjtools_plot_fancy_spectratype/fancy_spectratype_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'VDJTools spectratype'

  vdj_diversity_stats_file:
    type: File
    label: "VDJTools diversity statistics"
    doc: "VDJTools diversity statistics"
    outputSource: vdjtools_calc_diversity_stats/diversity_stats_file

  vdj_quantile_stats_file:
    type: File
    label: "VDJTools quantile statistics"
    doc: "VDJTools quantile statistics"
    outputSource: vdjtools_plotquantilestats/quantile_stats_file

  vdj_quantile_stats_plot:
    type: File
    label: "VDJTools quantile statistics"
    doc: "VDJTools quantile statistics"
    outputSource: vdjtools_plotquantilestats/quantile_stats_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'VDJTools quantile statistics'


steps:

  extract_fastq_file:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
      output_prefix:
        default: "read_1"
    out:
    - fastq_file

  fastqc_fastq_file:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_file/fastq_file
      threads: threads
    out:
    - html_file

  trim_adapters:
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq_file/fastq_file
      dont_gzip:
        default: true
      length: min_trim_length
    out:
    - trimmed_file
    - report_file

  bypass_trim:
    run: ../tools/bypass-trimgalore-se.cwl
    in:
      original_fastq_file: extract_fastq_file/fastq_file
      trimmed_fastq_file: trim_adapters/trimmed_file
      trimming_report_file: trim_adapters/report_file
      min_reads_count:
        default: 100
    out:
    - selected_fastq_file
    - selected_report_file

  mixcr_shotgun:
    run: ../tools/mixcr-shotgun.cwl
    in:
      species: species
      starting_material: starting_material
      receptor_type: receptor_type
      contig_assembly: contig_assembly
      impute_germline_on_export: impute_germline_on_export
      only_productive: only_productive
      assemble_partial_rounds: assemble_partial_rounds
      do_not_extend_alignments: do_not_extend_alignments
      threads: threads
      fastq_file: bypass_trim/selected_fastq_file
      output_prefix:
        default: "mixcr"
    out:
    - alignment_file
    - assembled_clonotypes_file
    - all_clonotypes_file
    - report_file

  vdjtools_convert_mixcr:
    run: ../tools/vdjtools-convert-mixcr.cwl
    in:
      clonotypes_file: mixcr_shotgun/all_clonotypes_file
    out:
    - vdj_file

  vdjtools_calc_basic_stats:
    run: ../tools/vdjtools-calc-basic-stats.cwl
    in:
      vdj_file: vdjtools_convert_mixcr/vdj_file
      unweighted: unweighted
    out:
    - basic_stats_file

  vdjtools_calc_spectratype:
    run: ../tools/vdjtools-calc-spectratype.cwl
    in:
      vdj_file: vdjtools_convert_mixcr/vdj_file
      unweighted: unweighted
      amino_acid: amino_acid
    out:
    - spectratype_insert_wt_file
    - spectratype_ndn_wt_file
    - spectratype_aa_wt_file
    - spectratype_nt_wt_file

  vdjtools_plot_fancy_spectratype:
    run: ../tools/vdjtools-plot-fancy-spectratype.cwl
    in:
      vdj_file: vdjtools_convert_mixcr/vdj_file
      top: top
    out:
    - fancy_spectratype_file
    - fancy_spectratype_plot

  vdjtools_calc_diversity_stats:
    run: ../tools/vdjtools-calc-diversity-stats.cwl
    in:
      vdj_file: vdjtools_convert_mixcr/vdj_file
      intersect_type: intersect_type
    out:
    - diversity_stats_file

  vdjtools_plotquantilestats:
    run: ../tools/vdjtools-plot-quantile-stats.cwl
    in:
      vdj_file: vdjtools_convert_mixcr/vdj_file
      top: top
    out:
    - quantile_stats_file
    - quantile_stats_plot


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Repertoire Sequencing Analysis (single-read)"
label: "Repertoire Sequencing Analysis (single-read)"
s:alternateName: "Runs MiXCR and VDJTools for analysis of T- and B- cell receptor repertoire sequencing data (single-read)"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/mixcr-se.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Repertoire Sequencing Analysis (single-read)
  ===========================================
  
  Runs MiXCR and VDJTools for analysis of T- and B- cell receptor
  repertoire sequencing data (single-read)