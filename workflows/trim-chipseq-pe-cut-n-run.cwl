cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

'sd:metadata':
  - "../metadata/chipseq-header.cwl"

'sd:upstream':
  base_experiment:
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-pe.cwl"
  control_experiment:
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-pe.cwl"
  genome_indices:
    - "genome-indices.cwl"


inputs:

  bambai_pair:
    type: File
    secondaryFiles:
      - .bai
    'sd:upstreamSource': "base_experiment/bambai_pair"
    'sd:localLabel': true
    format: "http://edamontology.org/format_2572"
    label: "ChIP-Seq paired-end experiment"
    doc: "Coordinate sorted BAM alignment and index BAI files"

  alignment_log:
    type: File
    'sd:upstreamSource': "base_experiment/bowtie_log"
    format: "http://edamontology.org/format_2330"
    label: "Read alignment log"
    doc: "Read alignment log file from Bowtie"

  rmdup_log:
    type: File
    'sd:upstreamSource': "base_experiment/samtools_rmdup_log"
    format: "http://edamontology.org/format_2330"
    label: "Remove duplicates log"
    doc: "Remove duplicates log file from Samtools"

  annotation_file:
    type: File
    'sd:upstreamSource': "genome_indices/annotation"
    label: "Genome annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Genome annotation file in TSV format"

  genome_size:
    type: string
    'sd:upstreamSource': "genome_indices/genome_size"
    label: "Effective genome size"
    doc: "The length of the mappable genome (hs, mm, ce, dm or number, for example 2.7e9)"

  chrom_length:
    type: File
    'sd:upstreamSource': "genome_indices/chrom_length"
    label: "Chromosome lengths file"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosome lengths file in TSV format"

  control_file:
    type: File?
    default: null
    'sd:upstreamSource': "control_experiment/bambai_pair"
    'sd:localLabel': true
    label: "Control ChIP-Seq paired-end experiment"
    format: "http://edamontology.org/format_2572"
    doc: "Indexed BAM file from the ChIP-Seq paired-end experiment to be used as a control for MACS2 peak calling"

  broad_peak:
    type: boolean?
    default: False
    label: "Call broad peaks"
    doc: "Make MACS2 call broad peaks by linking nearby highly enriched regions"

  min_fragment_size:
    type: int
    label: "Minimum fragment size"
    doc: "The minimum fragment size needed for read/pair inclusion"

  max_fragment_size:
    type: int
    label: "Maximum fragment size"
    doc: "The maximum fragment size needed for read/pair inclusion"

  promoter_dist:
    type: int?
    default: 1000
    'sd:layout':
      advanced: true
    label: "Max distance from gene TSS (in both direction) overlapping which the peak will be assigned to the promoter region"
    doc: "Max distance from gene TSS (in both direction) overlapping which the peak will be assigned to the promoter region"

  upstream_dist:
    type: int?
    default: 20000
    'sd:layout':
      advanced: true
    label: "Max distance from the promoter (only in upstream direction) overlapping which the peak will be assigned to the upstream region"
    doc: "Max distance from the promoter (only in upstream direction) overlapping which the peak will be assigned to the upstream region"

  threads:
    type: int?
    default: 2
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "Genome coverage"
    doc: "Genome coverage in bigWig format"
    outputSource: bam_to_bigwig/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Genome Coverage"
        height: 120

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Filtered aligned reads"
    doc: "Coordinate sorted filtered BAM alignment and index BAI files"
    outputSource: samtools_sort_index/bam_bai_pair
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        optional: true
        type: 'alignment'
        format: 'bam'
        name: "Filtered Nucleotide Sequence Alignments"
        displayMode: "SQUISHED"

  iaintersect_result:
    type: File
    label: "Gene annotated peaks"
    format: "http://edamontology.org/format_3475"
    doc: "MACS2 peak file annotated with nearby genes"
    outputSource: island_intersect/result_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Peak Calling'
        Title: 'Peak Coordinates'

  atdp_result:
    type: File
    label: "Average Tag Density Plot"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density Plot file in TSV format"
    outputSource: average_tag_density/result_file
    'sd:visualPlugins':
    - scatter:
        tab: 'QC Plots'
        Title: 'Average Tag Density Plot'
        xAxisTitle: 'Distance From TSS (bp)'
        yAxisTitle: 'Average Tag Density (per bp)'
        colors: ["#b3de69"]
        height: 500
        data: [$1, $2]
        comparable: "atdp"

  macs2_called_peaks:
    type: File
    label: "Called peaks"
    format: "http://edamontology.org/format_3468"
    doc: "Called peaks file with 1-based coordinates in XLS format"
    outputSource: macs2_callpeak/peak_xls_file

  macs2_narrow_peaks:
    type: File?
    label: "Narrow peaks"
    format: "http://edamontology.org/format_3613"
    doc: "Called peaks file in ENCODE narrow peak format"
    outputSource: macs2_callpeak/narrow_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Narrow peaks"
        displayMode: "COLLAPSE"
        height: 40

  macs2_broad_peaks:
    type: File?
    label: "Broad peaks"
    format: "http://edamontology.org/format_3614"
    doc: "Called peaks file in ENCODE broad peak format"
    outputSource: macs2_callpeak/broad_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Broad peaks"
        displayMode: "COLLAPSE"
        height: 40

  alignmentsieve_log:
    type: File
    label: "Alignment filtering log"
    format: "http://edamontology.org/format_2330"
    doc: "Alignment filtering log from deepTool's alignmentSieve"
    outputSource: filter_bam/alignmentsieve_log


steps:

  get_statistics:
    run: ../tools/python-get-stat-chipseq.cwl
    in:
      bowtie_log: alignment_log
      rmdup_log: rmdup_log
    out: [mapped_reads]

  filter_bam:
    run: ../tools/deeptools-alignmentsieve.cwl
    in:
      bambai_pair: bambai_pair
      min_fragment_length: min_fragment_size
      max_fragment_length: max_fragment_size
      threads: threads
    out:
      - filtered_bam_file
      - alignmentsieve_log

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: filter_bam/filtered_bam_file
      threads: threads
    out: [bam_bai_pair]

  macs2_callpeak:
    label: "Peak detection"
    doc: |
      Identifies enriched with aligned reads genome areas. Those areas correspond to the
      transcription factor binding sites.
    run: ../tools/macs2-callpeak.cwl
    in:
      treatment_file: samtools_sort_index/bam_bai_pair
      control_file: control_file
      nolambda:
        source: control_file
        valueFrom: $(!self)
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format_mode:
        default: BAMPE
      buffer_size:
        default: 10000
    out:
      - peak_xls_file
      - narrow_peak_file
      - broad_peak_file

  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      chrom_length_file: chrom_length
      mapped_reads_number: get_statistics/mapped_reads
      pairchip:
        default: true
    out: [bigwig_file]

  island_intersect:
    label: "Peak annotation"
    doc: |
      Assigns nearest genes to peaks to explore the biological implication of the open
      chromatin binding sites.
    run: ../tools/iaintersect.cwl
    in:
      input_filename: macs2_callpeak/peak_xls_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out: [result_file, log_file]

  average_tag_density:
    label: "Read enrichment around genes TSS"
    doc: |
      Generates average tag density plot around genes TSS as a lot of cis-regulatory
      elements are close to the TSS of their targets.
    run: ../tools/atdp.cwl
    in:
      input_file: samtools_sort_index/bam_bai_pair
      annotation_filename: annotation_file
      avd_window_bp:
        default: 5000
      avd_smooth_bp:
        default: 50
      ignore_chr:
        default: chrM
      double_chr:
        default: "chrX chrY"
      avd_heat_window_bp:
        default: 200
      mapped_reads: get_statistics/mapped_reads
    out: [result_file, log_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Cut-n-Run pipeline paired-end"
label: "Cut-n-Run pipeline paired-end"
s:alternateName: "Cut-n-Run analysis workflow for paired-end data"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/datirium/trim-chipseq-pe.cwl
s:codeRepository: https://github.com/datirium/workflows
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


# doc:
#   $include: ../descriptions/trim-chipseq-pe-cut-n-run.md


doc: |
  Experimental pipeline for Cut-n-Run analysis. Uses mapping results from the following experiment types:

  - `chipseq-pe.cwl`
  - `trim-chipseq-pe.cwl`
  - `trim-atacseq-pe.cwl`

  Note, the upstream analyses should not have duplicates removed