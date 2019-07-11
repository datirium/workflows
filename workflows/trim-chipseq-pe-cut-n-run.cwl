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
    - "chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-pe.cwl"
  control_experiment:
    - "chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-pe.cwl"
  genome_indices:
    - "genome-indices.cwl"


inputs:

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

  bambai_pair:
    type: File
    secondaryFiles:
      - .bai
    'sd:upstreamSource': "base_experiment/bambai_pair"
    format: "http://edamontology.org/format_2572"
    label: "ChIP-Seq paired-end experiment"
    doc: "Coordinate sorted BAM alignment and index BAI files"

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
    type: boolean
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

  exp_fragment_size:
    type: int?
    default: 150
    label: "Expected fragment size"
    doc: "Expected fragment size for read extenstion towards 3' end if force_fragment_size was set to True or if calculated by MACS2 fragment size was less that 80 bp"

  force_fragment_size:
    type: boolean?
    default: false
    label: "Force peak calling with expected fragment size"
    doc: "Make MACS2 don't build the shifting model and use expected fragment size for read extenstion towards 3' end"

  threads:
    type: int?
    default: 2
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

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
        type: 'alignment'
        format: 'bam'
        name: "Filtered Nucleotide Sequence Alignments"
        displayMode: "SQUISHED"

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
        type: 'bed'
        name: "Narrow peaks"
        height: 120

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
        type: 'bed'
        name: "Broad peaks"
        height: 120

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
      output_filename:
        source: bambai_pair
        valueFrom: $(self.basename.split('.').slice(0,-1).join('.')+'_filtered.bam')
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
    run: ../tools/macs2-callpeak-biowardrobe-only.cwl
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
      nomodel: force_fragment_size
      extsize: exp_fragment_size
      bw: exp_fragment_size
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
      - macs2_fragments_calculated

  bam_to_bigwig:
    run: ../subworkflows/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      chrom_length_file: chrom_length
      mapped_reads_number: get_statistics/mapped_reads
      fragment_size: macs2_callpeak/macs2_fragments_calculated
      pairchip:
        default: true
    out: [bigwig_file]

  island_intersect:
      run: ../tools/iaintersect.cwl
      in:
        input_filename: macs2_callpeak/peak_xls_file
        annotation_filename: annotation_file
        promoter_bp:
          default: 1000
      out: [result_file, log_file]

  average_tag_density:
      run: ../tools/atdp.cwl
      in:
        input_file: samtools_sort_index/bam_bai_pair
        annotation_filename: annotation_file
        fragmentsize_bp: macs2_callpeak/macs2_fragments_calculated
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
- http://schema.org/docs/schema_org_rdfa.html

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

doc: |
  Experimental pipeline for Cut-n-Run analysis. Uses mapping results from the following experiment types:
  
  - `chipseq-pe.cwl`
  - `trim-chipseq-pe.cwl`
  - `trim-atacseq-pe.cwl`
  
  Note, the upstream analyses should not have duplicates removed
