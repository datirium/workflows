cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement


'sd:upstream':
  genome_indices:
    - "bismark-index.cwl"


'sd:metadata':
  - "../metadata/chipseq-header.cwl"


inputs:

  fastq_file_r1:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "FASTQ file(s)"
    doc: "Uncompressed or gzipped FASTQ file(s), single-end"

  indices_folder:
    type: Directory
    label: "Bismark indices folder"
    doc: "Path to Bismark generated indices folder"
    'sd:upstreamSource': "genome_indices/indices_folder"

  chrom_length:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Chromosomes length file"
    doc: "Chromosomes length file"
    'sd:upstreamSource': "genome_indices/chrom_length"

  threads:
    type: int?
    default: 2
    label: "Number of Bowtie2/Trimmomatic/Samtools threads to use"
    doc: "Set the number of threads for Bowtie2, Trimmomatic, Samtools"
    'sd:layout':
      advanced: true


outputs:

  alignment_file:
    type: File
    label: "BAM alignment file"
    doc: "Bismark generated not sorted BAM alignment file"
    format: "http://edamontology.org/format_2572"
    outputSource: bismark_align/bam_file

  bambai_pair:
    type: File
    label: "BAM alignment and BAI index files"
    doc: "Bismark generated coordinate sorted BAM alignment and BAI index files"
    format: "http://edamontology.org/format_2572"
    outputSource: samtools_sort_index/bam_bai_pair
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        # optional: true
        type: 'alignment'
        format: 'bam'
        name: "BAM Track"
        displayMode: "SQUISHED"

  bismark_alignment_report:
    type: File
    label: "Bismark alignment and methylation report"
    doc: "Bismark generated alignment and methylation summary report"
    format: "http://edamontology.org/format_2330"
    outputSource: bismark_align/alignment_report

  chg_context_file:
    type: File
    label: "CHG methylation call"
    doc: "CHG methylation call"
    format: "http://edamontology.org/format_3475"
    outputSource: bismark_extract_methylation/chg_context_file

  chh_context_file:
    type: File
    label: "CHH methylation call"
    doc: "CHH methylation call"
    format: "http://edamontology.org/format_3475"
    outputSource: bismark_extract_methylation/chh_context_file

  cpg_context_file:
    type: File
    label: "CpG methylation call"
    doc: "CpG methylation call"
    format: "http://edamontology.org/format_3475"
    outputSource: bismark_extract_methylation/cpg_context_file

  mbias_plot:
    type: File
    label: "Methylation bias plot"
    doc: "QC data showing methylation bias across read lengths"
    format: "http://edamontology.org/format_3475"
    outputSource: bismark_extract_methylation/mbias_plot

  mbias_plot_png:
    type: File[]
    label: "Methylation bias plot (PNG)"
    doc: "QC data showing methylation bias across read lengths"
    format: "http://edamontology.org/format_3603"
    outputSource: bismark_extract_methylation/mbias_plot_png
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Methylation bias plot'

  bedgraph_coverage_file:
    type: File
    label: "Methylation statuses bedGraph coverage file"
    doc: "Coverage text file summarising cytosine methylation values in bedGraph format (tab-delimited; 0-based start coords, 1-based end coords)"
    format: "http://edamontology.org/format_3583"
    outputSource: bismark_extract_methylation/bedgraph_coverage_file
    # 'sd:visualPlugins':
    # - igvbrowser:
    #     tab: 'IGV Genome Browser'
    #     id: 'igvbrowser'
    #     type: 'annotation'
    #     name: "Methylation statuses"
    #     displayMode: "COLLAPSE"
    #     height: 40

  bigwig_coverage_file:
    type: File
    label: "Methylation statuses bigWig coverage file"
    doc: "Coverage text file summarising cytosine methylation values in bigWig format (tab-delimited; 0-based start coords, 1-based end coords)"
    format: "http://edamontology.org/format_3006"
    outputSource: sorted_bedgraph_to_bigwig/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Methylation statuses"
        height: 120

  bismark_coverage_file:
    type: File
    label: "Methylation statuses Bismark coverage file"
    doc: "Coverage text file summarising cytosine methylation values in Bismark format (tab-delimited, 1-based genomic coords)"
    format: "http://edamontology.org/format_3583"
    outputSource: bismark_extract_methylation/bismark_coverage_file

  genome_wide_methylation_report:
    type: File
    label: "Genome-wide cytosine methylation report"
    doc: "Genome-wide methylation report for all cytosines in the genome"
    format: "http://edamontology.org/format_3475"
    outputSource: bismark_extract_methylation/genome_wide_methylation_report

  splitting_report:
    type: File
    label: "Methylation extraction log"
    doc: "Log file giving summary statistics about methylation extraction"
    format: "http://edamontology.org/format_2330"
    outputSource: bismark_extract_methylation/splitting_report

  collected_report:
    type: File
    label: "Bismark Processing Report"
    doc: "Bismark generated graphical HTML report page"
    format: "http://edamontology.org/format_2331"
    outputSource: bismark_report/collected_report
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  collected_report_formatted:
    type: File
    label: "Combined Bismark alignment and splitting reports"
    doc: "Bismark generated alignment and splitting reports. Combined"
    format: "http://edamontology.org/format_3475"
    outputSource: format_bismark_report/collected_report_formatted
    'sd:visualPlugins':
    - tableView:
        vertical: true
        tab: 'Overview'
    'sd:preview':
      'sd:visualPlugins':
      - pie:
          colors: ['#b3de69', '#99c0db', '#fb8072', '#fdc381']
          data: [$2, $3, $4, $5]

  trim_adapters_report:
    type: File
    label: "TrimGalore report"
    doc: "TrimGalore generated report"
    format: "http://edamontology.org/format_2330"
    outputSource: trim_adapters/report_file


steps:

  extract_fastq_r1:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_r1
      output_prefix:
        default: "read_1"
    out:
    - fastq_file

  trim_adapters:
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq_r1/fastq_file
      dont_gzip:
        default: true
      length:
        default: 30
    out:
      - trimmed_file
      - report_file

  bismark_align:
    run: ../tools/bismark-align.cwl
    in:
      fastq_file_r1: trim_adapters/trimmed_file
      indices_folder: indices_folder
      processes:
        default: 1
      threads: threads
    out: [bam_file, alignment_report]

  bismark_extract_methylation:
    run: ../tools/bismark-extract-methylation.cwl
    in:
      genome_folder: indices_folder
      bam_file: bismark_align/bam_file
      processes:
        default: 1
    out:
      - chg_context_file
      - chh_context_file
      - cpg_context_file
      - mbias_plot
      - mbias_plot_png
      - bedgraph_coverage_file
      - bismark_coverage_file
      - genome_wide_methylation_report
      - splitting_report

  bismark_report:
    run: ../tools/bismark-report.cwl
    in:
      alignment_report: bismark_align/alignment_report
      splitting_report: bismark_extract_methylation/splitting_report
      mbias_report: bismark_extract_methylation/mbias_plot
    out: [collected_report]

  format_bismark_report:
    run: ../tools/python-get-stat-bismark.cwl
    in:
      alignment_report: bismark_align/alignment_report
      splitting_report: bismark_extract_methylation/splitting_report
    out: [collected_report_formatted]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bismark_align/bam_file
      sort_output_filename:
        default: "sorted_alignments.bam"
      threads: threads
    out: [bam_bai_pair]

  remove_metadata:
    run: ../tools/custom-bash.cwl
    in:
      input_file: bismark_extract_methylation/bedgraph_coverage_file
      script:
        default: |
          zcat "$0" | grep -v "track" > methylation_statuses.bedGraph
    out: [output_file]

  sort_bedgraph:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: remove_metadata/output_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  sorted_bedgraph_to_bigwig:
    run: ../tools/ucsc-bedgraphtobigwig.cwl
    in:
      bedgraph_file: sort_bedgraph/sorted_file
      chrom_length_file: chrom_length
      output_filename:
        source: sort_bedgraph/sorted_file
        valueFrom: $(self.basename.split('.').slice(0,-1).join('.') + ".bigWig")
    out: [bigwig_file]



$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Bismark Methylation SE"
label: "Bismark Methylation SE"
s:alternateName: "Bisulfite-Sequencing analysis pipeline to map bisulfite converted sequence reads and determine cytosine methylation states"

s:downloadUrl: https://github.com/Barski-lab/workflows/blob/master/workflows/bismark-methylation-se.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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
#   $include: ../descriptions/bismark-methylation-se.md


doc: |
  Sequence reads are first cleaned from adapters and transformed into fully bisulfite-converted forward (C->T) and reverse
  read (G->A conversion of the forward strand) versions, before they are aligned to similarly converted versions of the
  genome (also C->T and G->A converted). Sequence reads that produce a unique best alignment from the four alignment processes
  against the bisulfite genomes (which are running in parallel) are then compared to the normal genomic sequence and the
  methylation state of all cytosine positions in the read is inferred. A read is considered to align uniquely if an alignment
  has a unique best alignment score (as reported by the AS:i field). If a read produces several alignments with the same number
  of mismatches or with the same alignment score (AS:i field), a read (or a read-pair) is discarded altogether.

  On the next step we extract the methylation call for every single C analysed. The position of every single C will be written
  out to a new output file, depending on its context (CpG, CHG or CHH), whereby methylated Cs will be labelled as forward
  reads (+), non-methylated Cs as reverse reads (-). The output of the methylation extractor is then transformed into a bedGraph
  and coverage file. The bedGraph counts output is then used to generate a genome-wide cytosine report which reports the number
  on every single CpG (optionally every single cytosine) in the genome, irrespective of whether it was covered by any reads or not.
  As this type of report is informative for cytosines on both strands the output may be fairly large (~46mn CpG positions or >1.2bn
  total cytosine positions in the human genome).