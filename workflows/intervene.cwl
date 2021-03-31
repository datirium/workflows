cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement


'sd:upstream':
  filtered_experiment:
  - "filter-peaks-for-heatmap.cwl"
  genome_indices:
  - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  intervals_files:
    type: File[]
    format: "http://edamontology.org/format_3003"
    label: "Filtered ChIP/ATAC experiment"
    doc: |
      Filtered peaks file from ChIP/ATAC experiment formatted
      as headerless BED [chrom start end name] file
    'sd:upstreamSource': "filtered_experiment/filtered_file"
    'sd:localLabel': true

  intervals_aliases:
    type: string[]
    label: "Filtered ChIP/ATAC experiment"
    doc: "Filtered ChIP/ATAC experiment alias"
    'sd:upstreamSource': "filtered_experiment/filtered_file_alias"

  annotation_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Genome to use for the nearest gene assignment"
    doc: "Tab-separated annotation file"
    'sd:upstreamSource': "genome_indices/annotation"
    'sd:localLabel': true

  promoter_dist:
    type: int?
    default: 1000
    label: "Max distance from gene TSS (in both direction) to assign interval to promoter"
    doc: |
      Max distance from gene TSS (in both direction) overlapping
      which the interval will be assigned to the promoter region
    'sd:layout':
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Max distance from promoter (only in upstream direction) to assign interval to upstream"
    doc: |
      Max distance from the promoter (only in upstream direction)
      overlapping which the interval will be assigned to the
      upstream region
    'sd:layout':
      advanced: true


outputs:

  overlapped_intervals_files:
    type: File[]
    format: "http://edamontology.org/format_3475"
    outputSource: annotate_overlapped_intervals/result_file
    label: "Overlapped intervals files with the nearest genes assigned"
    doc: "Overlapped intervals files with the nearest genes assigned"

  overlapped_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    outputSource: overlap_intervals/overlapped_plot
    label: "Intervals overlap diagram"
    doc: "Intervals overlap diagram"
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Intervals overlap diagram'

  overlapped_combinations:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: overlap_intervals/overlapped_combinations
    label: "Overlapped combinations file"
    doc: "Overlapped combinations file"

  intervene_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "intervene stdout log"
    doc: "intervene stdout log"
    outputSource: overlap_intervals/stdout_log

  intervene_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "intervene stderr log"
    doc: "intervene stderr log"
    outputSource: overlap_intervals/stderr_log


steps:

  sort_intervals:
    run:
      cwlVersion: v1.0
      class: Workflow
      requirements:
      - class: ScatterFeatureRequirement
      inputs:
        not_sorted_intervals_files:
          type: File[]
      outputs:
        sorted_intervals_files:
          type: File[]
          outputSource: batch/output_file
      steps:
        batch:
          run: ../tools/custom-bash.cwl
          in:
            input_file: not_sorted_intervals_files
            script:
              default: |
                cat "$0" | tr -d '\r' | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n > `basename $0`
          scatter: input_file                
          out:
          - output_file
    in:
      not_sorted_intervals_files: intervals_files
    out:
    - sorted_intervals_files

  merge_intervals:
    run:
      cwlVersion: v1.0
      class: Workflow
      requirements:
      - class: ScatterFeatureRequirement
      inputs:
        bed_file:
          type: File[]
      outputs:
        merged_bed_file:
          type: File[]
          outputSource: batch/merged_bed_file
      steps:
        batch:
          run: ../tools/bedtools-merge.cwl
          in:
            bed_file: bed_file
          scatter: bed_file
          out:
          - merged_bed_file    
    in:
      bed_file: sort_intervals/sorted_intervals_files
    out:
    - merged_bed_file

  overlap_intervals:
    run: ../tools/intervene.cwl
    in:
      intervals_files: merge_intervals/merged_bed_file
      intervals_aliases: intervals_aliases
      figure_format:
        default: "png"
      diagram_type:
        source: merge_intervals/merged_bed_file
        valueFrom: |
          ${
            if (self.length > 6){
              return "upset";
            }
            return "venn";
          }
    out:
    - overlapped_intervals_files
    - overlapped_plot
    - overlapped_combinations
    - stdout_log
    - stderr_log

  refactore_overlapped_intervals:
    run:
      cwlVersion: v1.0
      class: Workflow
      requirements:
      - class: ScatterFeatureRequirement
      inputs:
        input_file:
          type: File[]
      outputs:
        output_file:
          type: File[]
          outputSource: batch/output_file
      steps:
        batch:
          run: ../tools/custom-bash.cwl
          in:
            input_file: input_file
            script:
              default: >
                cat $0 | awk
                'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
                {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t0\t0\t0\t0\t0"}' > `basename $0`
          scatter: input_file
          out:
          - output_file
    in:
      input_file: overlap_intervals/overlapped_intervals_files
    out:
    - output_file

  annotate_overlapped_intervals:
    run:
      cwlVersion: v1.0
      class: Workflow
      requirements:
      - class: ScatterFeatureRequirement
      inputs:
        input_filename:
          type: File[]
        annotation_filename:
          type: File
        promoter_bp:
          type: int?
        upstream_bp:
          type: int?
      outputs:
        result_file:
          type: File[]
          outputSource: batch/result_file
      steps:
        batch:
          run: ../tools/iaintersect.cwl
          in:
            input_filename: input_filename
            annotation_filename: annotation_filename
            promoter_bp: promoter_bp
            upstream_bp: upstream_bp
          scatter: input_filename
          out:
          - result_file
    in:
      input_filename: refactore_overlapped_intervals/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out:
    - result_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Genomic regions intersection and visualization"
label: "Genomic regions intersection and visualization"
s:alternateName: "Genomic regions intersection and visualization"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/intervene.cwl
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
  
  Genomic regions intersection and visualization
  ==============================================
  
  1. Merges intervals within each of the filtered peaks files from ChIP/ATAC experiments
  2. Overlaps merged intervals and assigns the nearest genes to them
