cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  filtered_experiment_a:
  - "filter-peaks-for-heatmap.cwl"
  - "super-enhancer.cwl"
  filtered_experiment_b:
  - "filter-peaks-for-heatmap.cwl"  
  - "super-enhancer.cwl"
  genome_indices:
  - "genome-indices.cwl"
  - "https://github.com/datirium/workflows/workflows/genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  intervals_file_a:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filtered ChIP/ATAC experiment A"
    doc: |
      Filtered peaks file from ChIP/ATAC experiment formatted
      as headerless BED [chrom start end name] file
    'sd:upstreamSource': "filtered_experiment_a/filtered_file"
    'sd:localLabel': true

  intervals_file_b:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filtered ChIP/ATAC experiment B"
    doc: |
      Filtered peaks file from ChIP/ATAC experiment formatted
      as headerless BED [chrom start end name] file
    'sd:upstreamSource': "filtered_experiment_b/filtered_file"
    'sd:localLabel': true

  annotation_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Genome type for intervals annotation"
    doc: "Tab-separated annotation file to assign the nearest genes"
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

  annotated_unique_from_a:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: groom_unique_from_a/output_file
    label: "Annotated intervals unique for experiment A"
    doc: "Annotated intervals unique for experiment A"
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Unique from A"
        Title: "Unique intervals from experiment A"

  unique_from_a_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    outputSource: sort_unique_from_a/sorted_file
    label: "BED file with unique for experiment A intervals"
    doc: "BED file with unique for experiment A intervals"
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Unique from A"
        displayMode: "COLLAPSE"
        height: 40

  annotated_unique_from_b:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: groom_unique_from_b/output_file
    label: "Annotated intervals unique for experiment B"
    doc: "Annotated intervals unique for experiment B"
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Unique from B"
        Title: "Unique intervals from experiment B"

  unique_from_b_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    outputSource: sort_unique_from_b/sorted_file
    label: "BED file with unique for experiment B intervals"
    doc: "BED file with unique for experiment B intervals"
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Unique from B"
        displayMode: "COLLAPSE"
        height: 40

  annotated_overlapped_from_a:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: groom_overlapped_from_a/output_file
    label: "Annotated intervals from experiment A overlapped with experiment B"
    doc: "Annotated intervals from experiment A overlapped with experiment B"
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Overlapped from A"
        Title: "Overlapped intervals from experiment A"

  overlapped_from_a_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    outputSource: sort_overlapped_from_a/sorted_file
    label: "BED file with intervals from experiment A overlapped with experiment B"
    doc: "BED file with intervals from experiment A overlapped with experiment B"
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Overlapped from A"
        displayMode: "COLLAPSE"
        height: 40

  annotated_overlapped_from_b:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: groom_overlapped_from_b/output_file
    label: "Annotated intervals from experiment B overlapped with experiment A"
    doc: "Annotated intervals from experiment B overlapped with experiment A"
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Overlapped from B"
        Title: "Overlapped intervals from experiment B"

  overlapped_from_b_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    outputSource: sort_overlapped_from_b/sorted_file
    label: "BED file with intervals from experiment B overlapped with experiment A"
    doc: "BED file with intervals from experiment B overlapped with experiment A"
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Overlapped from B"
        displayMode: "COLLAPSE"
        height: 40

  annotated_merged_overlapped_from_a_and_b:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: groom_merged_overlapped_from_a_and_b/output_file
    label: "Annotated merged overlapped intervals for experiments A and B"
    doc: "Annotated merged overlapped intervals for experiments A and B"
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Merged overlapped from A and B"
        Title: "Merged overlapped intervals from experiments A and B"

  merged_overlapped_from_a_and_b_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    outputSource: sort_merged_overlapped_from_a_and_b/sorted_file
    label: "BED file with merged overlapped intervals for experiments A and B"
    doc: "BED file with merged overlapped intervals for experiments A and B"
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Merged overlapped from A and B"
        displayMode: "COLLAPSE"
        height: 40

  collected_statistics:
    type: File
    format: "http://edamontology.org/format_3835"
    outputSource: collect_statistics/output_file
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'
    label: "Overlap statistics"
    doc: "Overlap statistics"


steps:

  dedup_and_sort_a:
    run: ../tools/custom-bash.cwl
    in:
      input_file: intervals_file_a
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n > unique_from_a.bed
    out:
    - output_file

  dedup_and_sort_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: intervals_file_b
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n > unique_from_b.bed
    out:
    - output_file

  get_unique_from_a:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: dedup_and_sort_a/output_file
      file_b: dedup_and_sort_b/output_file
      no_overlaps:
        default: true
    out:
    - intersected_file

  sort_unique_from_a:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: get_unique_from_a/intersected_file
      key:
        default: ["1,1","2,2n","3,3n"]
    out:
    - sorted_file

  convert_to_xls_unique_from_a:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sort_unique_from_a/sorted_file
      script:
        default: >
          cat $0 | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t0\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  annotate_unique_from_a:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: convert_to_xls_unique_from_a/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
      output_filename:
        default: "annotated_unique_from_a.tsv"
    out:
    - result_file

  groom_unique_from_a:
    run: ../tools/custom-bash.cwl
    in:
      input_file: annotate_unique_from_a/result_file
      script:
        default: >
          cat $0 | cut -f 1-9,15 > `basename $0`
    out:
    - output_file

  get_unique_from_b:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: dedup_and_sort_b/output_file
      file_b: dedup_and_sort_a/output_file
      no_overlaps:
        default: true
    out:
    - intersected_file

  sort_unique_from_b:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: get_unique_from_b/intersected_file
      key:
        default: ["1,1","2,2n","3,3n"]
    out:
    - sorted_file

  convert_to_xls_unique_from_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sort_unique_from_b/sorted_file
      script:
        default: >
          cat $0 | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t0\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  annotate_unique_from_b:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: convert_to_xls_unique_from_b/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
      output_filename:
        default: "annotated_unique_from_b.tsv"
    out:
    - result_file

  groom_unique_from_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: annotate_unique_from_b/result_file
      script:
        default: >
          cat $0 | cut -f 1-9,15 > `basename $0`
    out:
    - output_file

  get_overlapped_from_a:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: dedup_and_sort_a/output_file
      file_b: dedup_and_sort_b/output_file
      report_from_a_once:
        default: true
    out:
    - intersected_file

  sort_overlapped_from_a:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: get_overlapped_from_a/intersected_file
      key:
        default: ["1,1","2,2n","3,3n"]
      output_filename:
        default: "overlapped_from_a.bed"
    out:
    - sorted_file

  convert_to_xls_overlapped_from_a:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sort_overlapped_from_a/sorted_file
      script:
        default: >
          cat $0 | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t0\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  annotate_overlapped_from_a:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: convert_to_xls_overlapped_from_a/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
      output_filename:
        default: "annotated_overlapped_from_a.tsv"
    out:
    - result_file

  groom_overlapped_from_a:
    run: ../tools/custom-bash.cwl
    in:
      input_file: annotate_overlapped_from_a/result_file
      script:
        default: >
          cat $0 | cut -f 1-9,15 > `basename $0`
    out:
    - output_file

  get_overlapped_from_b:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: dedup_and_sort_b/output_file
      file_b: dedup_and_sort_a/output_file
      report_from_a_once:
        default: true
    out:
    - intersected_file

  sort_overlapped_from_b:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: get_overlapped_from_b/intersected_file
      key:
        default: ["1,1","2,2n","3,3n"]
      output_filename:
        default: "overlapped_from_b.bed"
    out:
    - sorted_file

  convert_to_xls_overlapped_from_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sort_overlapped_from_b/sorted_file
      script:
        default: >
          cat $0 | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t0\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  annotate_overlapped_from_b:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: convert_to_xls_overlapped_from_b/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
      output_filename:
        default: "annotated_overlapped_from_b.tsv"
    out:
    - result_file

  groom_overlapped_from_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: annotate_overlapped_from_b/result_file
      script:
        default: >
          cat $0 | cut -f 1-9,15 > `basename $0`
    out:
    - output_file

  combine_overlapped_from_a_and_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: [get_overlapped_from_a/intersected_file, get_overlapped_from_b/intersected_file]
      script:
        default: |
          cat "$0" > temp.tsv
          cat "$1" >> temp.tsv
          cat temp.tsv | sort -u -k1,1 -k2,2n -k3,3n > merged_overlapped_from_a_and_b.bed
          rm temp.tsv
    out:
    - output_file

  merge_overlapped_from_a_and_b:
    run: ../tools/bedtools-merge.cwl
    in:
      bed_file: combine_overlapped_from_a_and_b/output_file
    out:
    - merged_bed_file

  sort_merged_overlapped_from_a_and_b:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: merge_overlapped_from_a_and_b/merged_bed_file
      key:
        default: ["1,1","2,2n","3,3n"]
    out:
    - sorted_file

  convert_to_xls_merged_overlapped_from_a_and_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sort_merged_overlapped_from_a_and_b/sorted_file
      script:
        default: >
          cat $0 | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t0\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  annotate_merged_overlapped_from_a_and_b:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: convert_to_xls_merged_overlapped_from_a_and_b/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
      output_filename:
        default: "annotated_merged_overlapped_from_a_and_b.tsv"
    out:
    - result_file

  groom_merged_overlapped_from_a_and_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: annotate_merged_overlapped_from_a_and_b/result_file
      script:
        default: >
          cat $0 | cut -f 1-9,15 > `basename $0`
    out:
    - output_file

  collect_statistics:
    run: ../tools/custom-bash.cwl
    in:
      input_file:
      - groom_unique_from_a/output_file
      - groom_unique_from_b/output_file
      - groom_overlapped_from_a/output_file
      - groom_overlapped_from_b/output_file
      - groom_merged_overlapped_from_a_and_b/output_file
      script:
        default: |
          UNIQUE_A=$(($(cat $0 | wc -l)-1))
          UNIQUE_B=$(($(cat $1 | wc -l)-1))
          OVERLAPPED_A=$(($(cat $2 | wc -l)-1))
          OVERLAPPED_B=$(($(cat $3 | wc -l)-1))
          MERGED_OVERLAPPED_A_B=$(($(cat $4 | wc -l)-1))
          echo "| Intervals                      | Counts                    |" >> statistics.md
          echo "| -------------------------------| ------------------------- |" >> statistics.md
          echo "| Unique from A                  | ${UNIQUE_A}               |" >> statistics.md
          echo "| Unique from B                  | ${UNIQUE_B}               |" >> statistics.md
          echo "| Overlapped from A              | ${OVERLAPPED_A}           |" >> statistics.md
          echo "| Overlapped from B              | ${OVERLAPPED_B}           |" >> statistics.md
          echo "| Merged overlapped from A and B | ${MERGED_OVERLAPPED_A_B}  |" >> statistics.md
    out:
    - output_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Pairwise genomic regions intersection"
label: "Pairwise genomic regions intersection"
s:alternateName: "Pairwise genomic regions intersection"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/peak-intersect.cwl
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
  Pairwise genomic regions intersection
  =============================================

  Overlaps peaks from two ChIP/ATAC experiments