cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  genome_indices:
    - "genome-indices.cwl"
  peaklist_A_samples:
    - "cutandrun-macs2-pe.cwl"
    - "cutandrun-seacr-pe.cwl"
    - "trim-atacseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "diffbind.cwl"
    - "filter-peaks-for-heatmap.cwl"
  peaklist_B_samples:
    - "cutandrun-macs2-pe.cwl"
    - "cutandrun-seacr-pe.cwl"
    - "trim-atacseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "diffbind.cwl"
    - "filter-peaks-for-heatmap.cwl"

inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  reference_genome:
    type: Directory?
    'sd:upstreamSource': "genome_indices/bowtie_indices"
    label: "IGV Genome:"
    doc: "Genome index for IGV to use for visualization."

  annotation_file:
    type: File
    'sd:upstreamSource': "genome_indices/annotation"
    label: "Annotation file:"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file for igv"

  set_operator:
    type:
    - "null"
    - type: enum
      name: "Set operation user selection"
      symbols:
      - Intersection
      - Union
      - Difference
      - Complement
    label: "Select set operation"
    sd:preview:
      position: 3
    doc: "Intersection: Returns genomic regions common to all input lists.\n\n
      Union: Merges overlapping and/or adjoining regions from all lists, returns contiguous regions.\n\n
      Difference: Returns genomic regions found within list A, excluding regions from list group B.\n\n
      Complement: Returns genomic regions in the gaps between the contiguous ranges from all lists (any region not covered by a peak(s))."
    'sd:localLabel': true

  peak_list_A:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Peak list A sample:"
    doc: "Choose a sample to use for peak list A, from the ChIP, ATAC, C&R, or diffbind workflows."
    'sd:upstreamSource': "peaklist_A_samples/iaintersect_result"
    'sd:localLabel': true
    sd:preview:
      position: 8

  peak_list_B_group:
    type: File[]
    format: "http://edamontology.org/format_3475"
    label: "Peak list B sample(s):"
    doc: "Select 1 or more samples for peak list group B, from the ChIP, ATAC, C&R, or diffbind workflows."
    'sd:upstreamSource': "peaklist_B_samples/iaintersect_result"
    'sd:localLabel': true
    sd:preview:
      position: 9

  promoter_dist:
    type: int?
    default: 1000
    'sd:layout':
      advanced: true
    label: "Max distance from gene TSS for promoter region assignment:"
    doc: "Max distance from gene TSS (in both directions) for peak to be assigned to the promoter region."

  upstream_dist:
    type: int?
    default: 20000
    'sd:layout':
      advanced: true
    label: "Max distance from the promoter (only in 5' direction) for peak to be assigned to the upstream region:"
    doc: "Max distance from the promoter (only in 5' direction) for peak to be assigned to the upstream region."


outputs:

  overview_tab:
    type: File?
    label: "Markdown formatted log"
    format: "http://edamontology.org/format_3835"
    doc: "Markdown formatted log"
    outputSource: set_operation/overview_file
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  set_operation_peaks_for_igv:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Set peaks from operator, in simple bed format."
    doc: "Regions of interest formatted as headerless BED file with [chrom start end peak_<n>]"
    outputSource: set_operation/filtered_set_for_igv
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "Set Operation Result"
        displayMode: "COLLAPSE"
        height: 40

  peak_list_A_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "input list A bed file, in simple bed format for IGV."
    doc: "headerless BED file with [chrom start end] for IGV"
    outputSource: set_operation/list_A_bed_file_for_igv
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "Input List A Peaks"
        displayMode: "COLLAPSE"
        height: 40

  peak_list_B_group_bed:
    type: File[]
    format: "http://edamontology.org/format_3003"
    label: "input list B array of bed files, in simple bed format for IGV."
    doc: "headerless BED files with [chrom start end] for IGV"
    outputSource: set_operation/list_B_bed_array_for_igv
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "Input List B Group Peaks"
        displayMode: "COLLAPSE"
        height: 40

  iaintersect_input:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Set peaks from operator, formatted for input into iaintersect to find nearest gene per peak."
    doc: "Regions of interest formatted as headered BED file with [chrom start end length abssummit pileup log1-p foldenrich log10q region]. NOTE: pileup, log1-p, foldenrich, and log10q are padded with 0s due to non-carrythrough after set operations."
    outputSource: set_operation/filtered_set_for_iaintersect

  annotated_peaks_file:
    type: File?
    format: "http://edamontology.org/format_3475"
    label: "gene annotated set operated peaks file"
    doc: "nearest gene annotation per peak [refseq_id gene_id txStart txEnd strand chrom start end length abssummit pileup log1-p foldenrich log10q region]"
    outputSource: island_intersect/result_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Annotated Peak Set Results'
        Title: 'set operated peaks with nearest gene annotation'

  filtered_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filtered overlapped peaks in headerless bed format for homer motif analysis"
    doc: "Regions of interest formatted as headerless BED file with [chrom start end name score strand]"
    outputSource: formatting_bed_for_homer/headerless_bed

  filtering_stdout_log_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Filtering stdout log"
    doc: "Filtering stdout log"
    outputSource: set_operation/log_file_stdout

  filtering_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Filtering stderr log"
    doc: "Filtering stderr log"
    outputSource: set_operation/log_file_stderr


steps:

  set_operation:
    run: ../tools/filter-peaks-by-overlap.cwl
    in:
      peak_list_A: peak_list_A
      peak_list_B_group: peak_list_B_group
      set_operator:
        source: set_operator
        valueFrom: $(self)
    out:
      - list_A_bed_file_for_igv
      - list_B_bed_array_for_igv
      - overview_file
      - filtered_set_for_igv
      - filtered_set_for_iaintersect
      - log_file_stdout
      - log_file_stderr

  island_intersect:
    label: "Peak annotation"
    doc: |
      Assigns nearest genes to peaks to explore the biological implication of the open
      chromatin binding sites.
    run: ../tools/iaintersect.cwl
    in:
      input_filename: set_operation/filtered_set_for_iaintersect
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out: [result_file, log_file]

  formatting_bed_for_homer:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: ScatterFeatureRequirement
      - class: ShellCommandRequirement
      inputs:
        script:
          type: string?
          default: |
            # format for homer [chrom start end name score strand]
            awk -F'\t' '{printf("%s\t%.0f\t%.0f\t%s\t%s\t%s\n",$6,$3,$4,$2,"NA",$5)}' <(tail -n+2 $0) | sort | uniq > output-for-homer.tsv
          inputBinding:
            position: 1
        input_file:
          type: File
          inputBinding:
            position: 2
      outputs:
        headerless_bed:
          type: File
          outputBinding:
            glob: output-for-homer.tsv
      baseCommand: ["bash", "-c"]
    in:
      input_file: island_intersect/result_file
    out:
    - headerless_bed


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Set Operations for Called Peaks (ChIP/ATAC/C&R/diffbind)"
label: "Set Operations for Called Peaks (ChIP/ATAC/C&R/diffbind)"
s:alternateName: "Set Operations for Called Peaks (ChIP/ATAC/C&R/diffbind)"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/filter-peaks-by-overlap.cwl
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
  # Set Operations for Peaks

  This workflow takes as input multiple peak list TSV files (the `iaintersect_result.tsv` output under the
  "Files" output tab) from the ChIP, ATAC, C&R, or diffbind workflows and performs the user-selected set
  operation on the group. Set operations include intersection, union, difference, and complement. See the
  tooltip for the `set_operator` input for more details.