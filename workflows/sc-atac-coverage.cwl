cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


"sd:upstream":
  sc_tools_sample:
  - "sc-multiome-filter.cwl"
  - "sc-atac-reduce.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-ctype-assign.cwl"
  - "sc-rna-azimuth.cwl"
  sc_atac_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  - "cellranger-atac-count.cwl"
  - "cellranger-atac-aggr.cwl"
  genome_indices:
  - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Experiment run through any pipeline related Single-cell ATAC-Seq"
    doc: |
      Path to the RDS file to load Seurat object from. This file
      should include chromatin accessibility information stored
      in the ATAC assay with a proper seqinfo data.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ATAC or RNA+ATAC Sample"
    doc: |
      Any "Cell Ranger ATAC or RNA+ATAC Sample"
      for generating ATAC fragments coverage
      files. This sample can be obtained from
      one of the following pipelines: "Cell
      Ranger Count (RNA+ATAC)", "Cell Ranger
      Aggregate (RNA+ATAC)", "Cell Ranger Count
      (ATAC)", or "Cell Ranger Aggregate (ATAC)".
    "sd:upstreamSource": "sc_atac_sample/atac_fragments_file"
    "sd:localLabel": true

  chrom_length_file:                                                         # not used - need it only for IGV
    type: File
    label: "Genome"
    doc: |
      Reference genome
    "sd:upstreamSource": "genome_indices/chrom_length"
    "sd:localLabel": true

  splitby:
    type: string?
    default: "new.ident"
    label: "Column(s) from the Seurat object metadata to split cells into groups"
    doc: |
      Column from the Seurat object metadata to split cells into groups.
      May be one of the columns added with --metadata or --barcodes
      parameters. Default: split by dataset

  datasets_metadata:
    type: File?
    label: "Optional TSV/CSV file to extend metadata by dataset"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata with
      categorical values using samples identities. First column - 'library_id'
      should correspond to all unique values from the 'new.ident' column of the
      loaded Seurat object. If any of the provided in this file columns are already
      present in the Seurat object metadata, they will be overwritten. When combined
      with --barcodes parameter, first the metadata will be extended, then barcode
      filtering will be applied. Default: no extra metadata is added

  barcodes_data:
    type: File?
    label: "Optional TSV/CSV file to prefilter and extend metadata by barcodes. First column should be named as 'barcode'"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and extend Seurat object
      metadata be selected barcodes. First column should be named as 'barcode'.
      If file includes any other columns they will be added to the Seurat object
      metadata ovewriting the existing ones if those are present.
      Default: all cells used, no extra metadata is added

  flank_distance:
    type: int?
    default: 5
    label: "Distance in bp to flank both start and end of the each fragment in both direction"
    doc: |
      Distance in bp to flank both start and end of the each fragment in both
      direction to generate cut sites coverage. Default: 5
    "sd:layout":
      advanced: true

  export_html_report:
    type: boolean?
    default: true
    label: "Show HTML report"
    doc: |
      Export tehcnical report in HTML format.
      Default: true
    "sd:layout":
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
      - "5"
      - "6"
    default: "4"
    label: "Number of cores/cpus to use"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    "sd:layout":
      advanced: true


outputs:

  peaks_bigbed_file:
    type: File
    outputSource: sc_atac_coverage/peaks_bigbed_file
    label: "Locations of open-chromatin regions"
    doc: |
      Locations of open-chromatin regions ("peaks")
      in bigBed format
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Peaks"
        height: 40

  cut_sites_bigwig_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_coverage/cut_sites_bigwig_file
    label: "Genome coverage for Tn5 cut sites"
    doc: |
      Genome coverage calculated for Tn5 cut sites
      in bigWig format
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Cut sites coverage"
        height: 120

  fragments_bigwig_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_coverage/fragments_bigwig_file
    label: "Genome coverage for ATAC fragments"
    doc: |
      Genome coverage calculated for ATAC fragments
      in bigWig format
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "ATAC fragments coverage"
        height: 120

  sc_report_html_file:
    type: File?
    outputSource: sc_atac_coverage/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  experiment_info:
    type: File
    label: "IGV tracks order"
    doc: |
      Markdown file to explain the tracks order for IGV
    outputSource: create_metadata/output_file
    "sd:visualPlugins":
    - markdownView:
        tab: "Overview"

  sc_atac_coverage_stdout_log:
    type: File
    outputSource: sc_atac_coverage/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the sc_atac_coverage step.

  sc_atac_coverage_stderr_log:
    type: File
    outputSource: sc_atac_coverage/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the sc_atac_coverage step.


steps:

  sc_atac_coverage:
    run: ../tools/sc-atac-coverage.cwl
    in:
      query_data_rds: query_data_rds
      atac_fragments_file: atac_fragments_file
      splitby:
        source: splitby
        valueFrom: $(split_features(self))
      datasets_metadata: datasets_metadata
      barcodes_data: barcodes_data
      flank_distance: flank_distance
      verbose:
        default: true
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      export_html_report: export_html_report
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - peaks_bigbed_file
    - cut_sites_bigwig_file
    - fragments_bigwig_file
    - sc_report_html_file
    - stdout_log
    - stderr_log

  create_metadata:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sc_atac_coverage/fragments_bigwig_file
      script:
        default: |
          #!/bin/bash
          set -- "$0" "$@"
          echo "| Name | Index |" > experiment_info.md
          echo "| :-- | --: |" >> experiment_info.md
          j=1
          for i in "${@}"; do
            echo "| `basename $i` | $j |" >> experiment_info.md
            (( j++ ))
          done;
    out:
    - output_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell ATAC-Seq Genome Coverage"
s:name: "Single-Cell ATAC-Seq Genome Coverage"
s:alternateName: "Generates genome coverage tracks from chromatin accessibility data of selected cells"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-atac-coverage.cwl
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
  Single-Cell ATAC-Seq Genome Coverage

  Generates genome coverage tracks from chromatin
  accessibility data of selected cells