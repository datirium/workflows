cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/hopach:v0.0.4


inputs:

  genelist_files:
    type: File[]
    inputBinding:
      prefix: "--input"
    doc: "Genelist TSV or CSV files with RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm columns"

  legend_names:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--name"
    doc: "Genelist file aliases to make the legend for generated plots. Order corresponds to the genelist files. Default: basename of genelist files"

  target_column:
    type: string?
    inputBinding:
      prefix: "--target"
    doc: "Column name from genelist files to be used by Hopach. Default: Rpkm"

  group_by:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--combine"
    doc: "Column name list from genelist files to be used for grouping. Default: RefseqId GeneId Chrom TxStart TxEnd Strand"

  dist_metric:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor"]
    inputBinding:
      prefix: "--dist"
    doc: "Distance metric. Default: cosangle"

  logtransform:
    type: boolean?
    inputBinding:
      prefix: "--logtransform"
    doc: "Log2 transform input data prior running hopach. Default: false"

  keep_discarded:
    type: boolean?
    inputBinding:
      prefix: "--keep"
    doc: "Keep discarded by --min parameter rows at the end of the file. Default: false"

  export_heatmap:
    type: boolean?
    inputBinding:
      prefix: "--heatmap"
    doc: "Export heatmap to png. Default: false"

  export_distance_matrix:
    type: boolean?
    inputBinding:
      prefix: "--distmatrix"
    doc: "Export distance matrix to png. Default: false"

  export_variability_plot:
    type: boolean?
    inputBinding:
      prefix: "--variability"
    doc: "Export clsuter variability plot to png. Default: false"

  threshold:
    type: float?
    inputBinding:
      prefix: "--min"
    doc: "Min value for target column. Default: 0"
  
  palette:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--palette"
    doc: "Color list for custom palette. Default: black red yellow"

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: "Output file prefix"


outputs:


  ordered_genelist:
    type: File
    outputBinding:
      glob: "*_raw_data.tsv"
    doc: "Combined genelist file ordered by hopach clustering results"

  distance_matrix_png:
    type: File?
    outputBinding:
      glob: "*_dist_matrix.png"
    doc: "Distance matrix ordered by hopach clustering results"

  heatmap_png:
    type: File?
    outputBinding:
      glob: "*_heatmap.png"
    doc: "Heatmap ordered by hopach clustering results"

  variability_plot_png:
    type: File?
    outputBinding:
      glob: "*_variablility.png"
    doc: "Cluster variability plot"

  stderr_log:
    type: File
    outputBinding:
      glob: "hopach_stderr.log"
    doc: "Hopach stderr log"

  stdout_log:
    type: File
    outputBinding:
      glob: "hopach_stdout.log"
    doc: "Hopach stdout log"


baseCommand: ["hopach_order.R"]
stderr: hopach_stderr.log
stdout: hopach_stdout.log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/hopach-metadata.yaml

label: "HOPACH - Hierarchical Ordered Partitioning and Collapsing Hybrid"
s:name: "HOPACH - Hierarchical Ordered Partitioning and Collapsing Hybrid"
s:alternateName: |
  The HOPACH clustering algorithm builds a hierarchical tree of clusters by recursively
  partitioning a data set,while ordering and possibly collapsing clusters at each level.

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/hopach.cwl
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Runs hopach clustering algorithm with the combined by specific columns genelist files.
  Return ordered combined genelist file and optional distance matrix, cluster variability and heatmap plots.
  Works with minimum two genelist files.

s:about: |
  usage: hopach_order.R [-h] --input INPUT [INPUT ...] [--name NAME [NAME ...]]
                        [--target TARGET] [--combine COMBINE [COMBINE ...]]
                        [--dist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
                        [--logtransform] [--keep] [--heatmap] [--distmatrix]
                        [--variability] [--min MIN]
                        [--palette PALETTE [PALETTE ...]] [--output OUTPUT]

  Hopach Ordering

  optional arguments:
    -h, --help            show this help message and exit
    --input INPUT [INPUT ...]
                          Input CSV/TSV genelist files
    --name NAME [NAME ...]
                          Names, the order corresponds to input. Default:
                          basename of --input files
    --target TARGET       Target column name to be used by Hopach. Default: Rpkm
    --combine COMBINE [COMBINE ...]
                          Combine inputs by columns names. Default: RefseqId,
                          GeneId, Chrom, TxStart, TxEnd, Strand
    --dist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric. Default: cosangle
    --logtransform        Log2 transform input data prior running hopach.
                          Default: false
    --keep                Keep discarded values at the end of the file. Default:
                          false
    --heatmap             Export heatmap to png. Default: false
    --distmatrix          Export distance matrix to png. Default: false
    --variability         Export clsuter variability plot to png. Default: false
    --min MIN             Min value for target column value. Default: 0
    --palette PALETTE [PALETTE ...]
                          Palette color names. Default: black, red, yellow
    --output OUTPUT       Output prefix. Default: ordered_genelist

