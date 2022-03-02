cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat-export:v0.0.2


inputs:

  seurat_data_rds:
    type: File
    inputBinding:
      prefix: "--rds"
    doc: |
      Path to the RDS file to load Seurat object from. RDS file can be
      produced by run_seurat.R, assign_cell_types.R, or run_seurat_wnn.R
      scripts.

  source_column:
    type: string
    inputBinding:
      prefix: "--source"
    doc: |
      Column from metadata to extract clusters information from

  splitby:
    type: string?
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from metadata to split output files by. Barcode suffixes that define
      the origin of the cells will be removed. Values from this column will be
      appended to the names of output files.
      Default: do not split, export combined file

  reductions:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--reductions"
    doc: |
      List of reductions to extract metadata from. Each reduction
      name will be appended to the output files names.
      Default: umap, pca, atac_lsi, rnaumap, atacumap, wnnumap

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./metadata

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  memory:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  clusters_metadata:
    type: File[]
    outputBinding:
      glob: "*_clusters.tsv"

  reductions_metadata:
    type: File[]
    outputBinding:
      glob: "*[!_clusters].tsv"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["extract_seurat_metadata.R"]

stdout: extract_metadata_stdout.log
stderr: extract_metadata_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Extract Seurat Metadata"
s:name: "Single-cell Extract Seurat Metadata"
s:alternateName: "Extracts metadata from Seurat RDS file"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-extract-metadata.cwl
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
  Single-cell Extract Seurat Metadata
  ===================================
  Extracts metadata from Seurat RDS file


s:about: |
  usage: extract_seurat_metadata.R
        [-h] --rds RDS --source SOURCE [--splitby SPLITBY]
        [--reductions [REDUCTIONS ...]] [--output OUTPUT] [--verbose]
        [--cpus CPUS] [--memory MEMORY]

  Extracts cells, UMAP, and clusters from the Seurat RDS file

  optional arguments:
    -h, --help            show this help message and exit
    --rds RDS             Path to the RDS file to load Seurat object from. RDS
                          file can be produced by run_seurat.R,
                          assign_cell_types.R, or run_seurat_wnn.R scripts.
    --source SOURCE       Column from metadata to extract clusters information
                          from.
    --splitby SPLITBY     Column from metadata to split output files by. Barcode
                          suffixes that define the origin of the cells will be
                          removed. Values from this column will be appended to
                          the names of output files. Default: do not split,
                          export combined file
    --reductions [REDUCTIONS ...]
                          List of reductions to extract metadata from. Each
                          reduction name will be appended to the output files
                          names. Default: umap, pca, atac_lsi, rnaumap,
                          atacumap, wnnumap
    --output OUTPUT       Output prefix. Default: ./metadata
    --verbose             Print debug information. Default: false
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
