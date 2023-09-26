cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.30


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file
      should include chromatin accessibility information stored
      in the ATAC assay with a proper seqinfo data.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment used in the
      loaded Seurat object. File should be saved in TSV format and to be
      tbi-indexed.

  splitby:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split cells into groups.
      May be one of the columns added with --metadata or --barcodes
      parameters. Default: split by dataset

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
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
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and extend Seurat object
      metadata be selected barcodes. First column should be named as 'barcode'.
      If file includes any other columns they will be added to the Seurat object
      metadata ovewriting the existing ones if those are present.
      Default: all cells used, no extra metadata is added

  flank_distance:
    type: int?
    inputBinding:
      prefix: "--flank"
    doc: |
      Distance in bp to flank both start and end of the each fragment in both
      direction to generate cut sites coverage. Default: 5

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./sc

  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  peaks_bigbed_file:
    type: File
    outputBinding:
      glob: "*_peaks.bigBed"
    doc: |
      Locations of open-chromatin regions ("peaks")
      in bigBed format

  cut_sites_bigwig_file:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cut_cov.bigWig"
    doc: |
      Genome coverage calculated for Tn5 cut sites
      in bigWig format

  fragments_bigwig_file:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_frg_cov.bigWig"
    doc: |
      Genome coverage calculated for fragments
      in bigWig format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_atac_coverage.R"]

stdout: sc_atac_coverage_stdout.log
stderr: sc_atac_coverage_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell ATAC-Seq Genome Coverage"
s:name: "Single-cell ATAC-Seq Genome Coverage"
s:alternateName: "Creates genome coverage bigWig files from the provided fragments file and selected grouping parameters"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-atac-coverage.cwl
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
  Single-cell ATAC-Seq Genome Coverage

  Creates genome coverage bigWig files from the provided fragments file
  and selected grouping parameters.

  --tmpdir parameter is not exposed as input.


s:about: |
  usage: sc_atac_coverage.R
        [-h] --query QUERY --fragments FRAGMENTS [--splitby [SPLITBY ...]]
        [--metadata METADATA] [--barcodes BARCODES] [--flank FLANK] [--verbose]
        [--tmpdir TMPDIR] [--output OUTPUT] [--cpus CPUS] [--memory MEMORY]

  Single-cell ATAC-Seq Genome Coverage

  options:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include chromatin accessibility
                          information stored in the ATAC assay with a proper
                          seqinfo data.
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          used in the loaded Seurat object. File should be saved
                          in TSV format and to be tbi-indexed.
    --splitby [SPLITBY ...]
                          Column from the Seurat object metadata to split cells
                          into groups. May be one of the columns added with
                          --metadata or --barcodes parameters. Default: split by
                          dataset
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - 'library_id' should
                          correspond to all unique values from the 'new.ident'
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          When combined with --barcodes parameter, first the
                          metadata will be extended, then barcode filtering will
                          be applied. Default: no extra metadata is added
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata be selected barcodes.
                          First column should be named as 'barcode'. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --flank FLANK         Distance in bp to flank both start and end of the each
                          fragment in both direction to generate cut sites
                          coverage. Default: 5
    --verbose             Print debug information. Default: false
    --tmpdir TMPDIR       Directory to keep temporary files. Default: either
                          /tmp or defined by environment variables TMPDIR, TMP,
                          TEMP.
    --output OUTPUT       Output prefix. Default: ./sc
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32