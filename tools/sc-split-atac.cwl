cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-split-atac:v0.0.1


inputs:

  atac_fragments_file:
    type: File
    inputBinding:
      prefix: "--fragments"
    doc: |
      Path to GZIP compressed TSV file with ATAC fragments (from Cell Ranger ARC)

  clusters_metadata:
    type: File
    inputBinding:
      prefix: "--clusters"
    doc: |
      Path to headerless TSV file with barcodes (first column) and
      clusters (second column)

  log_level:
    type:
    - "null"
    - type: enum
      symbols:
      - "fatal"
      - "error"
      - "warning"
      - "info"
      - "debug"
    inputBinding:
      prefix: "--loglevel"
    doc: |
      Logging level.
      Default: info

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output file prefix.
      Default: ./split

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of processes to run in parallel.
      Default: 1


outputs:

  atac_fragments_per_cluster_file:
    type: File[]
    outputBinding:
      glob: "*.bed"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_split_atac.py"]

stdout: sc_split_atac_stdout.log
stderr: sc_split_atac_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Split ATAC Fragments"
s:name: "Single-cell Split ATAC Fragments"
s:alternateName: "Splits scATAC fragments produced by Cell Ranger ARC Count/Aggregate pipelines"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-split-atac.cwl
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
  Single-cell Split ATAC Fragments
  =============================================================================
  Splits scATAC fragments produced by Cell Ranger ARC Count/Aggregate pipelines


s:about: |
  usage: sc_split_atac.py [-h] --fragments FRAGMENTS --clusters CLUSTERS [--cpus CPUS] [--loglevel {fatal,error,warning,info,debug}] [--output OUTPUT]

  optional arguments:
    -h, --help            show this help message and exit
    --fragments FRAGMENTS
                          Path to GZIP compressed TSV file with ATAC fragments (from Cell Ranger ARC)
    --clusters CLUSTERS   Path to headerless TSV file with barcodes (first column) and clusters (second column)
    --cpus CPUS           Number of processes to run in parallel
    --loglevel {fatal,error,warning,info,debug}
                          Logging level. Default: info
    --output OUTPUT       Output file prefix
