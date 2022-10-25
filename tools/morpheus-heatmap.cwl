cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/deseq:v0.0.2


inputs:

  read_counts_gct:
    type: File
    inputBinding:
      prefix: "--gct"
    doc: |
      Path to the input GCT file.

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix for generated files


outputs:

  heatmap_html:
    type: File
    outputBinding:
      glob: "*.html"
    doc: |
      Morpheus heatmap in HTML format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [run_morpheus.R]
stdout: morpheus_stdout.log
stderr: morpheus_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "Morpheus Heatmap"
label:  "Morpheus Heatmap"
s:alternateName: "Generates Morpheus heatmap from input GCT file"


s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/morpheus-heatmap.cwl
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
  Morpheus Heatmap

  Generates Morpheus heatmap from input GCT file


s:about: |
  usage: run_morpheus.R
        [-h] --gct GCT [--cluster {row,column,both}] [--output OUTPUT]

  Morpheus heatmap from GCT file

  options:
    -h, --help            show this help message and exit
    --gct GCT             Path to the input GCT file.
    --output OUTPUT       Output prefix for generated files