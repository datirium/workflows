cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-deseq:v0.0.23


inputs:

  isoforms_file:
    type: File
    inputBinding:
      position: 5
      prefix: "--isoforms"
    doc: "Isoforms CSV file"

  genes_filename:
    type: string?
    inputBinding:
      position: 6
      prefix: "--gene"
    doc: "Output TSV gene expression filename"

  common_tss_filename:
    type: string?
    inputBinding:
      position: 7
      prefix: "--tss"
    doc: "Output TSV common tss expression filename"


outputs:

  genes_file:
    type: File
    outputBinding:
      glob: $(inputs.genes_filename?inputs.genes_filename:"*genes.tsv")
    doc: "Output TSV gene expression file"

  common_tss_file:
    type: File
    outputBinding:
      glob: $(inputs.common_tss_file?inputs.common_tss_file:"*common_tss.tsv")
    doc: "Output TSV common tss expression file"


baseCommand: ["get_gene_n_tss.R"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/deseq-metadata.yaml

s:name: "group-isoforms"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/group-isoforms.cwl
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
  Tool runs get_gene_n_tss.R script to group isoforms by gene and common TSS

s:about: |
  usage: get_gene_n_tss.R [-h] --isoforms ISOFORMS [--gene GENE] [--tss TSS]

  Group isoform expression data by gene and common TSS

  optional arguments:
    -h, --help           show this help message and exit
    --isoforms ISOFORMS  Input CSV isoform expression file
    --gene GENE          Output TSV gene expression file
    --tss TSS            Output TSV common tss expression file