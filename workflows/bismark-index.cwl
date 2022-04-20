cwlVersion: v1.0
class: Workflow


inputs:

  genome:
    type: string
    label: "Genome type"
    doc: "Genome type, such as mm10, hg19, hg38, etc"

  genome_label:
    type: string?
    label: "Genome label"
    sd:preview:
      position: 1

  genome_description:
    type: string?
    label: "Genome description"
    sd:preview:
      position: 2

  genome_details:
    type: string?
    label: "Genome details"
    sd:preview:
      position: 3

  genome_file:
    type: File
    format: "http://edamontology.org/format_3009"
    label: "Reference genome file (*.2bit, *.fasta, *.fa, *.fa.gz, *.fasta.gz)"
    doc: "Reference genome file (*.2bit, *.fasta, *.fa, *.fa.gz, *.fasta.gz). All chromosomes are included"

  chromosome_list:
    type:
      - "null"
      - string[]
    label: "Chromosome list to be included into the reference genome FASTA file"
    doc: "Filter chromosomes while extracting FASTA from 2bit"


outputs:

  indices_folder:
    type: Directory
    label: "Bismark indices folder"
    doc: "Bismark generated indices folder"
    outputSource: prepare_indices/indices_folder

  stdout_log:
    type: File
    label: "Bismark stdout log"
    doc: "Bismark generated stdout log"
    outputSource: prepare_indices/stdout_log

  stderr_log:
    type: File
    label: "Bismark stderr log"
    doc: "Bismark generated stderr log"
    outputSource: prepare_indices/stderr_log


steps:

  extract_fasta:
    run: ../tools/ucsc-twobit-to-fa.cwl
    in:
      reference_file: genome_file
      chr_list: chromosome_list
    out:
    - fasta_file

  fasta_to_folder:
    in:
      genome_fasta: extract_fasta/fasta_file
      genome_type: genome
    out: [genome_folder]
    run:
      cwlVersion: v1.0
      class: ExpressionTool
      requirements:
      - class: InlineJavascriptRequirement
      inputs:
        genome_fasta: File
        genome_type: string
      outputs:
        genome_folder: Directory
      expression: |
        ${
            return { "genome_folder": {
              "class": "Directory",
              "basename": inputs.genome_type,
              "listing": [inputs.genome_fasta]
            }};
        }

  prepare_indices:
    run: ../tools/bismark-prepare-genome.cwl
    in:
      genome_folder: fasta_to_folder/genome_folder
    out:
    - indices_folder
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Build Bismark indices"
label: "Build Bismark indices"
s:alternateName: "Build indices for Bismark Methylation Pipeline. Bowtie2 aligner is used by default"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/bismark-index.cwl
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

# doc:
#   $include: ../descriptions/bismark-index.md


doc: |
  Copy fasta_file file to the folder and run run bismark_genome_preparation script to prepare indices for Bismark Methylation Analysis.
  Bowtie2 aligner is used by default. The name of the output indices folder is equal to the genome input.