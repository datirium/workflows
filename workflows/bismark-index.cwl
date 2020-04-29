cwlVersion: v1.0
class: Workflow


'sd:metadata':
  - "../metadata/indices-header.cwl"


inputs:

  genome:
    type: string
    label: "Genome type"
    doc: "Genome type, such as mm10, hg19, hg38, etc"

  fasta_file:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file (*.fa or *.fasta, or *.gz)"
    doc: "Reference genome FASTA file (*.fa or *.fasta, or *.gz). Includes all chromosomes"


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

  fasta_to_folder:
    in:
      genome_fasta: fasta_file
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
- http://schema.org/docs/schema_org_rdfa.html

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

doc:
  $include: ../descriptions/bismark-index.md


# doc: |
#   Copy fasta_file file to the folder and run run bismark_genome_preparation script to prepare indices for Bismark Methylation Analysis.
#   Bowtie2 aligner is used by default. The name of the output indices folder is equal to the genome input.