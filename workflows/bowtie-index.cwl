cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:

  genome:
    type: string
    label: "Genome"
    doc: "Used by BioWardrobe to set genome"

  fasta_input_file:
    type: File
    label: "FASTA input file"
    format: "http://edamontology.org/format_1929"
    doc: "Reference genome input FASTA file"

outputs:
  indices_folder:
    type: Directory
    label: "Bowtie indices folder"
    doc: "Folder which includes all Bowtie generated indices files"
    outputSource: files_to_folder/folder

  annotation_file:
    type: File?
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"

  genome_size:
    type: string?
    label: "Effective genome size"
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

  chrom_length:
    type: File?
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosome length file"

steps:
  bowtie_generate_indices:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: fasta_input_file
      index_base_name: genome
    out: [indices]

  files_to_folder:
    run: ../expressiontools/files-to-folder.cwl
    in:
      input_files: bowtie_generate_indices/indices
    out: [folder]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "bowtie-index"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/dummy/bowtie-index.cwl
s:codeRepository: https://github.com/SciDAP/workflows
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
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681

doc: >
  Workflow makes indices for bowtie v1.2.0 (12/30/2016).

s:about: |
  Workflow makes indices for [bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016).

  It performs the following steps:
  1. Executes `bowtie-index` to generate indices requires genome [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) file as input, returns results as a group of main and secondary files
  2. Transforms results from the previous step into Direcotry data type