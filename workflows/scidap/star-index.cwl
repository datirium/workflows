cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:

  fasta_input_files:
    type: File[]
    label: "FASTA input files"
    format: "http://edamontology.org/format_192"
    doc: "Reference genome input FASTA file(s)"

  annotation_input_file:
    type: File
    label: "GTF input file"
    format: "http://edamontology.org/format_2306"
    doc: "Annotation input file"

  threads:
    type: int?
    label: "Number of threads to run tools"
    doc: "Number of threads for those steps that support multithreading"

outputs:
  indices_folder:
    type: Directory
    label: "STAR indices folder"
    format: "http://edamontology.org/format_2330"
    doc: "Folder which includes all STAR generated indices files"
    outputSource: files_to_folder/folder

  chr_length:
    type: File
    label: "Chromosome lenth file"
    doc: "STAR generated chromosome length file"
    outputSource: get_chr_length_file/selected_file

steps:
  star_generate_indices:
    run: ../../tools/star-genomegenerate.cwl
    in:
      genomeFastaFiles: fasta_input_files
      sjdbGTFfile: annotation_input_file
      threads: threads
    out: [indices]

  files_to_folder:
    run: ../../expressiontools/files-to-folder.cwl
    in:
      input_files: star_generate_indices/indices
    out: [folder]

  get_chr_length_file:
    run: ../../expressiontools/get-file-by-name.cwl
    in:
      input_files: star_generate_indices/indices
      basename_regex:
        default: chrNameLength.txt
    out: [selected_file]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "star-index"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/scidap/star-index.cwl
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

s:about: >
  Current workflow should be used to generate STAR genome indices files. It performs the following steps:
  1. Use STAR to generate genome indices files on the base of input FASTA and GTF files, return results as an array of files
  2. Transform array of file to the Direcotry data type
  3. Get chrNameLength.txt file as a separate output