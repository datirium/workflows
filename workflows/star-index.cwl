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

  annotation_input_file:
    type: File?
    label: "GTF input file"
    format: "http://edamontology.org/format_2306"
    doc: "Annotation input file"

  threads:
    type: int?
    label: "Number of threads to run tools"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true

outputs:

  indices_folder:
    type: Directory
    label: "STAR indices folder"
    doc: "Folder which includes all STAR generated indices files"
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
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    outputSource: get_chr_length_file/selected_file
    doc: "Chromosome length file"

steps:
  star_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_fasta_files: fasta_input_file
      sjdb_gtf_file: annotation_input_file
      threads: threads
    out: [indices]

  files_to_folder:
    run: ../expressiontools/files-to-folder.cwl
    in:
      input_files: star_generate_indices/indices
    out: [folder]

  get_chr_length_file:
    run: ../expressiontools/get-file-by-name.cwl
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

doc: >
  Workflow makes indices for STAR v2.5.3a (03/17/2017).

s:about: |
  Workflow makes indices for [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886).

  It performs the following steps:
  1. Runs `STAR --runMode genomeGenerate` to generate indices, based on [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) and [GTF](http://mblab.wustl.edu/GTF2.html) input files, returns results as an array of files
  2. Transforms array of files into [Direcotry](http://www.commonwl.org/v1.0/CommandLineTool.html#Directory) data type
  3. Separates *chrNameLength.txt* file as an output

