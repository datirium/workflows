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
    doc: "Folder which include all STAR generated indices files"
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