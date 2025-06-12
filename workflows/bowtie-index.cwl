cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
sd:metadata:
- ../metadata/indices-header.cwl
inputs:
  genome:
    type: string
    label: Genome type
    doc: Genome type, such as mm10, hg19, hg38, etc
  fasta_file:
    type: File
    format: http://edamontology.org/format_1929
    label: Reference genome FASTA file
    doc: Reference genome FASTA file. Includes all chromosomes
outputs:
  indices_folder:
    type: Directory
    label: Bowtie indices
    doc: Bowtie generated indices folder
    outputSource: bowtie_build/indices_folder
  stdout_log:
    type: File
    label: Bowtie stdout log
    doc: Bowtie generated stdout log
    outputSource: bowtie_build/stdout_log
  stderr_log:
    type: File
    label: Bowtie stderr log
    doc: Bowtie generated stderr log
    outputSource: bowtie_build/stderr_log
steps:
  bowtie_build:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: fasta_file
      index_base_name: genome
    out:
    - indices_folder
    - stdout_log
    - stderr_log
label: Build Bowtie indices
doc: |-
  Workflow runs [Bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016) to build indices for reference
  genome provided in a single FASTA file as fasta_file input. Generated indices are saved in a folder with the name that
  corresponds to the input genome
sd:version: 100
