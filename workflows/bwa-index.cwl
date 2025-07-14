cwlVersion: v1.0
class: Workflow
requirements:
- class: StepInputExpressionRequirement
sd:upstream:
  genome_indices: genome-indices.cwl
inputs:
  alias:
    type: string
    label: Sample short name/Alias
    sd:preview:
      position: 1
  reference_fasta:
    type: File
    sd:upstreamSource: genome_indices/fasta_output
    label: 'Reference genome FASTA file to index:'
    sd:localLabel: true
    doc: |
      FASTA file of the reference genome that will be indexed.
    sd:preview:
      position: 2
outputs:
  index_directory:
    type: Directory
    label: Directory containing the original FASTA, faidx, dict, and bwa index files.
    outputSource: index_reference/bwa_index
  log_file_stdout:
    type: File
    format: http://edamontology.org/format_2330
    label: stdout logfile
    outputSource: index_reference/log_file_stdout
    sd:visualPlugins:
    - markdownView:
        tab: Overview
  log_file_stderr:
    type: File
    format: http://edamontology.org/format_2330
    label: stderr logfile
    outputSource: index_reference/log_file_stderr
steps:
  index_reference:
    run: ../tools/bwa-index.cwl
    in:
      ref_genome_fasta: reference_fasta
    out:
    - bwa_index
    - log_file_stdout
    - log_file_stderr
label: BWA index pipeline
doc: |
  This workflow indexes the input reference FASTA with bwa, and generates faidx and dict file using samtools.
  This index sample can then be used as input into the germline variant calling workflow, or others that may
  include this workflow as an upstream source.

  ### __Inputs__
   - FASTA file of the reference genome that will be indexed.

  ### __Outputs__
   - Directory containing the original FASTA, faidx, dict, and bwa index files.
   - stdout log file (output in Overview tab as well)
   - stderr log file

  ### __Data Analysis Steps__
  1. cwl calls dockercontainer robertplayer/scidap-gatk4 to index reference FASTA with bwa, and generates faidx and dict files using samtools

  ### __References__
    - Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 25(14), 1754–1760.
sd:version: 100
