#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: scidap/parse-fastqc-results:v0.0.1
  dockerFile: >
    $import: ./dockerfiles/parse-fastqc-results-Dockerfile

inputs:
  input_directory:
    type: Directory
    inputBinding:
      position: 1
    doc: >
      Path to folder which includes RSEM *.isoforms.results and *.genes.results

  prefix_name:
    type: string
    inputBinding:
      position: 2
    default: "./"
    doc: >
      Prefix name to export script results

baseCommand: [perl, /usr/local/bin/gen_matrix.pl]

outputs:
  isoforms_tpm_matrix:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.prefix_name === "."){
            return "./" + "Isoforms_TPM_matrix.txt";
          } else {
            return inputs.prefix_name + "Isoforms_TPM_matrix.txt";
          }
        }
  isoforms_counts_matrix:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.prefix_name === "."){
            return "./" + "Isoforms_Counts_matrix.txt";
          } else {
            return inputs.prefix_name + "Isoforms_Counts_matrix.txt";
          }
        }
  genes_tpm_matrix:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.prefix_name === "."){
            return "./" + "Genes_TPM_Matrix.txt";
          } else {
            return inputs.prefix_name + "Genes_TPM_Matrix.txt";
          }
        }
  genes_counts_matrix:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.prefix_name === "."){
            return "./" + "Genes_Counts_Matrix.txt";
          } else {
            return inputs.prefix_name + "Genes_Counts_Matrix.txt";
          }
        }