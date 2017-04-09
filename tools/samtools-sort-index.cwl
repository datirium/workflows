#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InitialWorkDirRequirement
  listing:
    - $(inputs.sort_input)
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() {
      if (inputs.csi && !inputs.bai){
        return '.csi';
      } else {
        return '.bai';
      }
    };
  - var default_bam = function() {
      if (inputs.trigger == true){
        return inputs.sort_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".sorted.bam";
      } else {
        return inputs.sort_input.location.split('/').slice(-1)[0];
      }
    };


hints:
- class: DockerRequirement
  dockerPull: scidap/samtools:v1.4
  dockerFile: >
    $import: ./dockerfiles/samtools-Dockerfile

inputs:

  bash_script_sort:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = True ]
      then
        echo "Run: samtools sort " ${@:1}
        samtools sort "${@:1}"
      else
        echo "Skip samtools sort " ${@:1}
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run samtools sort with all input parameters or skip it if trigger is false

  bash_script_index:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = True ]
      then
        echo "Run: samtools index " ${@:1}
        samtools index "${@:1}"
      else
        echo "Skip samtools index " ${@:1}
      fi
    inputBinding:
      position: 20
    doc: |
      Bash function to run samtools index with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    doc: |
      If true - run samtools, if false - return sort_input, previously staged into output directory

  sort_compression_level:
    type: int?
    inputBinding:
      position: 11
      prefix: -l
    doc: |
      SORT: desired compression level for the final output file, ranging from 0 (uncompressed)
      or 1 (fastest but minimal compression) to 9 (best compression but slowest to write),
      similarly to gzip(1)'s compression level setting.
      If -l is not used, the default compression level will apply.

  sort_by_name:
    type: boolean?
    inputBinding:
      position: 14
      prefix: -n
    doc: |
      Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates

  sort_output_filename:
    type: string?
    inputBinding:
      position: 12
      prefix: -o
      valueFrom: |
        ${
            if (self == null || inputs.trigger == false){
              return default_bam();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Write the final sorted output to FILE, rather than to standard output.
      Only out.bam|out.cram

  threads:
    type: int?
    doc: |
      Set number of sorting and compression threads [1] (Only for sorting)

  sort_input:
    type: File
    inputBinding:
      position: 16
      valueFrom: $(self.basename)           # Why do we need this
    doc: |
      Input only in.sam|in.bam|in.cram

  interval:
    type: int?
    inputBinding:
      position: 24
      prefix: -m
    doc: |
      Set minimum interval size for CSI indices to 2^INT [14]

  csi:
    type: boolean?
    doc: |
      Generate CSI-format index for BAM files. If input isn't cram.

  bai:
    type: boolean?
    doc: |
      Generate BAI-format index for BAM files [default]. If input isn't cram.

outputs:
  bam_bai_pair:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.sort_output_filename == null || inputs.trigger == false){
              return default_bam();
            } else {
              return inputs.sort_output_filename;
            }
        }
    secondaryFiles:
      ${
          if (inputs.sort_input.secondaryFiles && inputs.trigger == false){
            return inputs.sort_input.secondaryFiles;
          } else {
            return self.location + ext();
          }
        }

baseCommand: [bash, '-c']

arguments:
#   run script sort position 5
  - valueFrom: $(inputs.trigger)
    position: 6
  # -l - position 11
  # -o sort_output_filename - position 12
  - valueFrom: bam
    position: 13
    prefix: -O
    # -n - position 14
  - valueFrom: $(inputs.threads?inputs.threads:1)
    position: 15
    prefix: -@
  # sort_input - position 16
  - valueFrom: ";"
    position: 17
    shellQuote: false

  - valueFrom: "bash"
    position: 18
  - valueFrom: "-c"
    position: 19
#   run script index position 20
  - valueFrom: $(inputs.trigger)
    position: 21
  - valueFrom: $(inputs.bai?'-b':inputs.csi?'-c':[])
    position: 23
    # -m - position 24
  - valueFrom: $(inputs.threads?inputs.threads:1)
    position: 25
    prefix: -@
  - valueFrom: |
      ${
          if (inputs.sort_output_filename == null || inputs.trigger == false){
            return default_bam();
          } else {
            return inputs.sort_output_filename;
          }
      }
    position: 26
  - valueFrom: |
      ${
          if (inputs.sort_output_filename == null || inputs.trigger == false){
            return default_bam() + ext();
          } else {
            return inputs.sort_output_filename + ext();
          }
      }
    position: 27