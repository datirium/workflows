cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { return inputs.bam_file.location.split('/').slice(-1)[0]; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
inputs:
  bash_script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        samtools rmdup "${@:1}"
      else
        echo "Skip samtools rmdup " ${@:1}
        SOURCE=${@:(-2):1}
        TARGET=${@:(-1):1}
        cp $SOURCE $TARGET
        if [ -f $SOURCE.bai ]; then
          cp $SOURCE.bai $TARGET.bai
        fi
      fi
    inputBinding:
      position: 1
    doc: |
      Bash function to run samtools rmdup with all input parameters or skip it if trigger is false
  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 2
      valueFrom: |
        ${ return self ? "true" : "false" }
    doc: |
      If true - run samtools rmdup, if false - return input bam_file and optional index file
      Use valueFrom to return string instead of boolean, because if return boolean False, argument is not printed
  single_end:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 3
      prefix: -s
    doc: |
      rmdup for SE reads
  force_single_end:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 4
      prefix: -S
    doc: |
      treat PE reads as SE in rmdup (force -s)
  bam_file:
    type: File
    inputBinding:
      position: 10
    doc: |
      Input sorted bam file (index file is optional)
  output_filename:
    type:
    - 'null'
    - string
    inputBinding:
      position: 11
      valueFrom: |
        ${
            if (self == "" || inputs.trigger == false){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: ''
    doc: |
      Writes the output bam file to output_filename if set,
      otherwise generates output_filename on the base of bam_file
outputs:
  rmdup_output:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == "" || inputs.trigger == false){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }
    secondaryFiles: |
      ${
          if (inputs.bam_file.secondaryFiles && inputs.trigger == false){
            return inputs.bam_file.secondaryFiles;
          } else {
            return "null";
          }
        }
    doc: File with removed duplicates or bam_file with optional secondaryFiles
  rmdup_log:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == "" || inputs.trigger == false){
            return default_output_filename().split('.').slice(0,-1).join('.') + '.rmdup';
          } else {
            return inputs.output_filename.split('.').slice(0,-1).join('.') + '.rmdup';
          }
        }
baseCommand:
- bash
- -c
arguments:
- valueFrom: |
    ${
      if (inputs.output_filename == "" || inputs.trigger == false){
        return " > " + default_output_filename().split('.').slice(0,-1).join('.') + ".rmdup 2>&1";
      } else {
        return " > " + inputs.output_filename.split('.').slice(0,-1).join('.') + ".rmdup 2>&1";
      }
    }
  position: 100000
  shellQuote: false
doc: "Tool to remove duplicates from coordinate sorted BAM file set as input `bam_file`.\nIf input `trigger` is set to `true` or isn't set at all (`true` is used by default), run `samtools rmdup`, return\nnewly generated BAM file without duplicates and log as outputs `rmdup_output` and `rmdup_log`.\nIf input `trigger` is set to `false`, return unchanged BAM and index files (if provided in secondaryFiles).\n\nTrigger logic is implemented in bash script set by default in input `bash_script`. If first argment $0 (which is `trigger` input)\nis true, run `samtools rmdup` with the rest of the arguments. If $0 is not true, skip samtools running and copy\ninto output directory BAM file. Input BAM file and the filename which we are waiting to be the\noutput of tool are always set as two last arguments of the script.\n\nInput `trigger` is Boolean, but returns String, because of `valueFrom` field. The `valueFrom` is used, because if `trigger`\nis false, cwl-runner doesn't append this argument at all to the the `baseCommand` - new feature of CWL v1.0.2. Alternatively,\n`prefix` field could be used, but it causes changing logic in bash script saved in `bash_script` input.\n\n`default_output_filename` function is used for generating output filename if input `output_filename` is not set or in\ncase when `trigger` is false and we need to return original BAM and index files.\n\nOutput `rmdup_output` returns `secondaryFiles` only in case when `trigger` was set to false (we need to rerun index).\n  \n"
label: samtools-rmdup
