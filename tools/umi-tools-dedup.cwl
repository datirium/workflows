cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { var ext = inputs.bam_file.basename.split('.').slice(-1)[0]; var root = inputs.bam_file.basename.split('.').slice(0,-1).join('.'); return inputs.output_filename?inputs.output_filename:root+"_dedup."+ext; };
hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/umi_tools:1.0.1--py38h0213d0e_2
inputs:
  bash_script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]; then
        umi_tools dedup --random-seed=12345 "${@:1}"
      else
        echo "Skip umi_tools dedup " ${@:1}
        cp $2 $4
        cp $2.bai $4.bai
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run umi_tools dedup with all input parameters or skip it if trigger is false
  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 6
      valueFrom: $(self?"true":"false")
    doc: |
      If true - run umi_tools dedup, if false - return bam_file input, previously staged into the output directory
  bam_file:
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 7
      prefix: -I
    doc: Input BAM file
  output_filename:
    type: string?
    inputBinding:
      position: 8
      prefix: -S
      valueFrom: $(default_output_filename())
    default: ''
    doc: Output filename
  paired_end:
    type: boolean?
    inputBinding:
      position: 9
      prefix: --paired
    doc: |
      Inputs BAM file is paired end - output both read pairs.
      This will also force the use of the template length to
      determine reads with the same mapping coordinates.
  output_stats:
    type: string?
    inputBinding:
      position: 10
      prefix: --output-stats=
      separate: false
    default: umi_tools_stats
    doc: Specify location to output stats
  multimapping_detection_method:
    type: string?
    inputBinding:
      position: 11
      prefix: --multimapping-detection-method=
      separate: false
    doc: |
      Some aligners identify multimapping using bam tags.
      Setting this option to NH, X0 or XT will use these
      tags when selecting the best read amongst reads with
      the same position and umi [default=none]
outputs:
  dedup_bam_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    secondaryFiles:
    - .bai
  output_stats:
    type:
    - 'null'
    - File[]
    outputBinding:
      glob: $(inputs.output_stats + "*")
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bash
- -c
stdout: umi_tools_dedup_stdout_file.log
stderr: umi_tools_dedup_stderr_file.log
doc: |
  Deduplicate BAM files based on the first mapping co-ordinate and the UMI attached to the read
  Only -I, --paired and -S parameters are implemented.
label: umi-tools-dedup
