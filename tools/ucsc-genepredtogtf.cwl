cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/ucscuserapps:v358_2
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ var root = inputs.annotation_tsv_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.annotation_tsv_file.basename+".gtf":root+".gtf"; } else { return inputs.output_filename; } };
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      TSV_FILE=$0
      GTF_FILE=$1
      cut -f 2- $TSV_FILE | grep -v "exonCount" | genePredToGtf file stdin $GTF_FILE
    inputBinding:
      position: 5
    doc: |
      Bash function to run cut -f 2- in.gp | genePredToGtf file stdin out.gp
  annotation_tsv_file:
    type: File
    inputBinding:
      position: 6
    doc: TSV annotation file from UCSC Genome Browser
  output_filename:
    type: string?
    default: ''
    inputBinding:
      valueFrom: $(default_output_filename())
      position: 7
    doc: Output file name
outputs:
  annotation_gtf_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: GTF annotation file
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bash
- -c
stdout: genepredtogtf_stdout.log
stderr: genepredtogtf_stderr.log
doc: |
  genePredToGtf - Convert genePred table or file to gtf
label: ucsc-genepredtogtf
