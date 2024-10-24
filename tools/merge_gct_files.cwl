cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: biowardrobe2/gctmerge:v0.0.1
  - class: InlineJavascriptRequirement

inputs:
  gct_files:
    type: File[]
    doc: "List of GCT files to merge"

  output_filename:
    type: string
    default: "merged_counts.gct"
    inputBinding:
      position: 2
      prefix: "-o"
    doc: "Name of the output merged GCT file"

arguments:
  - position: 1
    valueFrom: |
      ${ return ["-i"].concat(inputs.gct_files.map(function(f){ return f.path; })); }
    shellQuote: false

outputs:
  merged_gct_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

baseCommand: [ "Rscript", "/usr/src/app/merge_gct_files.R" ]

stdout: stdout.log
stderr: stderr.log