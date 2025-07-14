cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/morpheus:v0.0.2
inputs:
  read_counts_gct:
    type: File
    inputBinding:
      prefix: --gct
    doc: |
      Path to the input GCT file.
  output_prefix:
    type: string?
    inputBinding:
      prefix: --output
    doc: |
      Output prefix for generated files
outputs:
  heatmap_html:
    type: File
    outputBinding:
      glob: '*.html'
    doc: |
      Morpheus heatmap in HTML format
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- run_morpheus.R
stdout: morpheus_stdout.log
stderr: morpheus_stderr.log
label: Morpheus Heatmap
doc: |
  Morpheus Heatmap

  Generates Morpheus heatmap from input GCT file
