cwlVersion: v1.0
class: CommandLineTool
label: "Generate Counts Barplot"
doc: |
  This tool generates a bar plot from markdown reports of read counts.

requirements:
  - class: DockerRequirement
    dockerPull: "yourdockerimage:latest"  # Replace with the actual image name
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.read_counts_gct)  # Stage markdown files into the working directory

baseCommand: ["Rscript", "/usr/local/bin/generate_barplot.R"]

inputs:
  read_counts_gct:
    type: File[]
    inputBinding: {}
    doc: "Markdown files containing read count statistics"

  filenames:
    type: string[]
    doc: "Sample aliases corresponding to the markdown files"

outputs:
  reads_barplot:
    type: File
    outputBinding:
      glob: "alignment_stats_barchart.png"
    doc: "Generated bar plot of read counts"

arguments:
  - position: 1
    prefix: "--markdown_files"
    valueFrom: $(inputs.read_counts_gct.map(function(file) { return file.basename; }))
    separate: true
  - position: 2
    prefix: "--aliases"
    valueFrom: $(inputs.filenames)
    separate: true
  - position: 3
    prefix: "--output"
    valueFrom: "alignment_stats_barchart.png"
    separate: true