cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: ResourceRequirement
  coresMin: 8
  coresMax: 16
  ramMin: 32000
hints:
- class: DockerRequirement
  dockerPull: scidap/clipper:v0.2.0
baseCommand:
- clipper
inputs:
  species:
    type: string
    inputBinding:
      position: 0
      prefix: --species
    doc: 'species: one of ce10 ce11 dm3 hg19 GRCh38 mm9 mm10'
  input_file:
    type: File
    inputBinding:
      position: 1
      prefix: --bam
  gene:
    type: string?
    inputBinding:
      position: 8
      prefix: --gene
  savepickle:
    type: boolean
    default: true
    inputBinding:
      position: 9
      prefix: --save-pickle
  outfile:
    type: string
    default: ''
    inputBinding:
      position: 10
      prefix: --outfile
      valueFrom: |
        ${
          if (inputs.outfile == "") {
            return inputs.input_file.nameroot + ".peakClusters.bed";
          }
          else {
            return inputs.outfile;
          }
        }
outputs:
  output_tsv:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.outfile == "") {
            return inputs.input_file.nameroot + ".peakClusters.bed.tsv";
          }
          else {
            return inputs.outfile + ".tsv";
          }
        }
  output_bed:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.outfile == "") {
            return inputs.input_file.nameroot + ".peakClusters.bed";
          }
          else {
            return inputs.outfile;
          }
        }
  output_pickle:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.outfile == "") {
            return inputs.input_file.nameroot + ".peakClusters.bed.pickle";
          }
          else {
            return inputs.outfile + ".pickle";
          }
        }
doc: |-
  CLIPper is a tool to define peaks in your CLIP-seq dataset.
  CLIPper was developed in the Yeo Lab at the University of California, San Diego.
    Usage: clipper --bam CLIP-seq_reads.srt.bam --species hg19 --outfile CLIP-seq_reads.srt.peaks.bed
