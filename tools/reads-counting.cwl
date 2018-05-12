cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: scidap/reads-counting:v0.0.1
  dockerFile: >
    $import: ./dockerfiles/reads-counting-Dockerfile

inputs:
  aligned:
    type: File
    inputBinding:
      prefix: -in=
      separate: false
      position: 1
    doc: |
      Input bam file

  annotation:
    type: File
    inputBinding:
      prefix: -annotation=
      separate: false
      position: 2
    doc: |
      Input tab delimited annotation file

  rpkm-cutoff:
    type: double?
    inputBinding:
      prefix: -rpkm-cutoff=
      separate: false
      position: 3
    doc: |
      Set cutoff value for rpkm, below which everything will be changed to -rpkm-cutoff-val

  rpkm-cutoff-val:
    type: double?
    inputBinding:
      prefix: -rpkm-cutoff-val=
      separate: false
      position: 4
    doc: |
      Set minimum value for rpkm which will be set instead of all the values below -rpkm-cutoff

  math-converging:
    type:
      name: "average"
      type: enum
      symbols: ["arithmetic","geometric"]
    inputBinding:
      prefix: -math-converging=
      separate: false
      position: 5
    default: arithmetic
    doc: |
      Type of converging for density array processing: arithmetic, geometric

  threads:
    type: double?
    inputBinding:
      prefix: -threads=
      separate: false
      position: 6
    default: 1
    doc: |
      Number of threads to run the program

  log:
    type: string?
    inputBinding:
      prefix: -log=
      separate: false
      position: 7
    default: "output.log"
    doc: |
      Path to output log file

  rnaSeqType:
    type:
      name: "seqType"
      type: enum
      symbols: ["dUTP","RNA"]
    inputBinding:
      prefix: -rna_seq=
      separate: false
      position: 8
    doc: |
      Type of the analysis


  outputFilename:
    type: string?
    inputBinding:
      prefix: -out=
      separate: false
      position: 9
    default: rpkm.csv
    doc: |
      Path to the output rpkm file

  ignore_chrom:
    type: string?
    inputBinding:
      prefix: -sam_ignorechr=
      separate: false
      position: 10
    doc: |
      Ignore chromosome (in a case of using RNA spike-in)


outputs:
  rpkmFile:
    type: File
    outputBinding:
      glob: $(inputs.outputFilename)

  logFile:
    type: File
    outputBinding:
      glob: $(inputs.log)

baseCommand: [reads-counting]


