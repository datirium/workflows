cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_prefix = function() { var ext = '.'; var root = inputs.bam_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.bam_file.basename+ext:root+ext; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/geep:v0.0.5
inputs:
  bam_file:
    type:
    - File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 10
      prefix: --bam
    doc: |
      Set the path to the coordinate sorted BAM file. Required
  annotation_file:
    type:
    - File
    inputBinding:
      position: 11
      prefix: --annotation
    doc: |
      Set the path to the GTF or TAB-delimited file. Required
  log_filename:
    type: string?
    inputBinding:
      position: 12
      prefix: --log
    doc: |
      Set the path to save LOG file. Default: /dev/null
  output_prefix:
    type: string?
    inputBinding:
      position: 13
      prefix: --output
      valueFrom: |
        ${
            if (self == ""){
              return default_output_prefix();
            } else {
              return self;
            }
        }
    default: ''
    doc: |
      Set the prefix to save output files. Default: ""
  rpkm_threshold:
    type: double?
    inputBinding:
      position: 14
      prefix: --threshold
    doc: |
      Set rpkm cutoff threshold, below which everything will be changed to value set with --cutoff. Default: 0
  rpkm_cutoff:
    type: double?
    inputBinding:
      position: 15
      prefix: --cutoff
    doc: |
      Set rpkm cutoff value to be used for all rpkms below --threshold. Default: 0
  exclude_chr:
    type: string?
    inputBinding:
      position: 16
      prefix: --exclude
    doc: |
      Coma separated list of chromosomes to be ignored. Default: ""
  min_interval_length:
    type: int?
    inputBinding:
      position: 17
      prefix: --minIntLen
    doc: |
      Set the minimal interval length. All shorter intervals will be discarded. Default: 0
  min_read_length:
    type: int?
    inputBinding:
      position: 18
      prefix: --minReadLen
    doc: |
      Set the minimal read length. All parts of spliced reads that intersect with exon in
      less than minReadLen nucleotides will be discarded. Default: 0
  keep_unique:
    type: boolean?
    inputBinding:
      position: 19
      prefix: --keepUnique
    doc: |
      Set this flag if you want prevent distributing the isoform unique reads among other isoforms. Default: False
  dutp:
    type: boolean?
    inputBinding:
      position: 20
      prefix: --dutp
    doc: |
      Set this dutp flag if strand specific analysys should be made. Default: False
  max_cycles:
    type: int?
    inputBinding:
      position: 21
      prefix: --cycles
    doc: |
      Set the maximum number of cycles used for read balancing. Default: 2000
  threads:
    type: int?
    inputBinding:
      position: 22
      prefix: --threads
    doc: |
      Set the number of threads. Default: 1
outputs:
  isoforms_file:
    type: File
    outputBinding:
      glob: '*isoforms.csv'
  genes_file:
    type: File
    outputBinding:
      glob: '*genes.csv'
  raw_file:
    type: File
    outputBinding:
      glob: '*raw.txt'
  log_file:
    type:
    - 'null'
    - File
    outputBinding:
      glob: $(inputs.log_filename)
baseCommand:
- geep
label: geep
doc: |
  Tool calculates RPKM values grouped by isoforms or genes.

  `default_output_prefix` function returns default prefix based on `bam_file` basename, if `output_prefix` is not
  provided.
