cwlVersion: v1.0
class: Workflow
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  bam_file:
    type: File[]
    label: BAM files
    doc: Array of input BAM files
  output_folder:
    type: string[]
    label: BAM file names
    doc: Array of names for output folders
  fragment_size:
    type: int[]
    label: Fragment sizes
    doc: Array of fragment sizes
  total_reads:
    type: int[]
    label: Total reads numbers
    doc: Array of total reads number for downstream normalization
outputs:
  tag_folder:
    type: Directory[]
    label: Tag directories
    doc: Array of tag directories
    outputSource: make_tag_directory/output_tag_folder
steps:
  make_tag_directory:
    run: ../tools/homer-make-tag-directory.cwl
    in:
      bam_file: bam_file
      output_folder: output_folder
      fragment_size: fragment_size
      total_reads: total_reads
    scatter:
    - bam_file
    - output_folder
    - fragment_size
    - total_reads
    scatterMethod: dotproduct
    out:
    - output_tag_folder
doc: |
  Workflow runs homer-make-tag-directory.cwl tool using scatter for the following inputs
    - bam_file
    - fragment_size
    - total_reads

  `dotproduct` is used as a `scatterMethod`, so one element will be taken from each array to construct each job:
    1) bam_file[0] fragment_size[0] total_reads[0]
    2) bam_file[1] fragment_size[1] total_reads[1]
       ...
    N) bam_file[N] fragment_size[N] total_reads[N]

  `bam_file`, `fragment_size` and `total_reads` arrays should have the identical order.
sd:version: 100
label: heatmap-prepare
