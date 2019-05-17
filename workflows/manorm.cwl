cwlVersion: v1.0
class: Workflow


inputs:

  peak_file_first:
    type: File
  peak_file_second:
    type: File
  bam_file_first:
    type: File
  bam_file_second:
    type: File
  fragment_size_first:
    type: int
  fragment_size_second:
    type: int


outputs:

  common_peak_merged_file:
    type: File
    outputSource: manorm/common_peak_merged_file

steps:

  manorm:
    in:
      peak_file_first: peak_file_first
      peak_file_second: peak_file_second
      bam_file_first: bam_file_first
      bam_file_second: bam_file_second
      fragment_size_first: fragment_size_first
      fragment_size_second: fragment_size_second
    out: [common_peak_merged_file]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/manorm:v0.0.1
      inputs:
        peak_file_first:
          type: File
          inputBinding:
            position: 5
        peak_file_second:
          type: File
          inputBinding:
            position: 6
        bam_file_first:
          type: File
          inputBinding:
            position: 7
        bam_file_second:
          type: File
          inputBinding:
            position: 9
        fragment_size_first:
          type: int
          inputBinding:
            position: 10
        fragment_size_second:
          type: int
          inputBinding:
            position: 11
      outputs:
        common_peak_merged_file:
          type: File
          outputBinding:
            glob: "MAnorm_result_commonPeak_merged.xls"
      baseCommand: ["run_manorm.sh"]
