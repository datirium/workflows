cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for samtools-filter-fragmentlengths.cwl tool:\n\n"
      bam=$0; filter=$1
      # filter based on user input fragment length type selection
      if [[ "$filter" == "default_below_1000" ]]; then
        samtools view -h $bam | awk 'substr($0,1,1)=="@" || ($9>=0 && $9<=1000) || ($9<=0 && $9>=-1000)' | samtools view -b > filtered.bam
      elif [[ "$filter" == "histones_130_to_300" ]]; then
        samtools view -h $bam | awk 'substr($0,1,1)=="@" || ($9>=130 && $9<=300) || ($9<=-130 && $9>=-300)' | samtools view -b > filtered.bam
      elif [[ "$filter" == "TF_below_130" ]]; then
        samtools view -h $bam | awk 'substr($0,1,1)=="@" || ($9>=0 && $9<=130) || ($9<=0 && $9>=-130)' | samtools view -b > filtered.bam
      fi
    inputBinding:
      position: 4
  bam_file:
    type: File
    label: Input BAM file
    inputBinding:
      position: 5
    doc: Input BAM file, does not have to be coordinates sorted
  fragment_length_filter:
    type: string
    inputBinding:
      position: 13
    doc: |
      Fragment length filter type, retains fragments in ranges.
        Default_Range <1000 bp
        Histone_Binding_Library range 130-300 bp
        Transcription_Factor_Binding_Library range <130 bp
outputs:
  log_file_stdout:
    type: File
    outputBinding:
      glob: filter_fragment_lengths.log.stdout
    doc: |
      log for stdout
  log_file_stderr:
    type: File
    outputBinding:
      glob: filter_fragment_lengths.log.stderr
    doc: |
      log for stderr
  filtered_bam:
    type: File
    outputBinding:
      glob: filtered.bam
baseCommand:
- bash
- -c
stdout: filter_fragment_lengths.log.stdout
stderr: filter_fragment_lengths.log.stderr
doc: |
  Tool processes BAM file and returns only reads with insert length based on fragment length filter input.
label: samtools-filter-fragmentlengths
