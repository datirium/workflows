cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for bash script in seacr.cwl tool:\n\n"
      printf "INPUTS:\n"
      echo "\$0 - $0"
      echo "\$1 - $1"
      echo "\$2 - $2"
      echo "\$3 - $3"
      echo "\$4 - $4"
      # inputs
      bam="$0"; scale="$1"; prefix="$2"; genome="$3"; filter="$4"
      # generate paired end (pe) bed from bam
      bedtools bamtobed -i $bam -bedpe 2> /dev/null > mapped_pairs.bed
      # filter based on user input fragment length type selection
      if [[ "$filter" == "Default_Range" ]]; then
        awk '$1==$4 && $6-$2 < 1000 {print $0}' mapped_pairs.bed | grep -v "^\." > mapped_pairs.clean.bed
      elif [[ "$filter" == "Histone_Binding_Library" ]]; then
        awk '$1==$4 && $6-$2 > 130 && $6-$2 < 300 {print $0}' mapped_pairs.bed | grep -v "^\." > mapped_pairs.clean.bed
      elif [[ "$filter" == "Transcription_Factor_Binding_Library" ]]; then
        awk '$1==$4 && $6-$2 < 130 {print $0}' mapped_pairs.bed | grep -v "^\." > mapped_pairs.clean.bed
      fi
      # sort and count per base coverage per genome
      cut -f 1,2,6 mapped_pairs.clean.bed | sort -k1,1 -k2,2n -k3,3n  > mapped_pairs.fragments.bed
      bedtools genomecov -bg -i mapped_pairs.fragments.bed -g $genome > $prefix.fragmentcounts.bed
      # apply scale for normalization to spike-in (or ecoli mapped read count by default)
      awk -F'\t' -v scale=$scale '{
        printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4*scale)
        }' $prefix.fragmentcounts.bed > $prefix.fragmentcounts_scaled.bed
    inputBinding:
        position: 1

  bam_file:
    type: File?
    inputBinding:
      position: 2
    doc: |
      A bam formatted sequence read alignment file from a cut
      and run/tag experiment.

  scale:
    type: float?
    inputBinding:
      position: 3
    doc: "Coefficient to scale the genome coverage by a constant factor"

  output_prefix:
    type: string
    inputBinding:
      position: 9
    doc: |
      Basename of input file that SEACR will use to name the
      output tsv file: <output_prefix>.<peakcalling_mode>.bed

  chrom_length_file:
    type: File?
    inputBinding:
      position: 11
    doc: |
      A genome file, where col1 are chromosome names and col2
      contains an integer length (bp) of the chromosome.

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

  sorted_bed:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.fragmentcounts.bed')
    doc: |
      SEACR prepped bed file formatted from PE bam file.

  sorted_bed_scaled:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.fragmentcounts_scaled.bed')
    doc: |
      Scaled SEACR prepped bed file formatted from PE bam file.

  log_file_stderr:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.fragmentcounts.stderr')
    doc: |
      log for BASH stderr

  log_file_stdout:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.fragmentcounts.stdout')
    doc: |
      log for BASH stdout


baseCommand: ["bash", "-c"]
stderr: $(inputs.output_prefix + '.fragmentcounts.stderr')
stdout: $(inputs.output_prefix + '.fragmentcounts.stdout')


doc: |
  Tool runs bedtools and custom BASH commands to process and 
  format a paired end bam file into a sorted bed file ready
  for peak calling by SEACR. Also outputs a spike-in scaled
  bed file for normalized peak calling by SEACR.

  Reference protocol, STEP 14:
  https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=14