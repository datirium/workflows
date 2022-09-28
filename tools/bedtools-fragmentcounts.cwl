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
      # inputs
      bam="$0"; scale="$1"; prefix="$2"
      # clean bam/bed
      bedtools bamtobed -i $bam -bedpe > mapped_pairs.bed
      awk '$1==$4 && $6-$2 < 1000 {print $0}' mapped_pairs.bed | grep -v "^\." > mapped_pairs.clean.bed
      cut -f 1,2,6 mapped_pairs.clean.bed | sort -k1,1 -k2,2n -k3,3n  > mapped_pairs.fragments.bed
      sort mapped_pairs.fragments.bed | uniq -c | awk -F'\t' '{
        split($1,count_chr," ")
        if($3>$2){
          printf("%s\t%s\t%s\t%s\n",count_chr[2],$2,$3,count_chr[1])
        }else{
          printf("%s\t%s\t%s\t%s\n",count_chr[2],$3,$2,count_chr[1])
        }
      }' > $prefix.fragmentcounts.bed
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