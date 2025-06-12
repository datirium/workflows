cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.11
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      echo "Copy $0 to temp.bam"
      cp $0 temp.bam
      samtools sort temp.bam -o temp_sorted.bam
      samtools index temp_sorted.bam
      echo "Filtering BAM file"
      echo "samtools idxstats temp_sorted.bam | cut -f 1 | grep -v -E \"`echo $1 | sed -e 's/ /$|/g'`$|\*\" | xargs samtools view -q $2 -F $3 -o temp_filtered.bam temp_sorted.bam"
      samtools idxstats temp_sorted.bam | cut -f 1 | grep -v -E "`echo $1 | sed -e 's/ /$|/g'`$|\*" | xargs samtools view -q $2 -F $3 -o temp_filtered.bam temp_sorted.bam
      echo "Sorting BAM file"
      echo "samtools sort temp_filtered.bam -o $4"
      samtools sort temp_filtered.bam -o $4
      echo "Indexing BAM file"
      echo "samtools index $4"
      samtools index $4
      rm -f temp*
    inputBinding:
      position: 5
    doc: Script to exclude chromosomes from the BAM file and filter reads by quality
  bam_bai_pair:
    type: File
    inputBinding:
      position: 6
    secondaryFiles:
    - .bai
    doc: Indexed BAM+BAI files
  exclude_chromosome:
    type: string
    inputBinding:
      position: 7
    doc: Space separated list of the chromosemes to exclude
  quality:
    type: int?
    inputBinding:
      position: 8
    default: 0
    doc: Skip alignments with MAPQ smaller than INT. Default 0
  negative_flag:
    type: int?
    inputBinding:
      position: 9
    default: 0
    doc: Do not output alignments with any bits set in INT present in the FLAG field. Default 0
  output_filename:
    type: string?
    inputBinding:
      position: 10
      valueFrom: |
        ${
          return (self == "")?inputs.bam_bai_pair.basename:self;
        }
    default: ''
    doc: Output filename for the filtered BAM file
outputs:
  filtered_bam_bai_pair:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles:
    - .bai
    doc: Filtered BAM+BAI files
baseCommand:
- bash
- -c
doc: |
  Excludes chromosomes from the input BAM file.
  Optionally filters reads by quality and flags
label: samtools-filter
