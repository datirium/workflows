cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.bam_bai_pair,
                  "entryname": inputs.bam_bai_pair.basename,
                  "writable": true
                }
              ]
    }
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.11
inputs:
  script:
    type: string?
    default: |2

      #!/bin/bash
      echo "Rename $0 to temp.bam"
      mv $0 temp.bam
      if [ -f $0.bai ]; then
        echo "Rename $0.bai to temp.bam.bai"
        mv $0.bai temp.bam.bai
      fi

      echo "Sorting BAM file by name"
      echo "samtools sort -n -@ $3 -o namesorted.bam temp.bam"
      samtools sort -n -@ $3 -o namesorted.bam temp.bam

      echo "Filling in mate coordinates and inserting size fields"
      echo "samtools fixmate -m -@ $3 namesorted.bam fixed.bam"
      samtools fixmate -m -@ $3 namesorted.bam fixed.bam

      echo "Sorting BAM file by coordinates"
      echo "samtools sort -@ $3 -o positionsorted.bam fixed.bam"
      samtools sort -@ $3 -o positionsorted.bam fixed.bam

      if [ "$1" = "true" ]
      then
        echo "Only marking PCR duplicates"
        echo "samtools markdup -c -s -@ $3 positionsorted.bam markduped.bam"
        samtools markdup -c -s -@ $3 positionsorted.bam markduped.bam 2> markdup_report.tsv
      else
        echo "Removing PCR duplicates"
        echo "samtools markdup -c -r -s -@ $3 positionsorted.bam markduped.bam"
        samtools markdup -c -r -s -@ $3 positionsorted.bam markduped.bam 2> markdup_report.tsv
      fi

      echo "Sorting BAM file"
      echo "samtools sort -@ $3 markduped.bam -o $2"
      samtools sort -@ $3 markduped.bam -o $2

      echo "Indexing BAM file"
      echo "samtools index $2"
      samtools index $2

      echo "Removing temporary files"
      rm -f namesorted.bam fixed.bam positionsorted.bam markduped.bam temp.bam*
    inputBinding:
      position: 5
    doc: Script to remove PCR duplicates
  bam_bai_pair:
    type: File
    inputBinding:
      position: 6
    doc: BAM (optionally BAI) files
  keep_duplicates:
    type: boolean?
    default: false
    inputBinding:
      position: 7
      valueFrom: $(self?"true":"false")
    doc: |
      If true duplicates will be only
      marked, oterwise - removed
  output_filename:
    type: string?
    inputBinding:
      position: 8
      valueFrom: |
        ${
          return (self == "")?inputs.bam_bai_pair.basename:self;
        }
    default: ''
    doc: Output filename for the filtered BAM file
  threads:
    type: int?
    inputBinding:
      position: 9
    default: 1
    doc: Number of threads to use
outputs:
  deduplicated_bam_bai_pair:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles:
    - .bai
    doc: BAM+BAI files with PCR duplicates removed
  markdup_report:
    type: File
    outputBinding:
      glob: markdup_report.tsv
    doc: Markdup report
baseCommand:
- bash
- -c
doc: |
  Removes or only marks PCR duplicates from
  coordinate sorted and indexed BAM file.
  Returns coordinate sorted and indexed BAM
  files. Stages input bam_bai_pair to workdir.
  Otherwise samtools sort fails.
label: samtools-markdup
