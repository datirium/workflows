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
    default: |

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
    doc: "Script to remove PCR duplicates"

  bam_bai_pair:
    type: File
    inputBinding:
      position: 6
    doc: BAM (optionally BAI) files

  keep_duplicates:
    type: boolean?
    default: false                        # somehow when omitted, valueFrom is not evaluated
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
    default: ""
    doc: "Output filename for the filtered BAM file"

  threads:
    type: int?
    inputBinding:
      position: 9
    default: 1
    doc: "Number of threads to use"


outputs:

  deduplicated_bam_bai_pair:
    type: File
    outputBinding:
      glob: "*.bam"
    secondaryFiles:
    - .bai
    doc: "BAM+BAI files with PCR duplicates removed"

  markdup_report:
    type: File
    outputBinding:
      glob: "markdup_report.tsv"
    doc: "Markdup report"


baseCommand: [bash, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-markdup"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/samtools-markdup.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Removes or only marks PCR duplicates from
  coordinate sorted and indexed BAM file.
  Returns coordinate sorted and indexed BAM
  files. Stages input bam_bai_pair to workdir.
  Otherwise samtools sort fails.

s:about: |
  Removes or only marks PCR duplicates from
  coordinate sorted and indexed BAM file.
  Returns coordinate sorted and indexed BAM
  files. Stages input bam_bai_pair to workdir.
  Otherwise samtools sort fails.