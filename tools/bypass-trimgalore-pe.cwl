cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ResourceRequirement
    ramMin: 7024
    coresMin: 1


hints:
- class: DockerRequirement
  dockerPull: scidap/scidap:v0.0.4


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      ORIGINAL_FASTQ_1=$0
      TRIMMED_FASTQ_1=$1
      TRIMMING_REPORT_1=$2

      ORIGINAL_FASTQ_2=$3
      TRIMMED_FASTQ_2=$4
      TRIMMING_REPORT_2=$5

      MIN_COUNTS=$6
      TRIMMED_COUNTS=$(( `cat $TRIMMED_FASTQ_1 | wc -l` / 4 )) 
      
      echo "ORIGINAL_FASTQ_1: ${ORIGINAL_FASTQ_1}"
      echo "TRIMMED_FASTQ_1: ${TRIMMED_FASTQ_1}"
      echo "TRIMMING_REPORT_1: ${TRIMMING_REPORT_1}"

      echo "ORIGINAL_FASTQ_2: ${ORIGINAL_FASTQ_2}"
      echo "TRIMMED_FASTQ_2: ${TRIMMED_FASTQ_2}"
      echo "TRIMMING_REPORT_2: ${TRIMMING_REPORT_2}"

      echo "TRIMMED_COUNTS: ${TRIMMED_COUNTS}"

      if (( $TRIMMED_COUNTS < $MIN_COUNTS )); then
        echo "Bypassing adapter trimming"
        cp $ORIGINAL_FASTQ_1 `basename $TRIMMED_FASTQ_1`
        cp $ORIGINAL_FASTQ_2 `basename $TRIMMED_FASTQ_2`
      else
        echo "Using adapter trimming results"
        cp $TRIMMED_FASTQ_1 .
        cp $TRIMMED_FASTQ_2 .
        cp $TRIMMING_REPORT_1 .
        cp $TRIMMING_REPORT_2 .
      fi
    inputBinding:
      position: 5

  original_fastq_file_1:
    type: File
    inputBinding:
      position: 6

  trimmed_fastq_file_1:
    type: File
    inputBinding:
      position: 7

  trimming_report_file_1:
    type: File
    inputBinding:
      position: 8

  original_fastq_file_2:
    type: File
    inputBinding:
      position: 9

  trimmed_fastq_file_2:
    type: File
    inputBinding:
      position: 10

  trimming_report_file_2:
    type: File
    inputBinding:
      position: 11

  min_reads_count:
    type: int?
    default: 100000
    inputBinding:
      position: 12


outputs:

  selected_fastq_file_1:
    type: File
    outputBinding:
      glob: $(inputs.trimmed_fastq_file_1.basename)
      
  selected_report_file_1:
    type: File?
    outputBinding:
      glob: $(inputs.trimming_report_file_1.basename)

  selected_fastq_file_2:
    type: File
    outputBinding:
      glob: $(inputs.trimmed_fastq_file_2.basename)
      
  selected_report_file_2:
    type: File?
    outputBinding:
      glob: $(inputs.trimming_report_file_2.basename)


baseCommand: ["bash", "-c"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "bypass-trimgalore-pe"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/bypass-trimgalore-pe.cwl
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
  If the number of reads in the trimmed_fastq_file_1 is less then min_reads_count, tool
  will return original_fastq_file_1/2 and nulls as selected_report_file_1/2. Otherwise,
  the trimmed_fastq_file_1/2 and trimming_report_file_1/2 will be returned. Might be
  usefull in case of trimgalore removed all reads from the original_fastq_file_1/2.

s:about: |
  If the number of reads in the trimmed_fastq_file_1 is less then min_reads_count, tool
  will return original_fastq_file_1/2 and nulls as selected_report_file_1/2. Otherwise,
  the trimmed_fastq_file_1/2 and trimming_report_file_1/2 will be returned. Might be
  usefull in case of trimgalore removed all reads from the original_fastq_file_1/2.