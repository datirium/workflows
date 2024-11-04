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
      ORIGINAL_FASTQ=$0
      TRIMMED_FASTQ=$1
      TRIMMING_REPORT=$2
      MIN_COUNTS=$3
      TRIMMED_COUNTS=$(( `cat $TRIMMED_FASTQ | wc -l` / 4 )) 
      
      echo "ORIGINAL_FASTQ: ${ORIGINAL_FASTQ}"
      echo "TRIMMED_FASTQ: ${TRIMMED_FASTQ}"
      echo "TRIMMING_REPORT: ${TRIMMING_REPORT}"
      echo "TRIMMED_COUNTS: ${TRIMMED_COUNTS}"

      if (( $TRIMMED_COUNTS < $MIN_COUNTS )); then
        echo "Bypassing adapter trimming"
        cp $ORIGINAL_FASTQ `basename $TRIMMED_FASTQ`
      else
        echo "Using adapter trimming results"
        cp $TRIMMED_FASTQ .
        cp $TRIMMING_REPORT .
      fi
    inputBinding:
      position: 5

  original_fastq_file:
    type: File
    inputBinding:
      position: 6

  trimmed_fastq_file:
    type: File
    inputBinding:
      position: 7

  trimming_report_file:
    type: File
    inputBinding:
      position: 8

  min_reads_count:
    type: int?
    default: 100000
    inputBinding:
      position: 9


outputs:

  selected_fastq_file:
    type: File
    outputBinding:
      glob: $(inputs.trimmed_fastq_file.basename)
      
  selected_report_file:
    type: File?
    outputBinding:
      glob: $(inputs.trimming_report_file.basename)


baseCommand: ["bash", "-c"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "bypass-trimgalore-se"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/bypass-trimgalore-se.cwl
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
  If the number of reads in the trimmed_fastq_file is less then min_reads_count, tool
  will return original_fastq_file and null as selected_report_file. Otherwise, the
  trimmed_fastq_file and trimming_report_file will be returned. Might be usefull in
  case of trimgalore removed all reads from the original_fastq_file

s:about: |
  If the number of reads in the trimmed_fastq_file is less then min_reads_count, tool
  will return original_fastq_file and null as selected_report_file. Otherwise, the
  trimmed_fastq_file and trimming_report_file will be returned. Might be usefull in
  case of trimgalore removed all reads from the original_fastq_file