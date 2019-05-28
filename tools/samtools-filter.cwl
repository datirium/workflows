cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      echo "Filtering BAM file"
      samtools idxstats $0 | cut -f 1 | grep -v -E "`echo $1 | sed -e 's/ /$|/g'`$|\*" | xargs samtools view -o $2 $0
      echo "Sorting BAM file"
      samtools sort $2 -o $2
      echo "Indexing BAM file"
      samtools index $2
    inputBinding:
      position: 5
    doc: Script to exclude chromosomes from the BAM file

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
    doc: "Space separated list of the chromosemes to exclude"

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


outputs:

  filtered_bam_bai_pair:
    type: File
    outputBinding:
      glob: "*.bam"
    secondaryFiles:
      - .bai
    doc: "Filtered BAM+BAI files"


baseCommand: [bash, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-filter"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/samtools-filter.cwl
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
  Excludes chromosomes from the input BAM file

s:about: |
  Excludes chromosomes from the input BAM file
