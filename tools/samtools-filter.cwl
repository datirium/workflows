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
      echo "Copy $0 to temp.bam and create new index"
      cp $0 temp.bam
      samtools index temp.bam
      echo "Filtering BAM file"
      echo "samtools idxstats temp.bam | cut -f 1 | grep -v -E \"`echo $1 | sed -e 's/ /$|/g'`$|\*\" | xargs samtools view -q $2 -F $3 -o $4 temp.bam"
      samtools idxstats temp.bam | cut -f 1 | grep -v -E "`echo $1 | sed -e 's/ /$|/g'`$|\*" | xargs samtools view -q $2 -F $3 -o $4 temp.bam
      echo "Sorting BAM file"
      echo "samtools sort $4 -o $4"
      samtools sort $4 -o $4
      echo "Indexing BAM file"
      echo "samtools index $4"
      samtools index $4
      rm -f temp.bam temp.bam.bai
    inputBinding:
      position: 5
    doc: "Script to exclude chromosomes from the BAM file and filter reads by quality"

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

  quality:
    type: int?
    inputBinding:
      position: 8
    default: 0
    doc: "Skip alignments with MAPQ smaller than INT. Default 0"
  
  negative_flag:
    type: int?
    inputBinding:
      position: 9
    default: 0
    doc: "Do not output alignments with any bits set in INT present in the FLAG field. Default 0"

  output_filename:
    type: string?
    inputBinding:
      position: 10
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
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

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
  Excludes chromosomes from the input BAM file. Filters reads by quality.
  If there is only one chromosome present, you cannot exclude it

s:about: |
  Excludes chromosomes from the input BAM file