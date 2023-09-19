cwlVersion: v1.0
class: Workflow


requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


inputs:

  bam_file:
    type: File[]
    label: "BAM files"
    doc: "Array of input BAM files"

  output_folder:
    type: string[]
    label: "BAM file names"
    doc: "Array of names for output folders"

  fragment_size:
    type: int[]
    label: "Fragment sizes"
    doc: "Array of fragment sizes"

  total_reads:
    type: int[]
    label: "Total reads numbers"
    doc: "Array of total reads number for downstream normalization"


outputs:

  tag_folder:
    type: Directory[]
    label: "Tag directories"
    doc: "Array of tag directories"
    outputSource: make_tag_directory/output_tag_folder


steps:

  make_tag_directory:
    run: ../tools/homer-make-tag-directory.cwl
    in:
      bam_file: bam_file
      output_folder: output_folder
      fragment_size: fragment_size
      total_reads: total_reads
    scatter:
      - bam_file
      - output_folder
      - fragment_size
      - total_reads
    scatterMethod: dotproduct
    out: [output_tag_folder]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "heatmap-prepare"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/heatmap-prepare.cwl
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
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681

doc: |
  Workflow runs homer-make-tag-directory.cwl tool using scatter for the following inputs
    - bam_file
    - fragment_size
    - total_reads

  `dotproduct` is used as a `scatterMethod`, so one element will be taken from each array to construct each job:
    1) bam_file[0] fragment_size[0] total_reads[0]
    2) bam_file[1] fragment_size[1] total_reads[1]
       ...
    N) bam_file[N] fragment_size[N] total_reads[N]

  `bam_file`, `fragment_size` and `total_reads` arrays should have the identical order.

s:about: |
  Runs homer-make-tag-directory.cwl with the scatter using dot product of the inputs bam_file, fragment_size and total_reads.
