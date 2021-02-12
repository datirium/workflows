cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            return inputs.bed_file.basename;
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  bed_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: |
      The input BED file

  chrom_length_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-g"
    doc: |
      Input genome file with chromosome lengths

  bi_direction:
    type:
      - "null"
      - int
      - float
    inputBinding:
      position: 7
      prefix: "-b"
    doc: |
      Increase the BED entry by the same number base pairs
      or a fraction of the feature's length in each direction

  left_direction:
    type:
      - "null"
      - int
      - float
    inputBinding:
      position: 8
      prefix: "-l"
    doc: |
      The number of base pairs or a fraction of the feature's length
      to subtract from the start coordinate

  right_direction:
    type:
      - "null"
      - int
      - float
    inputBinding:
      position: 9
      prefix: "-r"
    doc: |
      The number of base pairs or a fraction of the feature's length
      to add to the end coordinate

  strand_based:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "-s"
    doc: |
      Define -l and -r based on strand
      E.g. if used, -l 500 for a negative-stranded feature, 
      it will add 500 bp downstream
  
  percent_based:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "-pct"
    doc: |
      Define -l, -b and -r as a fraction of the feature's length

  save_header:
    type: boolean?
    inputBinding:
      position: 12
      prefix: "-header"
    doc: |
      Print the header from the input file prior to results

  output_filename:
    type: string?
    default: ""
    doc: "Output file name"


outputs:

  extended_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Extended BED file"


baseCommand: ["bedtools", "slop"]
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/bedtools-metadata.yaml

s:name: "bedtools-slop"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/bedtools-slop.cwl
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
  Increases the size of each feature in a feature file by a user-defined number of bases.
  If not using -b, then -l and -r should be used together

  
s:about: |
  Usage:   bedtools slop [OPTIONS] -i <bed/gff/vcf> -g <genome> [-b <int> or (-l and -r)]

  Options: 
    -b	Increase the BED/GFF/VCF entry -b base pairs in each direction.
      - (Integer) or (Float, e.g. 0.1) if used with -pct.

    -l	The number of base pairs to subtract from the start coordinate.
      - (Integer) or (Float, e.g. 0.1) if used with -pct.

    -r	The number of base pairs to add to the end coordinate.
      - (Integer) or (Float, e.g. 0.1) if used with -pct.

    -s	Define -l and -r based on strand.
      E.g. if used, -l 500 for a negative-stranded feature, 
      it will add 500 bp downstream.  Default = false.

    -pct	Define -l and -r as a fraction of the feature's length.
      E.g. if used on a 1000bp feature, -l 0.50, 
      will add 500 bp "upstream".  Default = false.

    -header	Print the header from the input file prior to results.