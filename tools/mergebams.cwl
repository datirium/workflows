cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/mergebams:v0.0.1


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        echo "Skip running mergebams"
      else
        echo "Running mergebams ${@:1}"
        mergebams "${@:1}"
      fi
    inputBinding:
      position: 1
    doc: "Bash script to optionally skip running mergebams"

  skip:
    type: boolean?
    default: false
    inputBinding:
      position: 2
      valueFrom: $(self?"true":"false")      # need to have it as a string, otherwise false is not propogated.
    doc: "If true, skip running the tool"

  possorted_genome_bam_bai:
    type: File[]
    secondaryFiles:
    - .bai
    inputBinding:
      position: 3
      prefix: "--inputs"
      itemSeparator: ","
    doc: |
      Indexed RNA BAM file containing
      position-sorted reads aligned to
      the genome and transcriptome, as
      well as unaligned reads.

  output_filename:
    type: string
    inputBinding:
      position: 4
      prefix: "--output"
    doc: |
      Output name for the merged
      coordinate sorted BAM file.

  threads:
    type: int?
    inputBinding:
      position: 5
      prefix: "--threads"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  merged_possorted_genome_bam_bai:
    type: File?
    outputBinding:
      glob: "*.bam"
    secondaryFiles:
    - .bai
    doc: |
      Merged coordinate sorted
      indexed BAM file.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]


stdout: mergebams_stdout.log
stderr: mergebams_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Merge 10x reads (RNA)"
label: "Merge 10x reads (RNA)"
s:alternateName: "Merges 10x-generated BAM files with RNA reads and updates barcode suffixes"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/mergebams.cwl
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
  Merge 10x reads (RNA)

  Merges 10x-generated BAM files with RNA reads and
  updates barcode suffixes to reflect the order of
  the input BAM files similar to Cell Ranger Aggregate.


s:about: |
  USAGE:
      mergebams [OPTIONS] --inputs <inputs> --output <output>

  FLAGS:
      -h, --help       Prints help information
      -V, --version    Prints version information

  OPTIONS:
          --inputs <inputs>      Input bam files, comma-separated list.
          --output <output>      Name for the merged BAM file.
          --threads <threads>    Number of threads to be used when sorting and indexing.