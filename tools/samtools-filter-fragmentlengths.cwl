cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for samtools-filter-fragmentlengths.cwl tool:\n\n"
      bam=$0; filter=$1
      # filter based on user input fragment length type selection
      if [[ "$filter" == "Default_Range" ]]; then
        samtools view -h $bam | awk 'substr($0,1,1)=="@" || ($9>=0 && $9<=1000) || ($9<=0 && $9>=-1000)' | samtools view -b > filtered.bam
      elif [[ "$filter" == "Histone_Binding_Library" ]]; then
        samtools view -h $bam | awk 'substr($0,1,1)=="@" || ($9>=130 && $9<=300) || ($9<=-130 && $9>=-300)' | samtools view -b > filtered.bam
      elif [[ "$filter" == "Transcription_Factor_Binding_Library" ]]; then
        samtools view -h $bam | awk 'substr($0,1,1)=="@" || ($9>=0 && $9<=130) || ($9<=0 && $9>=-130)' | samtools view -b > filtered.bam
      fi
    inputBinding:
        position: 4

  bam_file:
    type: File
    label: "Input BAM file"
    inputBinding:
      position: 5
    doc: "Input BAM file, does not have to be coordinates sorted"

  fragment_length_filter:
    type: string
    inputBinding:
      position: 13
    doc: |
      Fragment length filter type, retains fragments in ranges.
        Default_Range <1000 bp
        Histone_Binding_Library range 130-300 bp
        Transcription_Factor_Binding_Library range <130 bp


outputs:

  log_file_stdout:
    type: File
    outputBinding:
      glob: "filter_fragment_lengths.log.stdout"
    doc: |
      log for stdout

  log_file_stderr:
    type: File
    outputBinding:
      glob: "filter_fragment_lengths.log.stderr"
    doc: |
      log for stderr

  filtered_bam:
    type: File
    outputBinding:
      glob: "filtered.bam"


baseCommand: ["bash", "-c"]
stdout: 'filter_fragment_lengths.log.stdout'
stderr: 'filter_fragment_lengths.log.stderr'


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "samtools-filter-fragmentlengths"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/samtools-filter-fragmentlengths.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

doc: |
    Tool processes BAM file and returns only reads with insert length based on fragment length filter input.
