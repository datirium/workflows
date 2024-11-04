cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/ucscuserapps:v358_2


requirements:
  - class: ResourceRequirement
    ramMin: 7620
    coresMin: 1
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: |
      ${
        return  [
                  {
                    "entry": inputs.reference_file,
                    "entryname": inputs.reference_file.basename,
                    "writable": true
                  }
                ]
      }

inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      set -e
      if [[ $0 == *.fasta ]] || [[ $0 == *.fa ]]; then
          echo "Skip extract step"
      elif [[ $0 == *.gz ]]; then
          gunzip -c $0 > "${0%%.*}".fa
          rm $0
      else
          twoBitToFa $0 "${0%%.*}".fa
          rm $0
      fi
      if [ "$#" -ge 1 ]; then
          FILTER=${@:1}
          FILTER=$( IFS=$','; echo "${FILTER[*]}" )
          FILTER=(${FILTER//, / })
          echo "Filtering by" ${FILTER[*]}
          samtools faidx "${0%%.*}".fa ${FILTER[*]} > t.fa
          mv t.fa "${0%%.*}".fa
          rm "${0%%.*}".fa.fai
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run twoBitToFa or gunzip to extract and samtools to filter chromosomes

  reference_file:
    type: File
    inputBinding:
      position: 6
    doc: "Reference genome *.2bit, *.fasta, *.fa, *.fa.gz, *.fasta.gz file"

  chr_list:
    type:
      - "null"
      - string
      - string[]
    inputBinding:
      position: 8
    doc: "List of the chromosomes to be included into the output file. If pass as string, should be comma-separated"

outputs:

  fasta_file:
    type: File
    outputBinding:
      glob: "*"
    doc: "Reference genome FASTA file"


baseCommand: ["bash", "-c"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/ucsc-metadata.yaml

s:name: "ucsc-twobit-to-fa"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/ucsc-twobit-to-fa.cwl
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
  twoBitToFa - Convert all or part of .2bit file to fasta.
  Outputs only those chromosomes that are set in chr_list intput.
  Tool will fail if you include in chr_list those chromosomes that are absent in 2bit file.
  If gz is provided - use gunzip instead of twoBitToFa
  If FASTA file is provided, do nothing

s:about: |
  usage:
    twoBitToFa input.2bit output.fa
  options:
    -seq=name       Restrict this to just one sequence.
    -start=X        Start at given position in sequence (zero-based).
    -end=X          End at given position in sequence (non-inclusive).
    -seqList=file   File containing list of the desired sequence names 
                    in the format seqSpec[:start-end], e.g. chr1 or chr1:0-189
                    where coordinates are half-open zero-based, i.e. [start,end).
    -noMask         Convert sequence to all upper case.
    -bpt=index.bpt  Use bpt index instead of built-in one.
    -bed=input.bed  Grab sequences specified by input.bed. Will exclude introns.
    -bedPos         With -bed, use chrom:start-end as the fasta ID in output.fa.
    -udcDir=/dir/to/cache  Place to put cache for remote bigBed/bigWigs.

  Sequence and range may also be specified as part of the input
  file name using the syntax:
        /path/input.2bit:name
    or
        /path/input.2bit:name
    or
        /path/input.2bit:name:start-end