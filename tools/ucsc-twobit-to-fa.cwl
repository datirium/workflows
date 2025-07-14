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
    doc: Reference genome *.2bit, *.fasta, *.fa, *.fa.gz, *.fasta.gz file
  chr_list:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      position: 8
    doc: List of the chromosomes to be included into the output file. If pass as string, should be comma-separated
outputs:
  fasta_file:
    type: File
    outputBinding:
      glob: '*'
    doc: Reference genome FASTA file
baseCommand:
- bash
- -c
doc: |
  twoBitToFa - Convert all or part of .2bit file to fasta.
  Outputs only those chromosomes that are set in chr_list intput.
  Tool will fail if you include in chr_list those chromosomes that are absent in 2bit file.
  If gz is provided - use gunzip instead of twoBitToFa
  If FASTA file is provided, do nothing
label: ucsc-twobit-to-fa
