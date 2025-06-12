cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ResourceRequirement
  ramMin: 7620
  coresMin: 1
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.fasta_file,
                  "entryname": inputs.fasta_file.basename,
                  "writable": true
                }
              ]
    }
inputs:
  fasta_file:
    type: File
    inputBinding:
      position: 5
    doc: Genome FASTA file
outputs:
  fai_file:
    type: File
    outputBinding:
      glob: '*.fai'
    doc: FAI index file
baseCommand:
- samtools
- faidx
doc: |
  Generates FAI index file for input FASTA file
  Output file has the same basename, as input file, but with updated `.fai` extension. `samtools faidx` exports
  output file alognside the input file. To prevent tool from failing, `input_file` should be staged into output
  directory using `"writable": true`. Setting `writable: true` makes cwl-runner to make a copy of input file and
  mount it to docker container with `rw` mode as part of `--workdir` (if set to false, the file staged into output
  directory will be mounted to docker container separately with `ro` mode)
label: samtools-faidx
