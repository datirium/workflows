cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.input_file,
                  "entryname": inputs.input_file.basename,
                  "writable": true
                }
              ]
    }
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2
inputs:
  input_file:
    type:
    - File
    inputBinding:
      position: 1
    doc: |
      File to be compressed
  fast:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 2
    doc: |
      Set block size to 100k
  best:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 3
    doc: |
      Set block size to 900k
outputs:
  output_file:
    type:
    - File
    outputBinding:
      glob: |
        ${
          return inputs.input_file.basename + '.bz2';
        }
baseCommand:
- bzip2
doc: |
  Tool compresses `input_file` to `*.bz2`.
  Output file has the same basename, as input file, but with updated `.bz2` extension. `bzip2` exports compressed
  output file alognside the input file. To prevent tool from failing, `input_file` should be staged into output
  directory using `"writable": true`. Setting `writable: true` makes cwl-runner to make a copy of input file and
  mount it to docker container with `rw` mode as part of `--workdir` (if set to false, the file staged into output
  directory will be mounted to docker container separately with `ro` mode)
label: bzip2-compress
