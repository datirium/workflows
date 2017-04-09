#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2
  dockerFile: >
    $import: ./dockerfiles/scidap-Dockerfile

inputs:

  script:
    type: string?
    default: |
      #!/usr/bin/env python
      import sys, re
      fragments, islands = 0, 0
      with open(sys.argv[1], 'r') as infile:
        for line in infile:
            if re.match('^# d = ', line):
                fragments = int(line.split('d = ')[1])
                continue
            if re.match('^#', line):
                continue
            if line.strip() != "":
                islands = islands + 1
      islands = islands - 1
      print fragments, '\n', islands
    inputBinding:
      position: 5
    doc: |
      Python script to get ISLANDS and FRAGMENTS from MACS2 output

  input_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      Output file from MACS2 peak calling

outputs:

  fragments:
    type: int
    outputBinding:
      loadContents: true
      glob: "island_count.log"
      outputEval: $(parseInt(self[0].contents.split(/\r?\n/)[0]))

  islands:
    type: int
    outputBinding:
      loadContents: true
      glob: "island_count.log"
      outputEval: $(parseInt(self[0].contents.split(/\r?\n/)[1]))

baseCommand: [python, '-c']
arguments:
  - valueFrom: $(" > island_count.log")
    position: 100000
    shellQuote: false
