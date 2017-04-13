#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
        return inputs.bowtie_log.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".stat";
    };

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
      TOTAL, ALIGNED, SUPRESSED, USED = 100, 80, 0, 0
      with open(sys.argv[1], 'r') as bowtie_log:
        for line in bowtie_log:
          if 'processed:' in line:
            TOTAL = int(line.split('processed:')[1])
          if 'alignment:' in line:
            ALIGNED = int(line.split('alignment:')[1].split()[0])
          if 'due to -m:' in line:
            SUPRESSED = int(line.split('due to -m:')[1].split()[0])
      USED = ALIGNED
      with open(sys.argv[2], 'r') as rmdup_log:
        for line in rmdup_log:
          if '/' in line and 'Skip' not in line:
            splt = line.split('/')
            USED = int((splt[1].split('='))[0].strip()) - int((splt[0].split(']'))[1].strip())
      print TOTAL, ALIGNED, SUPRESSED, USED
    inputBinding:
      position: 5
    doc: |
      Python script to get TOTAL, ALIGNED, SUPRESSED, USED

  bowtie_log:
    type: File
    inputBinding:
      position: 6
    doc: |
      Log file from Bowtie

  rmdup_log:
    type: File
    inputBinding:
      position: 7
    doc: |
      Log file from samtools rmdup

outputs:

  output:
    type: File
    outputBinding:
      glob: $(default_output_filename())


baseCommand: [python, '-c']
arguments:
  - valueFrom: $(" > " + default_output_filename())
    position: 100000
    shellQuote: false
