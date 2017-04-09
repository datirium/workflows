#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement
  expressionLib:
  - var new_out_filename = function() {
      if (inputs.out_filename){
        return inputs.out_filename;
      }
      if (inputs.input.location.split('/').slice(-1)[0].split('.').slice(-1)[0] == 'cram'){
        return inputs.input.location.split('/').slice(-1)[0]+'.crai';
      } else if (inputs.csi && !inputs.bai){
        return inputs.input.location.split('/').slice(-1)[0]+'.csi';
      } else {
        return inputs.input.location.split('/').slice(-1)[0]+'.bai';
      }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
  dockerFile: >
    $import: ./dockerfiles/samtools-Dockerfile

inputs:
  input:
    type: File
    inputBinding:
      position: 8
    doc: |
      Input bam file.

  interval:
    type: int?
    inputBinding:
      position: 6
      prefix: -m
    doc: |
      Set minimum interval size for CSI indices to 2^INT [14]

  threads:
    type: int?
    inputBinding:
      position: 7
      prefix: -@
    doc: |
      Sets the number of threads

  out_filename:
    type: string?
    doc: |
      Output filename

  csi:
    type: boolean?
    doc: |
      Generate CSI-format index for BAM files

  bai:
    type: boolean?
    doc: |
      Generate BAI-format index for BAM files [default]

outputs:
  index:
    type: File
    outputBinding:
      glob: $(new_out_filename())

    doc: The index file

baseCommand: [samtools, index]
arguments:
- valueFrom: $(inputs.bai?'-b':inputs.csi?'-c':[])
  position: 5
- valueFrom: $(new_out_filename())
  position: 9

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/samtools-index.cwl
s:codeRepository: https://github.com/common-workflow-language/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: Andrey Kartashov
  s:email: mailto:Andrey.Kartashov@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab
doc: |
  samtools-index.cwl is developed for CWL consortium
