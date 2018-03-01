#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: EnvVarRequirement
  envDef:
  - envName: AL_USE_CONCATENATED_GENOME
    envValue: $(inputs.concat_genome?"1":"0")
  - envName: AL_BWA_ALN_PARAMS
    envValue: -k 0 -n 0 -t 4
  - envName: AL_DIR_TOOLS
    envValue: /usr/local/bin/
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/alea:v1.2.2


inputs:

  concat_genome:
    type: boolean
    default: false

  reference:
    type: File
    inputBinding:
      position: 2
    secondaryFiles:
    - .fai
    doc: |
      the reference genome fasta file

  phased:
    type: File
    inputBinding:
      position: 3
    secondaryFiles:
    - .tbi
    doc: |
      the phased variants vcf file (including SNPs and Indels)
      or the phased SNPs (should be specified first)

  phased_indels:
    type: File?
    inputBinding:
      position: 4
    secondaryFiles:
    - .tbi
    doc: |
      the phased Indels (should be specified second)

  strain1:
    type: string
    inputBinding:
      position: 5
    doc: |
      name of strain1 exactly as specified in the vcf file (e.g. hap1)

  strain2:
    type: string
    inputBinding:
      position: 6
    doc: |
      name of strain2 exactly as specified in the vcf file (e.g. hap2)


outputs:

  strain12_indices:
    type: File?
    outputBinding:
      glob: $(inputs.strain1 + "_" + inputs.strain2 + ".fasta")
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .sa

  strain1_indices:
    type: File?
    outputBinding:
      glob: $(inputs.strain1 + ".fasta")
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .refmap
    - .sa

  strain2_indices:
    type: File?
    outputBinding:
      glob: $(inputs.strain2 + ".fasta")
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .refmap
    - .sa

baseCommand: [alea, createGenome]
arguments:
- valueFrom: $(inputs.phased_indels?"-snps-indels-separately":[])
  position: 1
- valueFrom: $("./")
  position: 7

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/alea-metadata.yaml

s:downloadUrl: https://github.com/common-workflow-language/workflows/blob/master/tools/alea-createGenome.cwl
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
