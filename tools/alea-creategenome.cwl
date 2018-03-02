#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: EnvVarRequirement
  envDef:
  - envName: AL_USE_CONCATENATED_GENOME
    envValue: $(inputs.concat_genome?'1':'0')
  - envName: AL_DIR_TOOLS
    envValue: /usr/local/bin/
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/alea:v1.2.2


inputs:

  concat_genome:
    type: boolean?
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
    type:
      - "null"
      - File[]
    outputBinding:
      glob: $('*' + inputs.strain1 + '_' + inputs.strain2 + '*')

  strain1_indices:
    type:
      - "null"
      - File[]
    outputBinding:
      glob: $('*' + inputs.strain1 + '*')

  strain2_indices:
    type:
      - "null"
      - File[]
    outputBinding:
      glob: $('*' + inputs.strain2 + '*')


baseCommand: [alea, createGenome]
arguments:
- valueFrom: $(inputs.phased_indels?'-snps-indels-separately':'')
  position: 1
- valueFrom: $('./')
  position: 7



$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/alea-metadata.yaml

s:name: "alea-creategenome"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/alea-creategenome.cwl
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
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool runs alea createGenome

s:about: |
  Usage:
    alea createGenome reference.fasta phased.vcf.gz strain1 strain2 outputDir
    alea createGenome -snps-indels-separately reference.fasta phased_snps.vcf.gz phased_indels.vcf.gz strain1 strain2 outputDir
  Options:
    reference.fasta        the reference genome fasta file
    phased.vcf.gz          the phased variants vcf file (including SNPs and Indels)
    strain1                name of strain1 exactly as specified in the vcf file (e.g. hap1)
    strain2                name of strain2 exactly as specified in the vcf file (e.g. hap2)
    outputDir              location of the output directory
    -snps-indels-separately    use if SNPs and Indels are in two separate vcf files
    phased-snps.vcf.gz         the phased SNPs (should be specified first)
    phased-indels.vcf.gz       the phased Indels  (should be specified second)
  Output:
    Creates two parental insilico genomes strain1.fasta and strain2.fasta as well
    as alignment indeces.
    A concatenated genome strain1_strain2.fasta will be created if
    AL_USE_CONCATENATED_GENOME=1 is set in alea.config
  Note:
    It is possible to have SNPs and Indels in two separate vcf files. In that case
    use -snps-indels-separately option, and make sure you specify SNPs before Indels.


