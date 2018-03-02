#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: EnvVarRequirement
  envDef:
  - envName: AL_USE_CONCATENATED_GENOME
    envValue: $(inputs.genome2?'0':'1')
  - envName: AL_BWA_ALN_PARAMS
    envValue: $(inputs.al_bwa_aln_params)
  - envName: AL_DIR_TOOLS
    envValue: /usr/local/bin/
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/alea:v1.2.2


inputs:

  al_bwa_aln_params:
    type: string?
    default: '-k 0 -n 0 -t 4'

  input_reads:
    type: File[]
    inputBinding:
      position: 2
    doc: |
      input_reads_1 the 1st input reads file in fastq.
      input_reads_2 (paired end) the 2nd input reads file in fastq.
      (fastq.gz or bam is supported when using BWA)

  genome1:
    type: File
    inputBinding:
      position: 3
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .sa
    doc: |
      (when AL_USE_CONCATENATED_GENOME=0)
      path to the indexed reference for 1st insilico genome (of strain1).
      for BWA, specifiy the fasta file.
      for Bowtie, specify index filename prefix (minus trailing .X.ebwt or .X.bt2)

  genome2:
    type: File?
    inputBinding:
      position: 4
    secondaryFiles:
    - .amb
    - .ann
    - .bwt
    - .fai
    - .pac
    - .sa
    doc: |
      (when AL_USE_CONCATENATED_GENOME=0)
      path to the indexed reference for 2nd insilico genome (of strain2).
      for BWA, specifiy the fasta file.
      for Bowtie, specify basename of index files.

  strain1:
    type: string
    inputBinding:
      position: 5
    doc: name of strain1 exactly as specified in the vcf file (e.g. hap1)

  strain2:
    type: string
    inputBinding:
      position: 6
    doc: |
      name of strain2 exactly as specified in the vcf file (e.g. hap2)


outputs:

  all_data:
    type: File[]
    outputBinding:
      glob: "*"


baseCommand: [alea, alignReads]
arguments:
- valueFrom: $(inputs.input_reads.length==1?"-s":"-p")
  position: 1
- valueFrom: $('./')
  position: 7


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/alea-metadata.yaml

s:name: "alea-alignreads"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/alea-alignreads.cwl
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
  Tool runs alea alignReads

s:about: |
  Usage:
    using concatenated genome method (AL_USE_CONCATENATED_GENOME=1):
        alea alignReads <-s/-p> <input_reads_1 [input_reads_2]> <genome_concat> <strain1 strain2> <outputPrefix>

    using separate insilico genomes method (AL_USE_CONCATENATED_GENOME=0):
        alea alignReads <-s/-p> <input_reads_1 [input_reads_2]> <genome1  genome2> <strain1 strain2> <outputPrefix>

  Options:
    -s              to align single-end reads (requires one input file)
    -p              to align paired-end reads (requires two input files)

    input_reads_1   the 1st input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)

    input_reads_2   (paired end) the 2nd input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)

    genome_concat   (when AL_USE_CONCATENATED_GENOME=1)
                    path to the indexed reference for concatenated insilico genome.
                    for BWA, specifiy path to the fasta.
                    for Bowtie, specify basename of index file

    genome1         (when AL_USE_CONCATENATED_GENOME=0)
                    path to the indexed reference for 1st insilico genome (of strain1).
                    for BWA, specifiy the fasta file.
                    for Bowtie, specify index filename prefix (minus trailing .X.ebwt or .X.bt2)

    genome1         (when AL_USE_CONCATENATED_GENOME=0)
                    path to the indexed reference for 2nd insilico genome (of strain2).
                    for BWA, specifiy the fasta file.
                    for Bowtie, specify index filename prefix (minus trailing .X.ebwt or .X.bt2)

    strain1         name of strain1
                    (e.g. hap1 or CASTEiJ)

    strain2         name of strain2
                    (e.g. hap2 or C57BL6J)

    outputPrefix    prefix for output files, including the full path, without an extension
                    (e.g. ./TSC_H3K36me3 )

  Output:

    outputPrefix_strain1_starin2.sam   (when AL_USE_CONCATENATED_GENOME=1)
                                       all reads aligned to the concatenated insilico genome

    outputPrefix_strain1_all.sam       (when AL_USE_CONCATENATED_GENOME=0)
                                       all reads aligned to the first insilico genome

    outputPrefix_strain2_all.sam       (when AL_USE_CONCATENATED_GENOME=0)
                                       all reads aligned to the second insilico genome

    outputPrefix_strain1.bam           allelic reads for strain1 genome (sorted bam)

    outputPrefix_strain2.bam           allelic reads for strain2 genome (sorted bam)

  Examples:
    (AL_USE_CONCATENATED_GENOME=1, AL_USE_BWA=1)
    alea alignReads -s H3K36me3.fastq CASTEiJ_C57BL6J.fasta CASTEiJ C57BL6J ./H3K36me3
    alea alignReads -p H3K36me3_1.fastq H3K36me3_2.fastq CASTEiJ_C57BL6J.fasta CASTEiJ C57BL6J ./H3K36me3

    (AL_USE_CONCATENATED_GENOME=0, AL_USE_BWA=1)
    alea alignReads -s H3K36me3.fastq CASTEiJ.fasta C57BL6J.fasta CASTEiJ C57BL6J ./H3K36me3

    (AL_USE_CONCATENATED_GENOME=0, AL_USE_BOWTIE1=1)
    alea alignReads -s H3K36me3.fastq bowtie1-index/CASTEiJ bowtie1-index/C57BL6J CASTEiJ C57BL6J ./H3K36me3
