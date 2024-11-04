cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 30510
    coresMin: 8

hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-metaphlan:v1.0.0


inputs:

  script_command:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nStdout log file for wgs-metaphlan-pe.cwl tool:\n"
      R1=$0
      R2=$1
      printf "EXECUTION:\n"
      #   commands start
      printf "\trun metaphlan analysis for PE reads\n"
      echo "concat paired reads together (metaphlan does not take paired information into account)"
      cp $R1 combined.fq
      cat $R2 >> combined.fq
      metaphlan combined.fq --bowtie2out metagenome.bowtie2.bz2 --nproc 8 --input_type fastq -o profiled_metagenome.txt
      # clean header for scidap output tab
      tail -n+5 profiled_metagenome.txt > profiled_metagenome_scidap.txt
      echo "need to make a merged table for heatmap - just a single sample so not super interesting"
      merge_metaphlan_tables.py profiled_metagenome.txt > merged_abundance_table.txt
      echo "generate species-only abundance table"
      grep -E "s__|clade" merged_abundance_table.txt \
        | grep -v "t__" \
        | sed "s/^.*|//g" \
          > merged_abundance_table_species.txt
      echo "cleaning up tmp files"
      rm combined.fq
      printf "END OF SCRIPT\n"
    inputBinding:
        position: 1

  read1file:
    type: File
    label: "R1 fastq"
    inputBinding:
      position: 6
    doc: "FASTQ file 1 of paired end read data."

  read2file:
    type: File
    label: "R2 fastq"
    inputBinding:
      position: 7
    doc: "FASTQ file 2 of paired end read data."


outputs:

  classification_alignments_bowtie2:
    type: File
    outputBinding:
      glob: "metagenome.bowtie2.bz2"

  abundance_profile:
    type: File
    outputBinding:
      glob: "profiled_metagenome.txt"

  abundance_profile_scidap:
    type: File
    outputBinding:
      glob: "profiled_metagenome_scidap.txt"

  abundance_table:
    type: File
    outputBinding:
      glob: "merged_abundance_table.txt"

  abundance_table_species:
    type: File
    outputBinding:
      glob: "merged_abundance_table_species.txt"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]
stdout: log.stdout
stderr: log.stderr


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "metaphlan4"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/wgs-metaphlan-pe.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

doc: |
    Tool runs metaphlan4 for metagenomic classification of metagenomic data sets.

    MetaPhlAn runs reads against the "latest" database:
        The database contains genomes from about 1 million microbial genomes, including 236,600 microbial isolates
        and 771,500 metagenomic assembled genomes. The database is used by the MetaPhlAn 4 computational tool to
        profile the composition of microbial communities.



