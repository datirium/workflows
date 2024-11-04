cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 244130                                     # equal to ~264GB
    coresMin: 20


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kraken2:v1.0.0


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nStdout log file for k2-classify-pe.cwl tool:\n"
      DATABASE=$0; R1=$1; R2=$2
      printf "INPUTS:\n"
      printf "\$0 - $DATABASE\n"
      printf "\$1 - $R1\n"
      printf "\$2 - $R2\n"
      printf "EXECUTION:\n"
      #   commands start
      printf "\trun classification for PE reads\n"
      kraken2 --db $DATABASE --threads 20 --paired --classified-out k2_classified_reads#.fastq --unclassified-out k2_unclassified_reads#.fastq --output k2.output --report k2.report $R1 $R2 2> k2.stderr
      printf "\tformatting outputs\n"
      # make stderr output markdown compatible for overview tab view
      head -1 k2.stderr > parsed.stderr
      tail -n+2 k2.stderr | sed 's/^ *//' | awk '{printf(" - %s\n",$0)}' >> parsed.stderr
      # format report into tsv for table tab view
      printf "percent_classified\treads_assigned_at_and_below_taxid\treads_assigned_directly_to_taxid\ttaxonomic_rank\ttaxid\tname\n" > k2_report.tsv
      sed 's/^ *//' k2.report | sed 's/\t  */\t/' >> k2_report.tsv
      printf "\tgenerate krona plot\n"
      python3 /usr/local/src/KrakenTools/kreport2krona.py -r k2.report -o k2.krona --no-intermediate-ranks
      perl /usr/local/src/Krona/KronaTools/scripts/ImportText.pl -o krona.html k2.krona
      printf "END OF SCRIPT\n"
    inputBinding:
        position: 1

  k2db:
    type: Directory
    label: "Kraken2 database for taxonomic classification"
    inputBinding:
      position: 5
    doc: "Name of kraken2 database to use for kraken2 classify."

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

  k2_classified_R1:
    type: File
    outputBinding:
      glob: "k2_classified_reads_1.fastq"

  k2_classified_R2:
    type: File
    outputBinding:
      glob: "k2_classified_reads_2.fastq"

  k2_unclassified_R1:
    type: File
    outputBinding:
      glob: "k2_unclassified_reads_1.fastq"

  k2_unclassified_R2:
    type: File
    outputBinding:
      glob: "k2_unclassified_reads_2.fastq"

  k2_output:
    type: File
    outputBinding:
      glob: "k2.output"

  k2_report:
    type: File
    outputBinding:
      glob: "k2.report"

  k2_report_tsv:
    type: File
    outputBinding:
      glob: "k2_report.tsv"

  k2_stderr:
    type: File
    outputBinding:
      glob: "parsed.stderr"

  krona_html:
    type: File
    outputBinding:
      glob: "krona.html"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]
stdout: k2-stdout.log
stderr: k2-stderr.log





$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "k2-classify-pe"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/k2-classify-pe.cwl
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
    Tool downloads specified kraken2 database from https://benlangmead.github.io/aws-indexes/k2.
    Resulting directory is used as upstream input to kraken2 classify tools.
