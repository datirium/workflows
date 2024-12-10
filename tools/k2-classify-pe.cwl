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
      exec 1> error_msg.txt 2>&1
      printf "k2-classify-pe.cwl\n$(date)\n"
      DATABASE=$0; R1=$1; R2=$2
      printf "INPUTS:\n"
      printf "\$0 - $DATABASE\n"
      printf "\$1 - $R1\n"
      printf "\$2 - $R2\n\n"
      #   commands start
      printf "\trun classification for PE reads\n"
      kraken2 --db $DATABASE --threads 20 --paired --classified-out k2_classified_reads#.fastq --unclassified-out k2_unclassified_reads#.fastq --output k2.output --report k2.report $R1 $R2 2> k2.stderr
      # check that k2 output is not empty
      if [[ $(awk 'END{print(NR)}' k2.report) == "0" ]]; then
        echo "k2.report file is empty; exiting" >> error_report.txt
        exit 1
      fi
      printf "\tformatting outputs\n"
      # make stderr output markdown compatible for overview tab view
      head -1 k2.stderr > parsed.stderr
      tail -n+2 k2.stderr | sed 's/^ *//' | awk '{printf(" - %s\n",$0)}' >> parsed.stderr
      # format report into tsv for table tab view
      printf "percent_classified\treads_assigned_at_and_below_taxid\treads_assigned_directly_to_taxid\ttaxonomic_rank\ttaxid\tname\n" > k2_report.tsv
      sed 's/^ *//' k2.report | sed 's/\t  */\t/' >> k2_report.tsv
      printf "\tgenerate krona plot\n"
      python3 /usr/local/src/KrakenTools/kreport2krona.py -r k2.report -o k2.krona --no-intermediate-ranks
      # check that kreport2krona didn't produce no or empty output
      if [[ ! -f k2.krona ]]; then
        echo "k2.krona file was not produced; exiting" >> error_report.txt
        exit 1
      fi
      if [[ $(awk 'END{print(NR)}' k2.krona) == "0" ]]; then
        echo "k2.krona file is empty; exiting" >> error_report.txt
        exit 1
      fi
      perl /usr/local/src/Krona/KronaTools/scripts/ImportText.pl -o krona.html k2.krona
      # check that transform to html didn't produce no or empty output
      if [[ ! -f krona.html ]]; then
        echo "krona.html file was not produced; exiting" >> error_report.txt
        exit 1
      fi
      if [[ $(awk 'END{print(NR)}' krona.html) == "0" ]]; then
        echo "krona.html file is empty; exiting" >> error_report.txt
        exit 1
      fi
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

  error_msg:
    type: File?
    outputBinding:
      glob: "error_msg.txt"

  error_report:
    type: File?
    outputBinding:
      glob: "error_report.txt"

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


baseCommand: ["bash", "-c"]


label: "k2-classify-pe"
doc: |
    Tool downloads specified kraken2 database from https://benlangmead.github.io/aws-indexes/k2.
    Resulting directory is used as upstream input to kraken2 classify tools.
