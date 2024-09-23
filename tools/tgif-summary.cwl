cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-tgif:v1.0.0


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for bash script in tgif-summary.cwl tool:\n\n"
      # inputs
      log="$0"
      all="$1"
      filtered="$2"
      # get metrics
      reads_total=$(grep "mapping.*reads to*" $log | cut -d$' ' -f2)
      reads_uniquely_mapped=$(grep "mapped.*reads uniquely" $log | cut -d$' ' -f2)
      alignment_rate=$(grep "alignment rate" $log | sed 's/.*= //')
      mapped_to_r=$(grep "mapped.* to (-r)" $log | cut -d$' ' -f2)
      mapped_to_i=$(grep "mapped.* to (-i)" $log | cut -d$' ' -f2)
      reads_mapq30=$(grep ".*reads with MAPQ" $log | cut -d$' ' -f1 | sed 's/^\t\+//')
      single_best_aln=$(grep "found single best alignments " $log | cut -d$' ' -f6)
      insertions_all=$(wc -l <$all)
      insertions_filtered=$(tail -n+2 $filtered | wc -l)
      # format report tsv
      printf "%s\t%s\n" "TgIF Results" "Value" > reportsummary.tsv
      printf "Total read count from input FASTQ\t$reads_total\n" >> reportsummary.tsv
      printf "Uniquely mapped reads to genome+vector FASTA\t$reads_uniquely_mapped\n" >> reportsummary.tsv
      printf "%s\t%s\n" "Alignment rate" "$alignment_rate" >> reportsummary.tsv
      printf "Alignments to reference genome\t$mapped_to_r\n" >> reportsummary.tsv
      printf "Alignments to vector sequence\t$mapped_to_i\n" >> reportsummary.tsv
      printf "Alignments to both with MAPQ>=30\t$reads_mapq30\n" >> reportsummary.tsv
      printf "Single best alignments for MAPQ>=30\t$single_best_aln\n" >> reportsummary.tsv
      printf "Total insertion sites founds\t$insertions_all\n" >> reportsummary.tsv
      printf "Total probable insertion sites founds\t$insertions_filtered\n" >> reportsummary.tsv
      # generate formatted markdown table
      Rscript /usr/local/bin/run_tsv_to_kable.R reportsummary.tsv reportsummary.md

    inputBinding:
        position: 1

  tgif_ncats_log_file:
    type: File
    inputBinding:
      position: 2
    doc: "log file from tgif-ncats.sh script"

  insertions_all:
    type: File
    inputBinding:
      position: 3
    doc: "file with all insertion sites from tgif-ncats.sh script"

  insertions_filtered:
    type: File
    inputBinding:
      position: 4
    doc: "file with filtered insertion sites from tgif-ncats.sh script"


outputs:

  summary_file:
    type: File?
    outputBinding:
      glob: "reportsummary.md"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]

stdout: tgif-summary_stdout.log
stderr: tgif-summary_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "tgif-summary"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/tgif-summary.cwl
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
  Tool runs simple bash commands to format a report tsv file for translation to md with kableExtra in R.