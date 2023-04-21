cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-ubuntu22:v1.0.0


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for bash script in bedtools-fragmentcounts.cwl tool:\n\n"
      printf "INPUTS:\n"
      echo "\$0 preseq_stderr_log_file - $0"
      echo "\$1 estimates_file - $1"
      echo "\$2 mapped_reads - $2"
      # get actual distinct reads from preseq's verbose stderr log
      dr=$(grep "DISTINCT READS" $0 | sed 's/.*= //')
      # find index where mapped reads count should be along the x-axis
      yval_index=$(tail -n+2 $1 | awk -F'\t' -v tr=$2 '{if($1<tr){x++}}END{print(x)}')
      # generate new headers for formatted plot file
      head -1 $1 | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\n",$1,"MAPPED_READS",$3,$4,$2)}' > estimates_file_plot_data.tsv
      # need to switch col2 and 5, first column listed with be the "top" line
      tail -n+2 "$1" | awk -F'\t' -v yval_index=$yval_index -v yval=$dr '{if(NR==yval_index){printf("%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n",$1,yval,$3,$4,$2)}else{printf("%.0f\t%.0f\t%.0f\t%.0f\t%s\n",$1,"null",$3,$4,$2)}}' >> estimates_file_plot_data.tsv
    inputBinding:
        position: 1

  preseq_stderr_log_file:
    type: File
    inputBinding:
      position: 2
    doc: |
      preseq stderr verbose log file containing actual distinct read count

  estimates_file:
    type: File
    inputBinding:
      position: 3
    doc: |
      preseq standard output estimates file

  mapped_reads:
    type: int
    inputBinding:
      position: 4
    doc: |
      mapped read count


outputs:

  estimates_file_plot_data:
    type: File
    outputBinding:
      glob: estimates_file_plot_data.tsv
    doc: |
      formatted estimates file for preseq plotting

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]
stdout: preseq-plot-data_stdout.log
stderr: preseq-plot-data_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "bedtools-fragmentcounts"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/preseq-plot-data.cwl
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
  Tool runs custom bash and awk commands to format the preseq estimates
  file in order to plot the total reads count on the estimates curve.