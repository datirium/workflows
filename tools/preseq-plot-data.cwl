cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-ubuntu22:v1.0.0


inputs:

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for bash script in bedtools-fragmentcounts.cwl tool:\n\n"
      printf "INPUTS:\n"
      echo "\$0 trim_fastq_report_file - $0"
      echo "\$1 estimates_file - $1"
      # get total reads from fastx stats file
      tr=$(grep "sequences processed in total" "$0" | sed 's/ .*//')
      # find index where total reads should be placed along x-axis
      yval_index=$(tail -n+2 "$1" | awk -F'\t' -v tr=$tr '{if($1<tr){x++}}END{print(x)}')
      # store that y-value
      yval=$(tail -n+2 "$1" | head -$yval_index | tail -1 | cut -f2)
      # generate new headers for formatted plot file
      head -1 "$1" | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\n",$1,"ACTUAL_READ_COUNT",$3,$4,$2)}' > estimates_file_plot_data.tsv
      # need to switch col2 and 5, first column listed with be the "top" line
      tail -n+2 "$1" | awk -F'\t' -v yval_index=$yval_index -v yval=$yval '{if(NR==yval_index){printf("%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n",$1,yval,$3,$4,$2)}else{printf("%.0f\t%.0f\t%.0f\t%.0f\t%s\n",$1,"null",$3,$4,$2)}}' >> estimates_file_plot_data.tsv
    inputBinding:
        position: 1

  trim_fastq_report_file:
    type: File
    inputBinding:
      position: 2
    doc: |
      trimming report file to get total read counts from

  estimates_file:
    type: File
    inputBinding:
      position: 3
    doc: |
      preseq standard output estimates file


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