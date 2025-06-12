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
    type: File?
    inputBinding:
      position: 3
    doc: |
      preseq standard output estimates file. If this tool fails the step, it's possibly because the 'estimates_file' was not produced due to preseq lc_extrap:
        ERROR:  max count before zero is les than min required count (4), sample not sufficiently deep or duplicates removed
      therefore this input file is optional
  mapped_reads:
    type: int
    inputBinding:
      position: 4
    doc: |
      mapped read count
outputs:
  estimates_file_plot_data:
    type: File?
    outputBinding:
      glob: estimates_file_plot_data.tsv
    doc: |
      Formatted estimates file for preseq plotting.
      If estimates_file is not provided, this file cannot be produced, therefore it is optional so the entire workflow does not fail.
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bash
- -c
stdout: preseq-plot-data_stdout.log
stderr: preseq-plot-data_stderr.log
doc: |-
  Tool runs custom bash and awk commands to format the preseq estimates
  file in order to plot the total reads count on the estimates curve.
label: bedtools-fragmentcounts
