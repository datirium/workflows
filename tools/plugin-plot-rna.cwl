cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/plugin-plot-rna:v0.0.4
inputs:
  annotation_file:
    type: File
    inputBinding:
      position: 5
      prefix: --annotation
    doc: |
      Path to the annotation TSV/CSV file
  bambai_pair:
    type: File
    inputBinding:
      position: 6
      prefix: --bam
    secondaryFiles:
    - .bai
    doc: |
      Path to the indexed BAM file
  isoforms_file:
    type: File
    inputBinding:
      position: 7
      prefix: --isoforms
    doc: |
      Path to the isoforms TSV/CSV file
  mapped_reads_number:
    type: int
    inputBinding:
      position: 8
      prefix: --mapped
    doc: |
      Mapped reads number
  output_prefix:
    type: string?
    inputBinding:
      position: 9
      prefix: --output
    doc: |
      Output prefix. Default: ./coverage
  pair:
    type: boolean?
    inputBinding:
      position: 10
      prefix: --pair
    doc: |
      Run as paired end. Default: false
  strand_specificity:
    type:
    - 'null'
    - type: enum
      symbols:
      - 'yes'
      - 'no'
      - reverse
    inputBinding:
      position: 11
      prefix: --stranded
    doc: |
      Whether the data is from a strand-specific assay.
      --stranded no      - a read is considered overlapping with a feature regardless of whether
                           it is mapped to the same or the opposite strand as the feature.
      --stranded yes     - the read has to be mapped to the same strand as the feature.
      --stranded reverse - the read has to be mapped to the opposite strand than the feature.
  minimum_rpkm:
    type: float?
    inputBinding:
      position: 12
      prefix: --minrpkm
    doc: |
      Ignore isoforms with RPKM smaller than --minrpkm.
      Default: 10
  minimum_isoform_length:
    type: float?
    inputBinding:
      position: 13
      prefix: --minlength
    doc: |
      Ignore isoforms shorter than --minlength.
      Default: 1000
  threads:
    type: int?
    inputBinding:
      position: 14
      prefix: --threads
    doc: |
      Threads. Default: 1
outputs:
  gene_body_report_file:
    type: File?
    outputBinding:
      glob: '*gene_body_report.tsv'
  gene_body_plot_png:
    type: File?
    outputBinding:
      glob: '*gene_body_plot.png'
  gene_body_plot_pdf:
    type: File?
    outputBinding:
      glob: '*gene_body_plot.pdf'
  rpkm_distribution_plot_png:
    type: File?
    outputBinding:
      glob: '*rpkm_distribution_plot.png'
  rpkm_distribution_plot_pdf:
    type: File?
    outputBinding:
      glob: '*rpkm_distribution_plot.pdf'
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- plot_rna.R
successCodes:
- 1
stderr: gene_body_stderr.log
stdout: gene_body_stdout.log
label: Gene body average tag density plot and RPKM distribution histogram
doc: |
  Runs R script to produce gene body average tag density plot and RPKM distribution histogram
  Doesn't fail even when we couldn't produce any plots
