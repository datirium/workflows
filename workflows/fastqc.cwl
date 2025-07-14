cwlVersion: v1.0
class: Workflow
requirements:
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_root = function(basename) { return basename.split('.').slice(0,1).join('.'); };
inputs:
  alias:
    type: string
    label: Experiment short name/Alias
    sd:preview:
      position: 1
  fastq_file:
    type:
    - File
    - type: array
      items: File
    label: FASTQ file
    format: http://edamontology.org/format_1930
    doc: Reads data in a FASTQ format, optionally compressed
outputs:
  fastqc_report:
    type: File
    outputSource: rename_fastqc_report/target_file
    label: FastqQC HTML report
    doc: FastqQC HTML report
    sd:visualPlugins:
    - linkList:
        tab: Overview
        target: _blank
steps:
  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out:
    - fastq_file
  run_fastqc:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_fastq/fastq_file
    out:
    - html_file
  rename_fastqc_report:
    run: ../tools/rename.cwl
    in:
      source_file: run_fastqc/html_file
      target_filename:
        source: fastq_file
        valueFrom: $(get_root(self.basename)+"_fastqc_report.html")
    out:
    - target_file
label: FastQC - a quality control tool for high throughput sequence data
doc: |-
  FastQC - a quality control tool for high throughput sequence data
  =====================================

  FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming
  from high throughput sequencing pipelines. It provides a modular set of analyses which you can use
  to give a quick impression of whether your data has any problems of which you should be aware before
  doing any further analysis.

  The main functions of FastQC are:

  - Import of data from FastQ files (any variant)
  - Providing a quick overview to tell you in which areas there may be problems
  - Summary graphs and tables to quickly assess your data
  - Export of results to an HTML based permanent report
  - Offline operation to allow automated generation of reports without running the interactive application
sd:version: 100
