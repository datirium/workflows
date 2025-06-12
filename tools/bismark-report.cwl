cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2
inputs:
  alignment_report:
    type: File
    inputBinding:
      prefix: --alignment_report
    label: Bismark alignment and methylation report
    doc: Bismark generated alignment and methylation summary report
  splitting_report:
    type: File?
    inputBinding:
      prefix: --splitting_report
    label: Methylation extraction log
    doc: Log file giving summary statistics about methylation extraction
  mbias_report:
    type: File?
    inputBinding:
      prefix: --mbias_report
    label: Methylation bias plot
    doc: QC data showing methylation bias across read lengths
  deduplication_report:
    type: File?
    inputBinding:
      prefix: --dedup_report
    label: Bismark deduplication report
    doc: Bismark generated deduplication report
  nucleotide_report:
    type: File?
    inputBinding:
      prefix: --nucleotide_report
    label: Nucleotide coverage report
    doc: Bismark nucleotide coverage report
outputs:
  collected_report:
    type: File
    label: HTML report page
    doc: Bismark generated graphical HTML report page
    outputBinding:
      glob: '*'
baseCommand:
- bismark2report
doc: |
  This tool uses a Bismark alignment report to generate a graphical HTML report page. Optionally, further reports of
  the Bismark suite such as deduplication, methylation extractor splitting or M-bias reports can be specified as well.

  Skipped arguments:
    -o/--output
    --dir
label: bismark-report
