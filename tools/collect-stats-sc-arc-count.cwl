cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_prefix = function() { if (inputs.output_prefix) { return inputs.output_prefix; } var root = inputs.metrics_summary_report.basename.split('.').slice(0,-1).join('.'); var suffix = "_stats"; return (root == "")?inputs.metrics_summary_report.basename+suffix:root+suffix; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scstats:v0.0.2
inputs:
  metrics_summary_report:
    type: File
    inputBinding:
      position: 6
      prefix: --metrics
  output_prefix:
    type: string?
    inputBinding:
      position: 7
      prefix: --output
      valueFrom: $(get_output_prefix())
    default: ''
outputs:
  collected_statistics_yaml:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".yaml")
  collected_statistics_tsv:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".tsv")
  collected_statistics_md:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".md")
baseCommand:
- cell_ranger_arc_count_stats.py
label: Cell Ranger ARC Count Statistics
doc: |
  Cell Ranger ARC Count Statistics
  ================================

  Collects statistics from Cell Ranger ARC Count experiment
