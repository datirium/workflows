cwlVersion: v1.0
class: Workflow
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  isoforms_file:
    type: File[]
    doc: Isoforms CSV file(s)
outputs:
  genes_file:
    type: File[]
    outputSource: group_isoforms/genes_file
    doc: Gene expression file(s) in TSV format
  common_tss_file:
    type: File[]
    outputSource: group_isoforms/common_tss_file
    doc: Common tss expression file(s) in TSV format
steps:
  group_isoforms:
    run: ../tools/group-isoforms.cwl
    in:
      isoforms_file: isoforms_file
    scatter:
    - isoforms_file
    out:
    - genes_file
    - common_tss_file
doc: |
  Workflow runs group-isoforms.cwl tool using scatter for isoforms_file input. genes_filename and common_tss_filename inputs are ignored.
sd:version: 100
label: group-isoforms-batch
