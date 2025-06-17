cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.42


inputs:

  gct_file:
    type: File?
    inputBinding:
      position: 1

  tsv_file:
    type: File
    inputBinding:
      position: 2


outputs:

  tag_dnst_htmp_gct:
    type: File?
    outputBinding:
      glob: "*.gct"

  tag_dnst_htmp_html:
    type: File?
    outputBinding:
      glob: "*.html"

  tag_dnst_htmp_tsv:
    type: File?
    outputBinding:
      glob: "*.tsv"


baseCommand: ["sc_utils_extend_htmp.R"]