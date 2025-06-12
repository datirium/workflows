cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: scrowley1/scidap-deseq:v0.1.0
inputs:
  isoforms_file:
    type: File
    inputBinding:
      position: 5
      prefix: --isoforms
    doc: Isoforms CSV file
  genes_filename:
    type: string?
    inputBinding:
      position: 6
      prefix: --gene
    doc: Output TSV gene expression filename
  common_tss_filename:
    type: string?
    inputBinding:
      position: 7
      prefix: --tss
    doc: Output TSV common tss expression filename
outputs:
  genes_file:
    type: File
    outputBinding:
      glob: $(inputs.genes_filename?inputs.genes_filename:"*genes.tsv")
    doc: Output TSV gene expression file
  common_tss_file:
    type: File
    outputBinding:
      glob: $(inputs.common_tss_file?inputs.common_tss_file:"*common_tss.tsv")
    doc: Output TSV common tss expression file
baseCommand:
- get_gene_n_tss.R
doc: |
  Tool runs get_gene_n_tss.R script to group isoforms by gene and common TSS
label: group-isoforms
