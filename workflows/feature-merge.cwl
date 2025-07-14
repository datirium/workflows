cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement
sd:upstream:
  rnaseq_sample:
  - rnaseq-se.cwl
  - rnaseq-pe.cwl
  - rnaseq-se-dutp.cwl
  - rnaseq-pe-dutp.cwl
  - rnaseq-se-dutp-mitochondrial.cwl
  - rnaseq-pe-dutp-mitochondrial.cwl
  - trim-rnaseq-pe.cwl
  - trim-rnaseq-se.cwl
  - trim-rnaseq-pe-dutp.cwl
  - trim-rnaseq-pe-smarter-dutp.cwl
  - trim-rnaseq-se-dutp.cwl
inputs:
  alias:
    type: string
    label: Experiment short name/Alias
    sd:preview:
      position: 1
  rpkm_isoforms:
    type:
    - 'null'
    - File[]
    default: null
    format: http://edamontology.org/format_3752
    label: RNA-Seq experiments
    doc: CSV/TSV file with RPKM grouped by isoform
    sd:upstreamSource: rnaseq_sample/rpkm_isoforms
    sd:localLabel: true
  rpkm_genes:
    type:
    - 'null'
    - File[]
    default: null
    format: http://edamontology.org/format_3475
    label: RNA-Seq experiments
    doc: CSV/TSV file with RPKM grouped by gene
    sd:upstreamSource: rnaseq_sample/rpkm_genes
  rpkm_common_tss:
    type:
    - 'null'
    - File[]
    default: null
    format: http://edamontology.org/format_3475
    label: RNA-Seq experiments
    doc: CSV/TSV file with RPKM grouped by common TSS
    sd:upstreamSource: rnaseq_sample/rpkm_common_tss
  sample_names:
    type: string[]
    label: RNA-Seq experiments
    doc: |
      Aliases for RNA-Seq experiments to be used as the names for the reported
      columns in the merged file. Order corresponds to the rpkm_genes input
    sd:upstreamSource: rnaseq_sample/alias
  group_by:
    type:
    - 'null'
    - type: enum
      symbols:
      - isoforms
      - genes
      - common tss
    default: genes
    label: Group by
    doc: 'Grouping method for features: isoforms, genes or common tss'
  merge_by_columns:
    type:
    - 'null'
    - string[]
    default:
    - RefseqId
    - GeneId
    - Chrom
    - TxStart
    - TxEnd
    - Strand
    label: Columns to merge experiments by
    doc: Column names to merge feature files by
    sd:layout:
      advanced: true
outputs:
  merged_file:
    type: File
    label: Merged feature expression file
    format: http://edamontology.org/format_3475
    doc: |
      Merged by merge_by_columns feature expression file
    outputSource: feature_merge/merged_file
    sd:visualPlugins:
    - syncfusiongrid:
        tab: Merged Feature Expression
        Title: Merged Feature Expression
  feature_merge_stdout_log:
    type: File
    format: http://edamontology.org/format_2330
    label: Merging stdout log
    doc: Merging stdout log
    outputSource: feature_merge/stdout_log
  feature_merge_stderr_log:
    type: File
    format: http://edamontology.org/format_2330
    label: Merging stderr log
    doc: Merging stderr log
    outputSource: feature_merge/stderr_log
steps:
  feature_merge:
    run: ../tools/feature-merge.cwl
    in:
      feature_files:
        source:
        - group_by
        - rpkm_isoforms
        - rpkm_genes
        - rpkm_common_tss
        valueFrom: |
          ${
              if (self[0] == "isoforms") {
                return self[1];
              } else if (self[0] == "genes") {
                return self[2];
              } else {
                return self[3];
              }
          }
      feature_aliases: sample_names
      mergeby: merge_by_columns
    out:
    - merged_file
    - stdout_log
    - stderr_log
label: Feature expression merge - combines feature expression from several experiments
doc: |
  Feature expression merge - combines feature expression from several experiments
  =========================================================================

  Workflows merges RPKM (by default) gene expression from several
  experiments based on the values from GeneId, Chrom, TxStart, TxEnd and
  Strand columns (by default). Reported unique columns are renamed based on
  the experiments names.
sd:version: 100
